#' Generate Knockoff Data for SNPs
#'
#' This function takes a genotype matrix and generates knockoff data using
#' second-order knockoffs with semidefinite programming (SDP) method.
#'
#' @param genetic_variants Matrix of genetic variants (SNPs) where rows are
#'   samples and columns are SNPs.
#'
#' @return Knockoff genotype matrix with the same dimensions as the original
#'   genotype matrix (n_samples x n_SNPs).
#'
#' @export
generate_knockoff_data <- function(genetic_variants) {
  # Extract only the SNP columns (those that start with "chr") to create the genotype matrix
  G <- as.matrix(genetic_variants)
  
  # Basic operations for knockoff preparation
  sd.G <- apply(G, 2, sd)
  cor.G <- matrix(as.numeric(corpcor::cor.shrink(G, verbose = FALSE)), nrow = ncol(G))
  
  # Generate preliminary fit for knockoff creation
  # GhostKnockoff package is in Depends, so functions are available
  fit.prelim.M.sdp <- GhostKnockoff::GhostKnockoff.prelim(cor.G, M = 1, method = 'sdp')
  
  # Second order knockoffs - sdp
  mu.G <- colMeans(G)
  scale.G <- t((t(G) - mu.G) / sd.G)  # Center and scale each SNP column
  # Generate the knockoff matrix
  M <- 1
  E <- matrix(rnorm(nrow(scale.G) * M * ncol(scale.G)), nrow(scale.G), M * ncol(scale.G))
  # Convert V.left to a standard matrix as it’s matrix-like from the Matrix package
  V.left <- as.matrix(fit.prelim.M.sdp$V.left)
  scale.G_k1.M <- t(apply(scale.G %*% t(fit.prelim.M.sdp$P.each), 1, rep, times = M)) + E %*% t(V.left)
  
  # Transform back to original scale
  temp.G_k1.M <- t(rep(sd.G, M) * t(scale.G_k1.M) + rep(mu.G, M))
  G_k <- as.matrix(temp.G_k1.M)  # Knockoff matrix

  message("Genotype knockoff dimensions: ", dim(G_k)[1], " samples x ", dim(G_k)[2], " SNPs")

  return(G_k)
}


#' Prepare Interaction Terms
#'
#' Creates interaction terms between genetic variants and heterogeneity variables.
#'
#' @param genetic_variants Matrix of original or knockoff genetic variants
#'   (n_samples x n_variants).
#' @param Z Vector or matrix representing heterogeneity variable(s)
#'   (e.g., EUR or 4 PCs) (n_samples x n_Z_variables).
#'
#' @return Matrix of interaction terms (n_samples x (n_variants * n_Z_variables)).
#'
#' @keywords internal
prepare_interaction_terms <- function(genetic_variants, Z) { 
  # Ensure Z is a matrix for consistent processing
  Z <- as.matrix(Z)

  # Initialize an empty list to store interaction terms
  interaction_terms <- list()
  # Populate interaction_terms by applying each column of Z to genetic_variants
  for (k in seq_len(ncol(Z))) {
    interaction_terms[[k]] <- sweep(genetic_variants, 1, Z[, k], "*")
  }

  # Combine all interaction terms into a single matrix
  interaction_terms <- do.call(cbind, interaction_terms)

  return(interaction_terms)
}


#' Combine Predictors for Lasso Model
#'
#' Combines all predictor matrices into a single design matrix for model fitting.
#'
#' @param unpenalized_covariates Matrix of covariates (not penalized)
#'   (n_samples x n_covariates).
#' @param genetic_variants Original genetic variants matrix
#'   (n_samples x n_variants).
#' @param interaction_terms Interaction terms for genetic variants
#'   (n_samples x n_interaction_terms).
#' @param genetic_variants_knockoff Knockoff genetic variants matrix
#'   (n_samples x n_variants).
#' @param interaction_terms_knockoff Interaction terms for knockoff genetic variants
#'   (n_samples x n_interaction_terms).
#' 
#' @return A combined predictor matrix for model fitting
#'   (n_samples x total_predictors).
#'
#' @keywords internal
combine_predictors <- function(unpenalized_covariates, genetic_variants, interaction_terms, genetic_variants_knockoff, interaction_terms_knockoff) {
  message("--- Combining all variables and generating X_matrix ...")
  message("unpenalized_covariates dimensions: ", dim(unpenalized_covariates)[1], " samples x ", dim(unpenalized_covariates)[2], " unpenalized covariates")
  message("genetic_variants dimensions: ", dim(genetic_variants)[1], " samples x ", dim(genetic_variants)[2], " SNPs")
  message("interaction_terms dimensions: ", dim(interaction_terms)[1], " samples x ", dim(interaction_terms)[2], " terms")
  message("genetic_variants_knockoff dimensions: ", dim(genetic_variants_knockoff)[1], " samples x ", dim(genetic_variants_knockoff)[2], " SNPs knockoff")
  message("interaction_terms_knockoff dimensions: ", dim(interaction_terms_knockoff)[1], " samples x ", dim(interaction_terms_knockoff)[2], " terms knockoff")
  
  predictors <- cbind(unpenalized_covariates, genetic_variants, interaction_terms, genetic_variants_knockoff, interaction_terms_knockoff)
  message("Combined predictors dimensions: ", dim(predictors)[1], " samples x ", dim(predictors)[2], " all variables")
  return(as.matrix(predictors))
}

#' Prepare Groups Information for Group Lasso Model
#'
#' Prepares group structure information for group lasso modeling.
#'
#' @param num_snps Number of SNPs (equal to ncol(genetic_variants)).
#' @param num_z Number of interacting variables (number of PCs).
#' @param num_u Number of unpenalized covariates (sum of additional vars and Z).
#' @param num_all Total number of predictors (first dimension of matrix X).
#' 
#' @return List with components:
#'   \item{new_order}{Vector of new column order for predictors.}
#'   \item{recovery_index}{Index to recover original coefficient order.}
#'   \item{grp_for_grplasso}{Group indices for group lasso fitting.}
#'
#' @keywords internal
prepare_feat_grps <- function(num_snps, num_z, num_u, num_all) {
  try(if(num_all != num_u + 2*(1+num_z)*num_snps) stop("Dimensions don't match."))
  num_same_snp <- 2*(1+num_z) # number of features associated with the same snp 2(1+num_z)
  new_order <- list()
  for (i in 1:num_snps) {
    grp <- i + seq(0, (num_same_snp-1)*num_snps, by=num_snps)
    new_order[[i]] <- grp
  }
  new_order <- c(seq(1,num_u), unlist(new_order) + num_u) # 
  origin_order <- seq(1:num_all)
  recovery_index <- match(origin_order, new_order)
  
  # adelie
  grp_for_grplasso <- c(seq(1,num_u),seq(1, (1+num_z)*(2*num_snps-1)+1, by = 1+num_z)+num_u)
  
  return(list(new_order = new_order, 
              recovery_index = recovery_index,
              grp_for_grplasso = grp_for_grplasso))
}


#' Set Penalty Factors for Lasso
#'
#' Creates penalty factors vector for lasso regression where unpenalized
#' covariates have penalty factor 0, and all other predictors have penalty factor 1.
#'
#' @param unpenalized_covariates Matrix of covariates (n_samples x n_covariates).
#' @param genetic_variants Matrix of original genetic variants (n_samples x n_variants).
#' @param interaction_terms Interaction terms for genetic variants
#'   (n_samples x n_interaction_terms).
#' @param genetic_variants_knockoff Knockoff genetic variants
#'   (n_samples x n_variants).
#' @param interaction_terms_knockoff Interaction terms for knockoff genetic variants
#'   (n_samples x n_interaction_terms).
#'
#' @return Vector of penalty factors with length equal to total number of predictors.
#'
#' @keywords internal
set_penalty_factors <- function(unpenalized_covariates, genetic_variants, interaction_terms, genetic_variants_knockoff, interaction_terms_knockoff) {
  penalty_factors <- c(rep(0, ncol(unpenalized_covariates)), 
                       rep(1, ncol(genetic_variants) + ncol(interaction_terms) + 
                             ncol(genetic_variants_knockoff) + ncol(interaction_terms_knockoff)))
  return(penalty_factors)
}

#' Set Penalty Factors for Group Lasso
#'
#' Creates penalty factors vector for group lasso regression. Uses square-root
#' of group size as the default penalty for groups (except unpenalized covariates).
#'
#' @param num_snps Number of SNPs (equal to ncol(genetic_variants)).
#' @param num_z Number of interacting variables (number of PCs).
#' @param num_u Number of unpenalized covariates (sum of additional vars and Z).
#'
#' @return Vector of penalty factors with length equal to total number of groups.
#'
#' @keywords internal
set_penalty_factors_grpLasso <- function(num_snps, num_z, num_u) {
  # Default is square-root of group sizes for each group.
  penalty_for_grp <- sqrt(1+num_z) # penalty for groups except unpenalized covariates
  num_grp <- 2 * num_snps # number of groups except unpenalized covariates
  penalty_factors <- c(rep(0, num_u), rep(penalty_for_grp, num_grp))
  return(penalty_factors)
}

#' Fit Lasso Model with Cross-Validation
#'
#' Fits a lasso regression model with cross-validation using either MSE or
#' correlation as the metric.
#'
#' @param X_matrix Matrix of predictors (n_samples x n_predictors).
#' @param y_vector Response vector (length n_samples).
#' @param penalty_factors Vector of penalty factors (length n_predictors).
#' @param n_folds Number of folds for cross-validation (default: 10).
#' @param standardize Logical flag for standardization (default: FALSE).
#' @param metric Metric to optimize: "mse" for mean squared error or "r" for correlation (default: "mse").
#'
#' @return Cross-validated Lasso model object.
#'
#' @importFrom glmnet cv.glmnet
#' @keywords internal
fit_lasso_model <- function(X_matrix, y_vector, penalty_factors, n_folds = 10, standardize = FALSE, metric = "mse") {
  # Validate inputs
  if (!metric %in% c("mse", "r")) {
    stop("Invalid metric. Choose 'mse' or 'r'.")
  }
  
  if (metric == "mse") {
    cv_fit <- glmnet::cv.glmnet(x = X_matrix, 
                                 y = y_vector, 
                                 family = "gaussian", 
                                 alpha = 1, 
                                 penalty.factor = penalty_factors, 
                                 nfolds = n_folds, 
                                 standardize = standardize)
  } else if (metric == "r") {
    cv_fit <- cv_glmnet_correlation(X_matrix, y_vector, 
                                     nfolds = n_folds, 
                                     alpha = 1, 
                                     penalty.factor = penalty_factors, 
                                     standardize = standardize)
  }
  
  return(cv_fit)
}


#' Extract Coefficients from Lasso Model
#'
#' Extracts coefficients from a cross-validated lasso model at lambda.min.
#'
#' @param cv_fit Cross-validated Lasso model object (from cv.glmnet).
#'
#' @return Sparse matrix of coefficients at lambda.min.
#'
#' @keywords internal
extract_coefficients <- function(cv_fit) {
  # coef() is an S3 generic from stats package - R will dispatch to the correct method
  coefs <- coef(cv_fit, s = "lambda.min")
  message("Number of coefficients: ", length(coefs))
  return(coefs)
}


#' Fit Group Lasso Model with Cross-Validation
#'
#' Fits a group lasso regression model with cross-validation using the adelie package.
#'
#' @param X_matrix Matrix of predictors (n_samples x n_predictors).
#' @param y_vector Response vector (length n_samples).
#' @param penalty_factors Vector of penalty factors for each group.
#' @param groups Vector defining the starting column index of X for each group.
#' @param n_folds Number of folds for cross-validation (default: 10).
#' @param standardize Logical flag for standardization (default: FALSE).
#'
#' @return Cross-validated group lasso model object.
#'
#' @importFrom adelie cv.grpnet glm.gaussian
#' @keywords internal
fit_grp_lasso_model <- function(X_matrix, y_vector, penalty_factors, groups, n_folds = 10, standardize = FALSE) {
  # Using adelie package for group lasso
  # Reference: https://cran.r-project.org/web/packages/adelie/adelie.pdf
  cv_fit <- adelie::cv.grpnet(X = X_matrix, 
                               glm = adelie::glm.gaussian(y_vector), 
                               groups = groups,
                               alpha = 1, 
                               penalty = penalty_factors,
                               nfolds = n_folds, 
                               adev_tol = 0.999,
                               ddev_tol = 1e-5,
                               lmda_path_size = 100,
                               min_ratio = 0.0001,
                               standardize = standardize)
  
  # Alternative using gglasso (commented out):
  # Reference: https://cran.r-project.org/web/packages/gglasso/gglasso.pdf
  # cv_fit <- gglasso::cv.gglasso(x = X_matrix, 
  #                               y = y_vector, 
  #                               group = groups,
  #                               pf = penalty_factors,
  #                               loss = "ls",
  #                               pred.loss = "L2",
  #                               nfolds = n_folds)
  
  return(cv_fit)
}


#' Extract Coefficients from Group Lasso Model
#'
#' Extracts coefficients from a cross-validated group lasso model and reorders
#' them to match the original predictor order.
#'
#' @param cv_fit Cross-validated Group Lasso model object (from cv.grpnet).
#' @param recovery_index Index vector to recover original coefficient order.
#'
#' @return Vector of coefficients for original and knockoff SNPs and their interaction terms.
#'
#' @keywords internal
extract_grp_lasso_coefficients <- function(cv_fit, recovery_index) {
  # Note: if standardize=TRUE was used in fitting the grpnet object,
  # the coefficients reported are for the standardized inputs. However,
  # the predict function will apply the stored standardization to newx and give
  # the correct predictions.
  
  # Using adelie package
  # coef() is an S3 generic from stats package - R will dispatch to the correct method
  ordered_coefs <- coef(cv_fit, lambda = "lambda.min")$betas
  
  # Alternative using gglasso (commented out):
  # ordered_coefs <- gglasso::coef(cv_fit, s = "lambda.min")[-1,]
  
  coefs <- ordered_coefs[recovery_index]
  message("Number of coefficients: ", length(coefs))
  return(coefs)
}

#' Calculate Local Feature Importance
#'
#' Calculates local feature importance matrices for original and knockoff variants
#' based on extracted model coefficients.
#'
#' @param coefs Coefficients extracted from Lasso or Group Lasso model.
#' @param genetic_variants Genetic variants matrix (n_samples x n_variants).
#' @param genetic_variants_knockoff Knockoff genetic variants matrix (n_samples x n_variants).
#' @param interaction_terms Interaction terms for original genetic variants
#'   (n_samples x n_interaction_terms).
#' @param interaction_terms_knockoff Interaction terms for knockoff genetic variants
#'   (n_samples x n_interaction_terms).
#' @param unpenalized_covariates Matrix of covariates (not penalized)
#'   (n_samples x n_covariates).
#' @param Z Heterogeneity variable (e.g., EUR or PCs) (n_samples x n_Z_variables).
#' @param model Model type: "Lasso" or "Grplasso".
#'
#' @return List with components:
#'   \item{local_feature_importance}{Matrix of feature importance for original variants
#'     (n_samples x n_variants).}
#'   \item{local_feature_importance_knockoff}{Matrix of feature importance for knockoff variants
#'     (n_samples x n_variants).}
#'
#' @keywords internal
calculate_feature_importance <- function(coefs, genetic_variants, genetic_variants_knockoff, interaction_terms, interaction_terms_knockoff, unpenalized_covariates, Z, model) {
  n_genetic_variants <- ncol(genetic_variants)
  n_interaction_terms <- ncol(interaction_terms)
  n_genetic_variants_knockoff <- ncol(genetic_variants_knockoff) # actually, same with genetic_variants
  n_interaction_terms_knockoff <- ncol(interaction_terms_knockoff)
  n_covariates <- ncol(unpenalized_covariates)
  Z <- as.matrix(Z)
  n_heterogeneity <- ncol(Z)

  message("--- Calculating feature importance ...")
  message("n_genetic_variants: ", n_genetic_variants)
  message("n_interaction_terms: ", n_interaction_terms)
  message("n_genetic_variants_knockoff: ", n_genetic_variants_knockoff)
  message("n_interaction_terms_knockoff: ", n_interaction_terms_knockoff)
  message("n_unpenalized_covariates: ", n_covariates)
  message("n_heterogeneity: ", n_heterogeneity)
  
  # Get indices and coefficients for genetic variants
  
  beta_start <- if (model == "Lasso") {
    n_covariates + 2 # Accounting for intercept
  } else {
    n_covariates + 1
  }
  beta_end <- beta_start + n_genetic_variants - 1
  beta <- coefs[beta_start:beta_end]

  gamma_start <- beta_end + 1
  gamma_end <- gamma_start + n_interaction_terms - 1
  gamma <- matrix(0, nrow = n_genetic_variants, ncol = n_heterogeneity)
  for (k in seq_len(n_heterogeneity)) {
    gamma[, k] <- coefs[(gamma_start + (k - 1) * n_genetic_variants):(gamma_start + k * n_genetic_variants - 1)]
  }
 
  beta_knockoff_start <- gamma_end + 1
  beta_knockoff_end <- beta_knockoff_start + n_genetic_variants_knockoff - 1
  beta_knockoff <- coefs[beta_knockoff_start:beta_knockoff_end]

  gamma_knockoff_start <- beta_knockoff_end + 1
  gamma_knockoff_end <- gamma_knockoff_start + n_interaction_terms_knockoff - 1

  gamma_knockoff <- matrix(0, nrow = n_genetic_variants_knockoff, ncol = n_heterogeneity)
  for (k in seq_len(n_heterogeneity)) {
    gamma_knockoff[, k] <- coefs[(gamma_knockoff_start + (k - 1) * n_genetic_variants_knockoff):(gamma_knockoff_start + k * n_genetic_variants_knockoff - 1)]
  }

  message("Length of beta: ", length(beta))
  message("Length of gamma: ", length(gamma))
  message("Length of beta_knockoff: ", length(beta_knockoff))
  message("Length of gamma_knockoff: ", length(gamma_knockoff))

  # Compute local feature importance with defined heterogeneity Z
  local_feature_importance <- calculate_local_importance(beta, gamma, Z, genetic_variants)
  local_feature_importance_knockoff <- calculate_local_importance(beta_knockoff, gamma_knockoff, Z, genetic_variants_knockoff)

  message("Local feature importance dimensions: ", dim(local_feature_importance)[1], " samples x ", dim(local_feature_importance)[2], " SNPs")
  message("Local feature importance knockoff dimensions: ", dim(local_feature_importance_knockoff)[1], " samples x ", dim(local_feature_importance_knockoff)[2], " SNPs knockoff")

  return(list(local_feature_importance = local_feature_importance,
              local_feature_importance_knockoff = local_feature_importance_knockoff))
}



#' Compute Local Feature Importance with Matrix Operations
#'
#' Internal helper function to compute local feature importance using matrix operations.
#'
#' @param beta Coefficient vector for genetic variants (length n_variants).
#' @param gamma Coefficient matrix for interaction terms (n_variants x n_Z_variables).
#' @param Z Heterogeneity variable matrix (n_samples x n_Z_variables).
#' @param genetic_variants Genetic variants matrix (n_samples x n_variants), used for dimension checking.
#'
#' @return Matrix of local feature importance (n_samples x n_variants).
#'
#' @keywords internal
calculate_local_importance <- function(beta, gamma, Z, genetic_variants) {
  # Matrix multiplication for the interaction term (gamma * Z^T)
  interaction_effect <- gamma %*% t(Z)  # Result: n_genetic_variants x n_samples
  # Add beta to the interaction effect (vectorized addition)
  local_feature_importance <- abs(sweep(interaction_effect, 1, beta, "+"))  # Result: n_genetic_variants x n_samples
  # Transpose back to match the original shape (samples x SNPs)
  return(t(local_feature_importance))
}


#' Calculate MK (Model-X Knockoff) Statistics
#'
#' Internal helper function to calculate kappa and tau statistics for knockoff filtering.
#'
#' @param T_0 Matrix or vector of statistics for original features.
#' @param T_k Matrix or vector of statistics for knockoff features.
#' @param method Method for calculating tau: "max" or "median" (default: "median").
#'
#' @return Matrix with columns kappa and tau statistics.
#'
#' @keywords internal
MK.statistic <- function(T_0, T_k, method = 'median') {
  T_0 <- as.matrix(T_0)
  T_k <- as.matrix(T_k)
  T.temp <- cbind(T_0, T_k)
  T.temp[is.na(T.temp)] <- 0
  
  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1
  
  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max.nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}


MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));
  if(length(which(tau[b]>0))!=0){
    index_bound<-max(which(tau[b]>0))
    for(i in 1:length(b)){
      temp.index<-i:min(length(b),Rej.Bound,index_bound)
      if(length(temp.index)==0){next}
      q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
      if(i>Rej.Bound){break}
    }
    q[q>1]<-1
  }
  
  return(q)
}

#' Knockoff Filter Calculations
#'
#' Applies the knockoff filter to feature importance matrices to control FDR.
#'
#' @param local_feature_importance Matrix of feature importance for original variants
#'   (n_samples x n_variants).
#' @param local_feature_importance_knockoff Matrix of feature importance for knockoff variants
#'   (n_samples x n_variants).
#' @param FDR_rate False Discovery Rate threshold (default: 0.1).
#'
#' @return List containing:
#'   \item{S_ij}{Binary selection matrix (n_samples x n_variants).}
#'   \item{scaled_selection_matrix}{Scaled selection matrix (n_samples x n_variants).}
#'   \item{W}{W statistic matrix (original - knockoff) (n_samples x n_variants).}
#'   \item{q_values}{Matrix of q-values (n_samples x n_variants).}
#'
#' @keywords internal
knockoff_filter <- function(local_feature_importance, local_feature_importance_knockoff, FDR_rate = 0.1) {
  W <- local_feature_importance - local_feature_importance_knockoff
  q_values <- matrix(NA, nrow = nrow(W), ncol = ncol(W))  # Each element represents a q-value for individual i and SNP j
  
  # Calculate q-values for each individual across all SNPs
  for (i in seq_len(nrow(W))) {
    # Calculate MK statistics for individual i across all SNPs
    MK_stat <- MK.statistic(local_feature_importance[i, ], local_feature_importance_knockoff[i, ])
    q <- MK.q.byStat(MK_stat[,'kappa'], MK_stat[,'tau'], M = 1, Rej.Bound = 10000)
    q_values[i, ] <- q
  }
  # Create the selection matrix S_ij: S_ij = I(q_ij <= FDR_rate)
  S_ij <- q_values <= FDR_rate
  
  # Create scaled selection matrix using W matrix and selection matrix S_ij
  scaled_selection_matrix <- t((t(W) / apply(W, 2, max, na.rm = TRUE))) * S_ij
  
  return(list(S_ij = S_ij, scaled_selection_matrix = scaled_selection_matrix, W = W, q_values = q_values))
}



#' Get Importance Matrices
#'
#' Main function that generates scaled selection, selection, and W statistic matrices
#' using genetic variants and knockoff variants data with Lasso or Group Lasso models.
#'
#' @param genetic_variants Matrix of genetic variants (n_samples x n_variants).
#' @param genetic_variants_knockoff Matrix of knockoff genetic variants
#'   (n_samples x n_variants).
#' @param additional_covariates Matrix of additional covariates (n_samples x n_covariates),
#'   or NULL if not provided.
#' @param Z Vector or matrix representing heterogeneity variable(s) (e.g., EUR or PCs)
#'   (n_samples x n_Z_variables). Cannot be NULL.
#' @param y Outcome vector (length n_samples).
#' @param n_folds Number of folds for cross-validation (default: 10).
#' @param standardize Logical flag for x variable standardization prior to fitting
#'   the model sequence. The coefficients are always returned on the original scale
#'   (default: TRUE).
#' @param FDR_rate False Discovery Rate threshold (default: 0.1).
#' @param model Model type: "Lasso" or "Grplasso" (default: "Lasso").
#' @param metric Metric for Lasso model: "mse" for mean squared error or "r" for correlation
#'   (default: "mse"). Only used when model="Lasso".
#'
#' @return A list containing:
#'   \item{coefs}{Extracted model coefficients.}
#'   \item{scaled_selection_matrix}{Scaled selection matrix (n_samples x n_variants).}
#'   \item{selection_matrix}{Binary selection matrix (n_samples x n_variants).}
#'   \item{W_statistic_matrix}{W statistic matrix (n_samples x n_variants).}
#'   \item{q_values}{Matrix of q-values (n_samples x n_variants).}
#'
#' @export
get_importance_matrices <- function(genetic_variants, genetic_variants_knockoff, additional_covariates, Z, y, n_folds = 10, standardize = TRUE, FDR_rate = 0.1, model = "Lasso", metric = "mse") {
  if (is.null(Z)) {
    stop("The heterogeneity variable Z cannot be NULL. Please provide a valid input for Z.")
  }
  # Check for overlapping variables between additional_covariates and Z
  if (!is.null(additional_covariates)) {
    overlapping_vars <- intersect(colnames(as.data.frame(additional_covariates)), 
                                   colnames(as.data.frame(Z)))
    if (length(overlapping_vars) > 0) {
      stop("Overlap detected between additional_covariates and Z: ", 
           paste(overlapping_vars, collapse = ", "), 
           ". Please redefine Z and additional_covariates to avoid overlapping variables.")
    }
  }

  if (!is.null(additional_covariates)) {
    unpenalized_covariates <- cbind(additional_covariates, Z)
  } else {
    unpenalized_covariates <- Z
  }
  # Prepare interaction terms
  interaction_terms <- prepare_interaction_terms(genetic_variants, Z)
  interaction_terms_knockoff <- prepare_interaction_terms(genetic_variants_knockoff, Z)
  
  # Combine predictors
  X_matrix <- combine_predictors(unpenalized_covariates, genetic_variants, interaction_terms, genetic_variants_knockoff, interaction_terms_knockoff)

  # Prepare groups information for group lasso model
  grp_info <- prepare_feat_grps(ncol(genetic_variants), ncol(Z), ncol(unpenalized_covariates), ncol(X_matrix))
  new_order <- grp_info$new_order
  recovery_index <- grp_info$recovery_index
  grp_for_grplasso <- grp_info$grp_for_grplasso

  # Set penalty factors
  penalty_factors <- if (model == "Lasso") {
      set_penalty_factors(unpenalized_covariates, genetic_variants, interaction_terms, genetic_variants_knockoff, interaction_terms_knockoff)
    } else if (model == "Grplasso") {
      set_penalty_factors_grpLasso(ncol(genetic_variants), ncol(Z), ncol(unpenalized_covariates))
    } else {
      stop('Invalid model type. Choose either "Lasso" or "Grplasso".')
    }
  
  # Change column order for group lasso
  if (model == "Grplasso") {
    X_matrix <- X_matrix[, new_order]
  }
  y_vector <- as.numeric(y)
  
  # Scale X_matrix (including unpenalized factors) manually 
  if (standardize == TRUE && model == "Grplasso") {
    X_mean <- colMeans(X_matrix)
    X_sd <- apply(X_matrix, 2, sd)
    X_matrix <- t((t(X_matrix)-X_mean)/X_sd)
  }
  
  # Fit Lasso model, using n-fold cross validation
  message("\n--- Fitting ", model, " model...")
  
  fit_data_timing <- system.time({
    if (model == "Lasso") {
      cv_fit <- fit_lasso_model(X_matrix, y_vector, penalty_factors, n_folds, 
                                 standardize = standardize, metric = metric)
    } else if (model == "Grplasso") {
      cv_fit <- fit_grp_lasso_model(X_matrix, y_vector, penalty_factors, grp_for_grplasso, 
                                     n_folds, standardize = FALSE)
    } else {
      stop('Invalid model type. Choose either "Lasso" or "Grplasso".')
    }
  })
  
  message("--- Time taken for ", model, " model fitting: ", fit_data_timing[3], " seconds")
  
  # Extract coefficients
  coefs <- if (model == "Lasso") {
    extract_coefficients(cv_fit)
  } else if (model == "Grplasso") {
    extract_grp_lasso_coefficients(cv_fit, recovery_index)
  } else {
    stop('Invalid model type. Choose either "Lasso" or "Grplasso".')
  }
  
  # Calculate feature importance
  message("\n--- Calculating feature importance...")
  feature_timing <- system.time({
    feature_importances <- calculate_feature_importance(coefs, genetic_variants, 
                                                        genetic_variants_knockoff, 
                                                        interaction_terms, 
                                                        interaction_terms_knockoff, 
                                                        unpenalized_covariates, Z, model)
  })
  message("--- Time taken for feature importance calculation: ", feature_timing[3], " seconds")
  
  # Apply knockoff filter to get matrices
  message("\n--- Applying knockoff filter...")
  knockoff_timing <- system.time({
    matrices <- knockoff_filter(feature_importances$local_feature_importance, 
                                feature_importances$local_feature_importance_knockoff, FDR_rate)
  })
  message("--- Time taken for knockoff filter: ", knockoff_timing[3], " seconds")
  
  return(list(
    coefs = coefs,
    scaled_selection_matrix = matrices$scaled_selection_matrix,
    selection_matrix = matrices$S_ij,
    W_statistic_matrix = matrices$W,
    q_values = matrices$q_values
  ))
}


#' Plot PCs Visualization for a SNP
#'
#' Creates a scatter plot of two principal components colored by SNP feature importance.
#'
#' @param pcs Principal Components data (a matrix or data frame with at least two columns,
#'   such as PC1, PC2, etc.) (n_samples x n_pcs).
#' @param snp_importance Numeric vector representing feature importance values for the
#'   selected SNP (length n_samples).
#' @param snp_name Name of the SNP to visualize.
#' @param pc_x The name of the principal component to plot on the x-axis (e.g., "PC1").
#'   Must be a column name in pcs (default: "PC1").
#' @param pc_y The name of the principal component to plot on the y-axis (e.g., "PC2").
#'   Must be a column name in pcs (default: "PC2").
#' @param save_path The directory path where the plot will be saved
#'   (default: current working directory).
#'
#' @return A ggplot object of the specified PCs colored by feature importance.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient labs theme_minimal ggsave
#' @importFrom rlang .data
#' @export
plot_pcs <- function(pcs, snp_importance, snp_name, pc_x = "PC1", pc_y = "PC2", save_path = ".") {
  if (!(pc_x %in% colnames(pcs)) || !(pc_y %in% colnames(pcs))) {
    stop("The specified principal components (pc_x or pc_y) are not found in the pcs data.")
  }
  
  # Create plotting data
  pc_data <- data.frame(PC_X = pcs[[pc_x]], PC_Y = pcs[[pc_y]], Importance = snp_importance)
  
  plot <- ggplot(pc_data, aes(x = PC_X, y = PC_Y, color = Importance)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = paste("Scaled Selected Feature Importance of", snp_name, "in", pc_x, "&", pc_y, "Space"),
      x = pc_x,
      y = pc_y,
      color = "Importance"
    ) +
    theme_minimal()
  
  # Construct the file path and save the plot
  file_name <- file.path(save_path, paste0(snp_name, "_scaled_selected_feature_importance_", pc_x, "&", pc_y, ".png"))
  ggsave(file_name, plot = plot, width = 8, height = 8)
  
  message(snp_name, "'s PC plot is saved as ", file_name)
  return(plot)
}


#' Plot Heatmap of Selected SNPs
#'
#' Creates a heatmap of selected SNPs based on the scaled selection matrix,
#' sorted by EUR percentages for individuals.
#'
#' @param scaled_selection_matrix Scaled selection matrix from the knockoff filter
#'   (n_samples x n_variants).
#' @param genetic_variants Matrix of genetic variants (n_samples x n_variants).
#'   Column names are used to match with selected_snp_names.
#' @param selected_snp_names Character vector of names of selected SNPs to plot.
#' @param eur Numeric vector representing EUR percentages to sort individuals
#'   (length n_samples).
#' @param save_path The directory path where the plot will be saved
#'   (default: current working directory).
#'
#' @return A ggplot object representing the heatmap.
#'
#' @importFrom reshape2 melt
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme element_text labs ggsave
#' @importFrom rlang .data
#' @export
plot_heatmap <- function(scaled_selection_matrix, genetic_variants, selected_snp_names, eur, save_path = ".") {
  colnames(scaled_selection_matrix) <- colnames(genetic_variants)
  # Filter the scaled_selection_matrix to keep only the selected variants (SNPs)
  valid_snp_names <- intersect(colnames(scaled_selection_matrix), selected_snp_names)
  message("Number of valid SNP names: ", length(valid_snp_names))
  message("Number of selected SNP names: ", length(selected_snp_names))
  
  if (length(valid_snp_names) == 0) {
    stop("No valid SNP names found in the selection matrix.")
  }

  # Apply the filtering using only the valid SNP names
  scaled_selection_matrix_filtered <- scaled_selection_matrix[, valid_snp_names]

  # Sort the filtered selection matrix based on EUR
  sorted_indices <- order(eur)
  scaled_selection_matrix_sorted_filtered <- scaled_selection_matrix_filtered[sorted_indices, ]

  # Extract column(s) values for sorting SNPs
  chromosome_order <- as.numeric(gsub("chr([0-9]+)\\..*", "\\1", valid_snp_names))
  
  snp_data <- data.frame(SNP = valid_snp_names, Chromosome = chromosome_order)

  # Sort the SNP data frame by chromosome and then by base pair position
  snp_data <- snp_data[order(snp_data$Chromosome), ]
  valid_snp_names_sorted <- snp_data$SNP

  # Convert the filtered scaled selection matrix to a long-format data frame for ggplot2
  scaled_selection_matrix_long_filtered <- reshape2::melt(scaled_selection_matrix_sorted_filtered)
  colnames(scaled_selection_matrix_long_filtered) <- c("Individual", "SNP", "Importance")

  # Update SNP names after filtering, ensuring the SNP names match the number of rows in the melted data frame
  scaled_selection_matrix_long_filtered$SNP <- factor(scaled_selection_matrix_long_filtered$SNP, levels = valid_snp_names_sorted)

  # Create heatmap
  colors <- RColorBrewer::brewer.pal(11, "RdYlBu")
  # Use rlang::.data for proper variable scoping in ggplot2
  heatmap_plot <- ggplot2::ggplot(scaled_selection_matrix_long_filtered, 
                                    ggplot2::aes(x = rlang::.data$SNP, y = rlang::.data$Individual, fill = rlang::.data$Importance)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = c("grey", colors),
      values = scales::rescale(c(0, min(scaled_selection_matrix_long_filtered$Importance[scaled_selection_matrix_long_filtered$Importance > 0]), max(scaled_selection_matrix_long_filtered$Importance))),
      na.value = "grey", limits = c(0, max(scaled_selection_matrix_long_filtered$Importance))
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(
      title = "Heatmap of Selected SNPs (for At Least One Person)",
      x = "Genetic Variant (SNP)",
      y = "Individual (Ranked by EUR%)",
      fill = "Importance"
    )
  
  # Save the plot
  file_name <- file.path(save_path, "scaled_selected_feature_importance_heatmap.png")
  ggplot2::ggsave(file_name, plot = heatmap_plot, width = 30, height = 10)
  message("Heatmap plot is saved as ", file_name)
  
  return(heatmap_plot)
}
