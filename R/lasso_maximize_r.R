# Custom cross-validation for Lasso regression minimizing Pearson correlation
# Instead of minimizing MSE like cv.glmnet(), we minimize correlation between predictions and outcomes

# Function to calculate Pearson correlation (negative because we want to minimize this to maximize correlation)
calc_correlation <- function(y_true, y_pred) {
  # Remove any NA values
  valid_idx <- !is.na(y_true) & !is.na(y_pred)
  if (sum(valid_idx) < 2) return(Inf)
  
  y_true_clean <- y_true[valid_idx]
  y_pred_clean <- y_pred[valid_idx]
  
  # Check if predictions are constant (which would cause correlation issues)
  if (length(unique(y_pred_clean)) <= 1) {
    return(Inf)  # Return Inf for constant predictions
  }
  
  # Check if true values are constant
  if (length(unique(y_true_clean)) <= 1) {
    return(Inf)  # Return Inf for constant true values
  }
  
  # Calculate Pearson correlation
  cor_val <- cor(y_true_clean, y_pred_clean, method = "pearson")
  
  # Return negative correlation (so minimizing this maximizes positive correlation)
  # Handle cases where correlation is NA or infinite
  if (is.na(cor_val) || is.infinite(cor_val)) {
    return(Inf)
  }
  
  return(-cor_val)
}

# Custom cross-validation function for Lasso minimizing correlation
cv_glmnet_correlation <- function(x, y, 
                                  nfolds = 10, 
                                  lambda = NULL,
                                  alpha = 1,  # 1 for Lasso, 0 for Ridge, 0-1 for Elastic Net
                                  family = "gaussian",
                                  foldid = NULL,
                                  seed = 123,
                                  penalty.factor = rep(1, ncol(x)),
                                  standardize = FALSE,
                                  ...) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # If lambda not provided, generate a sequence similar to standard cv.glmnet
  if (is.null(lambda)) {
    # Use the same lambda generation logic as standard cv.glmnet
    n <- nrow(x)
    p <- ncol(x)
    
    # Default parameters similar to cv.glmnet
    lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
    nlambda <- 100
    
    # Calculate lambda sequence
    lambda_max <- max(abs(t(x) %*% y)) / n
    lambda_min <- lambda_max * lambda.min.ratio
    
    # Generate lambda sequence on log scale
    lambda <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
  }
  
  # If foldid not provided, create random folds
  if (is.null(foldid)) {
    n <- length(y)
    foldid <- sample(rep(seq(nfolds), length.out = n))
  }
  
  # Standardize predictors if requested
  if (standardize) {
    # Calculate standardization parameters from full data
    x_means <- colMeans(x)
    x_sds <- apply(x, 2, sd)
    
    # Avoid division by zero
    x_sds[x_sds == 0] <- 1
    
    # Store standardization parameters for later use
    standardization_params <- list(means = x_means, sds = x_sds)
  } else {
    standardization_params <- NULL
  }
  
  # Initialize storage for cross-validation results
  cv_cor <- matrix(NA, nrow = length(lambda), ncol = nfolds)
  cv_cor_mean <- numeric(length(lambda))
  cv_cor_se <- numeric(length(lambda))
  
  # Perform cross-validation
  for (fold in 1:nfolds) {
    # Split data into training and validation sets
    train_idx <- foldid != fold
    val_idx <- foldid == fold
    
    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    x_val <- x[val_idx, , drop = FALSE]
    y_val <- y[val_idx]
    
    # Apply standardization if requested
    if (standardize) {
      # Standardize training data
      x_train_std <- scale(x_train, center = x_means, scale = x_sds)
      
      # Standardize validation data using training parameters
      x_val_std <- scale(x_val, center = x_means, scale = x_sds)
      
      # Fit model on standardized training data
      fit <- glmnet(x_train_std, y_train, alpha = alpha, family = family, lambda = lambda, 
                    penalty.factor = penalty.factor, standardize = FALSE, ...)
      
      # Make predictions on standardized validation data
      pred_val <- predict(fit, newx = x_val_std)
    } else {
      # Fit model on original data
      fit <- glmnet(x_train, y_train, alpha = alpha, family = family, lambda = lambda, 
                    penalty.factor = penalty.factor, ...)
      
      # Make predictions on original validation data
      pred_val <- predict(fit, newx = x_val)
    }
    
    # Calculate correlation for each lambda
    for (i in 1:length(lambda)) {
      cv_cor[i, fold] <- calc_correlation(y_val, pred_val[, i])
    }
  }
  
  # Calculate mean and standard error across folds
  cv_cor_mean <- apply(cv_cor, 1, mean, na.rm = TRUE)
  cv_cor_se <- apply(cv_cor, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  
  # Find optimal lambda (minimum negative correlation, which maximizes positive correlation)
  lambda_min <- lambda[which.min(cv_cor_mean)]
  lambda_1se <- lambda[which(cv_cor_mean <= min(cv_cor_mean) + cv_cor_se[which.min(cv_cor_mean)])[1]]
  
  # Create result object similar to cv.glmnet
  result <- list(
    lambda = lambda,
    cvm = cv_cor_mean,      # Mean correlation across folds
    cvsd = cv_cor_se,       # Standard error of correlation
    cvup = cv_cor_mean + cv_cor_se,  # Upper bound
    cvlo = cv_cor_mean - cv_cor_se,  # Lower bound
    lambda.min = lambda_min, # Lambda with minimum correlation
    lambda.1se = lambda_1se, # Lambda within 1 SE of minimum
    nfolds = nfolds,
    foldid = foldid,
    alpha = alpha,
    family = family
  )
  
  # Fit final model with optimal lambda to get coefficients
  if (standardize) {
    # Standardize full data
    x_std <- scale(x, center = x_means, scale = x_sds)
    
    # Fit model on standardized data
    final_model <- glmnet(x_std, y, alpha = alpha, family = family, 
                          lambda = lambda_min, penalty.factor = penalty.factor, 
                          standardize = FALSE, ...)
    
    # Extract coefficients (intercept first, then variable coefficients)
    # coef() is an S3 generic from stats package - R will dispatch to the correct method
    coef_matrix <- coef(final_model)
    coefficients <- as.numeric(coef_matrix)
    names(coefficients) <- rownames(coef_matrix)
    
    # Store standardization parameters
    result$standardization_params <- standardization_params
  } else {
    # Fit model on original data
    final_model <- glmnet(x, y, alpha = alpha, family = family, 
                          lambda = lambda_min, penalty.factor = penalty.factor, ...)
    
    # Extract coefficients (intercept first, then variable coefficients)
    # coef() is an S3 generic from stats package - R will dispatch to the correct method
    coef_matrix <- coef(final_model)
    coefficients <- as.numeric(coef_matrix)
    names(coefficients) <- rownames(coef_matrix)
  }
  
  # Add coefficients to result
  result$coefficients <- coefficients
  result$glmnet.fit <- final_model
  
  # Add method to get coefficients at specific lambda
  result$get_coefficients <- function(lambda_val = NULL) {
    if (is.null(lambda_val)) {
      lambda_val <- lambda_min
    }
    # coef() is an S3 generic from stats package - R will dispatch to the correct method
    coef_at_lambda <- coef(result$glmnet.fit, s = lambda_val)
    return(as.numeric(coef_at_lambda))
  }
  
  # Add method to get unstandardized coefficients if standardization was used
  result$get_unstandardized_coefficients <- function(lambda_val = NULL) {
    if (is.null(lambda_val)) {
      lambda_val <- lambda_min
    }
    
    if (!standardize || is.null(standardization_params)) {
      return(result$get_coefficients(lambda_val))
    }
    
    # Get standardized coefficients
    std_coef <- result$get_coefficients(lambda_val)
    
    # Convert to unstandardized coefficients
    unstd_coef <- std_coef
    unstd_coef[-1] <- std_coef[-1] / standardization_params$sds  # Variable coefficients
    unstd_coef[1] <- std_coef[1] - sum(std_coef[-1] * standardization_params$means)  # Intercept
    
    return(unstd_coef)
  }
  
  class(result) <- "cv.glmnet.correlation"
  return(result)
}

# Plot function for correlation-based CV results
plot.cv.glmnet.correlation <- function(x, main = NULL, ...) {
  # Plot negative correlation directly (what we're minimizing)
  cvm_neg <- x$cvm
  cvup_neg <- x$cvup
  cvlo_neg <- x$cvlo
  
  # Set default title if not provided
  if (is.null(main)) {
    main <- "Cross-validation: Negative Pearson Correlation vs Lambda"
  }
  
  # Ensure we use the exact same lambda values and limits as standard CV
  lambda_plot <- x$lambda
  xlim_range <- range(lambda_plot)
  
  # Create plot with log scale for lambda, using exact same x-axis as standard CV
  plot(lambda_plot, cvm_neg, type = "p", 
       xlab = "Lambda (log scale)", ylab = "Negative Pearson Correlation",
       main = main,
       xlim = log(xlim_range), ...)
  
  # Add confidence intervals
  lines(lambda_plot, cvup_neg, lty = 2, col = "red")
  lines(lambda_plot, cvlo_neg, lty = 2, col = "red")
  
  # Add optimal lambda lines for custom CV
  abline(v = x$lambda.min, lty = 3, col = "blue", lwd = 2)
  abline(v = x$lambda.1se, lty = 3, col = "green", lwd = 2)
  
  # Add legend
  legend("topright", 
         legend = c("Mean Negative Correlation", "Upper/Lower CI", 
                    paste("lambda.min =", round(x$lambda.min, 4)),
                    paste("lambda.1se =", round(x$lambda.1se, 4))),
         lty = c(1, 2, 3, 3), 
         col = c("black", "red", "blue", "green"))
}

# Example usage function
example_usage <- function() {
  # Generate sample data
  set.seed(123)
  n <- 100
  p <- 20
  
  # Create correlated predictors
  x <- matrix(rnorm(n * p), n, p)
  # Add some correlation structure
  x[, 2] <- x[, 1] + rnorm(n, 0, 0.1)
  x[, 3] <- x[, 1] + rnorm(n, 0, 0.1)
  
  # True coefficients (sparse)
  beta_true <- c(2, -1.5, 1, rep(0, p-3))
  
  # Generate outcome with noise
  y <- x %*% beta_true + rnorm(n, 0, 0.5)
  
  cv_mse <- cv.glmnet(x, y, nfolds = 5, alpha = 1)
  
  # Use the lambda path from standard CV
  lambda <- cv_mse$lambda
  
  # Run correlation-based cross-validation
  cv_result <- cv_glmnet_correlation(x, y, nfolds = 5, alpha = 1, lambda = NULL)
  
  # Run with standardization
  cv_result_std <- cv_glmnet_correlation(x, y, nfolds = 5, alpha = 1, lambda = NULL, standardize = TRUE)
  
  # Plot results
  plot(cv_result)
  
  # Print optimal lambdas
  cat("Custom CV (Correlation):\n")
  cat("  Lambda with maximum correlation:", cv_result$lambda.min, "\n")
  cat("  Lambda within 1 SE of maximum:", cv_result$lambda.1se, "\n")
  
  cat("\nStandard CV (MSE):\n")
  cat("  Lambda with minimum MSE:", cv_mse$lambda.min, "\n")
  cat("  Lambda within 1 SE of minimum MSE:", cv_mse$lambda.1se, "\n")
  
  # Fit final models with both optimal lambdas
  final_model_custom <- glmnet(x, y, alpha = 1, lambda = cv_result$lambda.min)
  final_model_standard <- glmnet(x, y, alpha = 1, lambda = cv_mse$lambda.min)
  
  # Make predictions
  y_pred_custom <- predict(final_model_custom, newx = x)
  y_pred_standard <- predict(final_model_standard, newx = x)
  
  # Calculate final correlations
  final_cor_custom <- cor(y, y_pred_custom)
  final_cor_standard <- cor(y, y_pred_standard)
  
  cat("\nFinal correlations on full data:\n")
  cat("  Custom CV (correlation):", final_cor_custom, "\n")
  cat("  Standard CV (MSE):", final_cor_standard, "\n")
  
  # Add final mse
  final_mse_custom <- mean((y - y_pred_custom)^2)
  final_mse_standard <- mean((y - y_pred_standard)^2)
  
  cat("\nFinal MSE on full data:\n")
  cat("  Custom CV (MSE):", final_mse_custom, "\n")
  cat("  Standard CV (MSE):", final_mse_standard, "\n")
  
  # Show coefficients from custom CV
  cat("\nCoefficients from Custom CV (lambda.min):\n")
  cat("  Intercept:", round(cv_result$coefficients[1], 4), "\n")
  cat("  Non-zero coefficients:\n")
  non_zero_coef <- cv_result$coefficients[-1][cv_result$coefficients[-1] != 0]
  if (length(non_zero_coef) > 0) {
    for (i in 1:length(non_zero_coef)) {
      cat("    Variable", names(non_zero_coef)[i], ":", round(non_zero_coef[i], 4), "\n")
    }
  } else {
    cat("    All coefficients are zero\n")
  }
  
  # Show standardized vs unstandardized coefficients
  cat("\nStandardized vs Unstandardized Coefficients:\n")
  cat("  Original CV (no standardization):\n")
  cat("    Intercept:", round(cv_result$coefficients[1], 4), "\n")
  
  cat("  Standardized CV:\n")
  cat("    Standardized coefficients (intercept):", round(cv_result_std$coefficients[1], 4), "\n")
  cat("    Unstandardized coefficients (intercept):", round(cv_result_std$get_unstandardized_coefficients()[1], 4), "\n")
  
  
  
  return(list(cv_result = cv_result, 
              cv_result_std = cv_result_std,
              final_model_custom = final_model_custom, 
              final_model_standard = final_model_standard,
              final_correlation_custom = final_cor_custom,
              final_correlation_standard = final_cor_standard))
}

# Function to compare with standard MSE-based CV
compare_cv_methods <- function(x, y, nfolds = 5, alpha = 1) {
  # Standard MSE-based CV
  cv_mse <- cv.glmnet(x, y, nfolds = nfolds, alpha = alpha)
  
  # Correlation-based CV
  cv_cor <- cv_glmnet_correlation(x, y, nfolds = nfolds, alpha = alpha)
  
  # Get the exact same lambda range for both plots
  lambda_range <- range(cv_mse$lambda)
  
  # Create comparison plot with identical x-axes
  par(mfrow = c(1, 2))
  
  # Plot MSE-based CV
  plot(cv_mse, main = "Standard CV (MSE)", xlim = log(lambda_range))
  
  # Plot correlation-based CV
  plot.cv.glmnet.correlation(cv_cor, main = "Custom CV (Correlation)")
  
  par(mfrow = c(1, 1))
  
  # Compare optimal lambdas
  cat("MSE-based CV:\n")
  cat("  lambda.min:", cv_mse$lambda.min, "\n")
  cat("  lambda.1se:", cv_mse$lambda.1se, "\n")
  
  cat("\nCorrelation-based CV:\n")
  cat("  lambda.min:", cv_cor$lambda.min, "\n")
  cat("  lambda.1se:", cv_cor$lambda.1se, "\n")
  
  return(list(mse_cv = cv_mse, correlation_cv = cv_cor))
}

# # Print usage instructions
# cat("Custom Lasso CV with Pearson Correlation Maximization\n")
# cat("==================================================\n\n")
# cat("Main functions:\n")
# cat("- cv_glmnet_correlation(): Main CV function maximizing correlation\n")
# cat("- plot.cv.glmnet.correlation(): Plot CV results\n")
# cat("- example_usage(): Run example with sample data\n")
# cat("- compare_cv_methods(): Compare with standard MSE-based CV\n\n")
# cat("Usage example:\n")
# cat("result <- example_usage()\n")
# cat("comparison <- compare_cv_methods(x, y)\n\n")
# 
# # Automatically run the example when script is sourced
# cat("Running example automatically...\n")
# cat("================================\n")
# example_result <- example_usage()
# 
# # Also run the comparison
# cat("\n\nComparing with standard MSE-based CV...\n")
# cat("========================================\n")
# # Generate sample data for comparison
# set.seed(123)
# n <- 100
# p <- 20
# x_comp <- matrix(rnorm(n * p), n, p)
# x_comp[, 2] <- x_comp[, 1] + rnorm(n, 0, 0.1)
# x_comp[, 3] <- x_comp[, 1] + rnorm(n, 0, 0.1)
# beta_true <- c(2, -1.5, 1, rep(0, p-3))
# y_comp <- x_comp %*% beta_true + rnorm(n, 0, 0.5)
# 
# comparison_result <- compare_cv_methods(x_comp, y_comp, nfolds = 5, alpha = 1)
# 
# cat("\n\nExample completed! Results stored in 'example_result' and 'comparison_result' variables.\n")
