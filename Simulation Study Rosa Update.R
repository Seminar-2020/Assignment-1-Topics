library(mvtnorm)
# Contaminate observations with probability epsilon
simulateOutliers <- function (n, epsilon, R) {
  # control parameters
  set.seed(10)
  eps <- rep(epsilon, n)
  intercept <- 1
  b <- matrix(c(1, 1))
  Sigma <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  # generate probability of being outlying for each observation
  pout <- runif(n)
  
  sampleClean <- replicate(R, {
    # generate clean data
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    df <- data.frame(Y = Y, X1 = X[,1], X2 = X[,2])
    # generate 10% extra observations and cutoff to use oos
    X <- rmvnorm(ceiling(1.1*n), sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(ceiling(1.1*n), mean = 0, sd = 1)
    full <- data.frame(Y = Y, X1 = X[,1], X2 = X[,2])
    df <- full[1:n,]
    df_oos <- full[n+1:nrow(full),]  
    # compute estimators
    lm <- lm(df[,1] ~ df[,2] + df[,3], df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate coefficients using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    # evaluate prediction performance using RMSE
    y_true <- df_oos[,1]
    lm_pred <-     
    lts_pred <-
    plugin_pred <-
    RMSE_lm_oos <- RMSE(y_true, lm_pred)
    RMSE_lts_oos <- RMSE(y_true, lts_pred)
    RMSE_plugin_oos <- RMSE(y_true, plugin_pred)
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin)
  })
  
  sampleGLP <- replicate(R, {
    # generate good leverage points
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    x1 <- ifelse(pout < eps, rnorm(n, mean = 20, sd = 1), X[,1])
    x2 <- ifelse(pout < eps, rnorm(n, mean = 20, sd = 1), X[,2])
    X_new <- cbind(x1, x2)
    y <- intercept + X_new%*%b + rnorm(n, mean = 0, sd = 1)
    # compute estimators
    lm <- lm(y ~ x1 + x2)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X_new,y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X_new,y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin)
  })
  
  sampleBLP <- replicate(R, {
    # generate bad leverage points
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    x1 <- ifelse(pout < eps, rnorm(n, mean = -20, sd = 1), X[,1])
    x2 <- ifelse(pout < eps, rnorm(n, mean = -20, sd = 1), X[,2]) 
    X_new <- cbind(x1, x2)
    y <- ifelse(pout < eps, intercept + X_new%*%-b + rnorm(n, mean = 20, sd = 1), Y)
    # compute estimators
    lm <- lm(y ~ x1 + x2)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X_new,y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X_new,y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin)
  }) 
  
  sampleVO <- replicate(R, {
    # generate vertical outliers
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    y <- ifelse(pout < eps, intercept + X%*%b + rnorm(n, mean = 20, sd = 1), Y)
    # compute estimators
    lm <- lm(y ~ X[,1] + X[,2])
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin)
  })
  
  # generate sample contaminated with good leverage points, bad leverage points & vertical outliers
  sampleContam <- replicate(R,{
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    eps <- eps/3
    # add good leverage points
    x1_GLP <- ifelse(pout < eps, rnorm(n, mean = 20, sd = 1), X[,1])
    x2_GLP <- ifelse(pout < eps, rnorm(n, mean = 20, sd = 1), X[,2])
    X_GLP <- cbind(x1_GLP, x2_GLP)
    Y_GLP <- intercept + X_GLP%*%b + rnorm(n, mean = 0, sd = 1)
    # add bad leverage points
    x1_BLP <- ifelse(eps < pout & pout < 2*eps, rnorm(n, mean = -20, sd = 1), x1_GLP)
    x2_BLP <- ifelse(eps < pout & pout < 2*eps, rnorm(n, mean = -20, sd = 1), x2_GLP) 
    X_BLP <- cbind(x1_BLP, x2_BLP)
    Y_BLP <- ifelse(eps < pout & pout < 2*eps, intercept + X_BLP%*%-b + rnorm(n, mean = 20, sd = 1), Y_GLP)
    # add vertical outliers
    Y_VO <- ifelse(2*eps < pout & pout < 3*eps, intercept + X_BLP%*%b + rnorm(n, mean = 20, sd = 1), Y_BLP)
    # compute estimators
    Y <- Y_VO
    X <- X_BLP
    lm <- lm(Y ~ X[,1] + X[,2])
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin)
  })
  
  results <- list(
    "Clean" = sampleClean,
    "GLP" = sampleGLP,
    "BLP" = sampleBLP,
    "VO" = sampleVO,
    "allOutliers" = sampleContam
  )
  return(results)
}

RMSE <- function (true, estimated) {
  MSE <- mean(se(true, estimated))
  RMSE <- sqrt(MSE)
  return (RMSE)
}