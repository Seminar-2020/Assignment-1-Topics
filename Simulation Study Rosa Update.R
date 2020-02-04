library(mvtnorm)
# Contaminate observations with probability epsilon
# n         number of observations in-sample
# epsilon   contamination level
# R         number of repetitions
# perc      percentage of in-sample size to generate for out-of-sample evaluations
# p         numer of predictor variables

simulateOutliers <- function (n, epsilon, R, perc, p) {
  # control parameters
  set.seed(10)
  eps <- matrix(epsilon, ncol = p, nrow = n)
  intercept <- 1
  b <- rep(1,p)
  coefficients <- c(intercept,b)
  Sigma <- diag(p)
  # generate probability of being outlying for each observation
  poutvalue <- runif(n)
  pout <- matrix(poutvalue, ncol = p, nrow = n)
  
  sampleClean <- replicate(R, {
    # generate multivariate normal data
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)), sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    df <- data.frame(Y,X)
    # generate 10% extra observations to use oos
    X_oos <- rmvnorm(ceiling(n*perc), mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,ceiling(n*perc)),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(ceiling(n*perc), mean = 0, sd = 1)
    # compute estimators
    lm <- lm(Y ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate coefficients using RMSE
    beta_true <- coefficients
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    # evaluate prediction performance using RMSE
    lm_pred <- X_oos%*%beta_ols    
    lts_pred <- X_oos%*%beta_lts
    plugin_pred <- X_oos%*%beta_plugin
    RMSE_lm_oos <- RMSE(Y_oos, lm_pred)
    RMSE_lts_oos <- RMSE(Y_oos, lts_pred)
    RMSE_plugin_oos <- RMSE(Y_oos, plugin_pred)
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin,RMSE_lm_oos,RMSE_lts_oos,RMSE_plugin_oos)
  })
  
  sampleGLP <- replicate(R, {
    # generate multivariate normal data
    X <-  rmvnorm(n,  mean = rep(20,nrow(Sigma)), sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # generate 10% extra (regular) observations to use oos
    X_oos <- rmvnorm(ceiling(n*perc), mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,ceiling(n*perc)),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(ceiling(n*perc), mean = 0, sd = 1)
    # change some data points into good leverage points
    X_GLP <- ifelse(pout < eps, mvrnorm(n, mean = rep(40,p), sigma = diag(p)), X)
    Y_GLP <- intercept + X_GLP%*%b + rnorm(n, mean = 0, sd = 1)
    df <- data.frame(Y_GLP,X_GLP)
    # compute estimators
    lm <- lm(Y_GLP ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X_GLP,Y_GLP)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X_GLP,Y_GLP,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate estimates using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    # evaluate prediction performance using RMSE
    lm_pred <- X_oos%*%beta_ols    
    lts_pred <- X_oos%*%beta_lts
    plugin_pred <- X_oos%*%beta_plugin
    RMSE_lm_oos <- RMSE(Y_oos, lm_pred)
    RMSE_lts_oos <- RMSE(Y_oos, lts_pred)
    RMSE_plugin_oos <- RMSE(Y_oos, plugin_pred)
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin,RMSE_lm_oos,RMSE_lts_oos,RMSE_plugin_oos)
  })
  
  sampleBLP <- replicate(R, {
    # generate multivariate normal data
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # generate 10% extra observations to use oos
    X_oos <- rmvnorm(ceiling(n*perc), mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,ceiling(n*perc)),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(ceiling(n*perc), mean = 0, sd = 1)
    # change some points into bad leverage points
    X_BLP <- ifelse(pout < eps, mvrnorm(n, mean = rep(-40,p), sigma = diag(p)), X)
    Y_BLP <- intercept + X_BLP%*%-b + rnorm(n, mean = 40, sd = 1)
    df <- data.frame(Y_BLP,X_BLP)
    # compute estimators
    lm <- lm(Y_BLP ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X_BLP,Y_BLP)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X_BLP,Y_BLP,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate using RMSE
    beta_true <- rbind(intercept, b)
    RMSE_lm <- RMSE(beta_true, beta_ols)
    RMSE_lts <- RMSE(beta_true, beta_lts)
    RMSE_plugin <- RMSE(beta_true, beta_plugin)
    
    c(beta_ols,beta_lts,beta_plugin,RMSE_lm,RMSE_lts,RMSE_plugin)
  }) 
  
  sampleVO <- replicate(R, {
    # generate multivariate normal data
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # change some points into vertical outliers
    Y_VO <- ifelse(pout < eps, intercept + X%*%b + rnorm(n, mean = 40, sd = 1), Y)
    df <- data.frame(Y_VO,X)
    # compute estimators
    lm <- lm(Y_VO ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y_VO)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y_VO,alpha=0.5)
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
    # generate multivariate normal data
    X <- rmvnorm(n, sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    eps <- eps/3
    # replace some point by good leverage points
    X_GLP <- ifelse(pout < eps, mvrnorm(n, mean = rep(40,p), sigma = diag(p)), X)
    Y_GLP <- intercept + X_GLP%*%b + rnorm(n, mean = 0, sd = 1)
    # replace some points by bad leverage points
    X_BLP <- ifelse(eps < pout & pout < 2*eps, mvrnorm(n, mean = rep(-40,p), sigma = diag(p)), X_GLP)
    Y_BLP <- ifelse(eps < pout & pout < 2*eps, intercept + X_BLP%*%-b + rnorm(n, mean = 40, sd = 1), Y_GLP)
    # replace some points by vertical outliers
    Y_VO <- ifelse(2*eps < pout & pout < 3*eps, intercept + X_BLP%*%b + rnorm(n, mean = 20, sd = 1), Y_BLP)
    Y <- Y_VO
    X <- X_BLP
    df <- data.frame(Y,X)
    # compute estimators
    lm <- lm(Y ~ ., df)
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
  MSE <- mean((true-estimated)^2)
  RMSE <- sqrt(MSE)
  return (RMSE)
}
