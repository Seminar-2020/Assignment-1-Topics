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
  eps <- rep(epsilon, n)
  intercept <- 1
  b <- rep(1,p)
  coefficients <- c(intercept,b)
  Sigma <- diag(p)
  # generate probability of being outlying for each observation
  pout <- runif(n)
  
  sampleClean <- replicate(R, {
    # generate clean data
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
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
    # generate good leverage points
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
    # generate 10% extra (regular) observations to use oos
    X_oos <- rmvnorm(ceiling(n*perc), mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,ceiling(n*perc)),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(ceiling(n*perc), mean = 0, sd = 1)
    # change some points into good leverage points
    index <- pout < eps
    
    for (i in 1:sum(index)) {
      X[] <- rnorm(p, mean = 40, sd = 1)
    }
    #X[index,] <- rnorm(p, mean = 40, sd = 1) #gaat fout want vervangt allemaal met dezelfde waarde!!! 
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    df <- data.frame(Y,X)
    # compute estimators
    lm <- lm(Y ~.,df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y,alpha=0.5)
    beta_plugin <- plugin$coefficients
    # evaluate using RMSE
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
  
  sampleBLP <- replicate(R, {
    # generate bad leverage points
    # X <- rmvnorm(n, sigma = Sigma)
    # Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # change some points into bad leverage points
    index <- pout < eps
    X[index,] <- rnorm(p, mean = 40, sd = 1) 
    Y[index,] <- intercept + X[index,]%*%-b + rnorm(n, mean = 40, sd = 1)
    df <- data.frame(Y,X)
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
