library(mvtnorm)
# calculate the Root Mean Squared (Prediction) Error
RMSE <- function (true, estimated) {
  MSE <- mean((true-estimated)^2)
  RMSE <- sqrt(MSE)
  return (RMSE)
}

# Contaminate observations with probability epsilon
# Input:
# n ....... number of observations in-sample
# epsilon . contamination level
# R ....... number of repetitions
# alpha ... size initial subset used for MCD and LTS
# p ....... numer of predictor variables

# Output:
# Clean ........ the average estimated coefficients, RMSE and RMSPE for clean data
# GLP .......... the average estimated coefficients, RMSE and RMSPE for data with good leverage points
# BLP .......... the average estimated coefficients, RMSE and RMSPE for data with bad leverage points
# VO ........... the average estimated coefficients, RMSE and RMSPE for data with vertical outliers
# allOutliers .. the average estimated coefficients, RMSE and RMSPE for data containing all types of outliers
# Configuration  the values for n, epsilon, R, alpha and p used in the specific simulation
simulateOutliers <- function (n, epsilon, R, alpha, p) {
  # control parameters
  set.seed(10)
  eps <- matrix(epsilon, ncol = p, nrow = n)
  eps_vector <- rep(epsilon,n)
  intercept <- 1
  b <- c(2,rep(1,p-1))
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
    X_oos <- rmvnorm(n, mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,n),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(n, mean = 0, sd = 1)
    # compute estimators
    lm <- lm(Y ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y,alpha=alpha)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y,alpha=alpha)
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
  
  average_sampleClean <- apply(sampleClean,1,mean)
  
  sampleGLP <- replicate(R, {
    # generate multivariate normal data
    X <-  rmvnorm(n,  mean = rep(20,nrow(Sigma)), sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # generate 10% extra (regular) observations to use oos
    X_oos <- rmvnorm(n, mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,n),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(n, mean = 0, sd = 1)
    # change some data points into good leverage points
    X_GLP <- ifelse(pout < eps, rmvnorm(n, mean = rep(30,p), sigma = diag(p)), X)
    Y_GLP <- intercept + X_GLP%*%b + rnorm(n, mean = 0, sd = 1)
    df <- data.frame(Y_GLP,X_GLP)
    # compute estimators
    lm <- lm(Y_GLP ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X_GLP,Y_GLP,alpha=alpha)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X_GLP,Y_GLP,alpha=alpha)
    beta_plugin <- plugin$coefficients
    # evaluate estimates using RMSE
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
  
  average_sampleGLP <- apply(sampleGLP,1,mean)
  
  sampleBLP <- replicate(R, {
    # generate multivariate normal data
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # generate 10% extra observations to use oos
    X_oos <- rmvnorm(n, mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,n),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(n, mean = 0, sd = 1)
    # change some points into bad leverage points
    X_BLP <- ifelse(pout < eps, rmvnorm(n, mean = rep(30,p), sigma = diag(p)), X)
    Y_BLP <- ifelse(poutvalue < eps_vector, intercept + X_BLP%*%b + rnorm(n, mean = -40, sd = 1),Y)
    df <- data.frame(Y_BLP,X_BLP)
    # compute estimators
    lm <- lm(Y_BLP ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X_BLP,Y_BLP,alpha=alpha)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X_BLP,Y_BLP,alpha=alpha)
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
  
  average_sampleBLP <- apply(sampleBLP,1,mean)
  
  sampleVO <- replicate(R, {
    # generate multivariate normal data
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
    Y <- intercept + X%*%b + rnorm(n, mean = 0, sd = 1)
    # generate 10% extra observations to use oos
    X_oos <- rmvnorm(n, mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,n),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(n, mean = 0, sd = 1)
    # change some points into vertical outliers
    Y_VO <- ifelse(poutvalue < eps_vector, intercept + X%*%b + rnorm(n, mean = 10, sd = 1), Y)
    df3 <- data.frame(Y=Y_VO,X)
    # compute estimators
    lm <- lm(Y_VO ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y_VO,alpha=alpha)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y_VO,alpha=alpha)
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
  
  average_sampleVO <- apply(sampleVO,1,mean)
  
  # generate sample contaminated with good leverage points, bad leverage points & vertical outliers
  sampleContam <- replicate(R,{ 
    # generate multivariate normal data
    X <- rmvnorm(n,  mean = rep(20,nrow(Sigma)),sigma = Sigma)
    # generate 10% extra observations to use oos
    X_oos <- rmvnorm(n, mean = rep(20,nrow(Sigma)),sigma = Sigma)
    X_oos <- cbind(rep(1,n),X_oos)
    Y_oos <- X_oos%*%coefficients+rnorm(n, mean = 0, sd = 1)
    # all types of outlying points occur at same rate: epsilon/3
    eps <- eps/3
    eps_vector <- eps_vector/3
    # replace some point by good leverage points
    X_GLP <- ifelse(pout <= eps, rmvnorm(n, mean = rep(30,p), sigma = diag(p)), X)
    Y_GLP <- intercept + X_GLP%*%b + rnorm(n, mean = 0, sd = 1)
    # replace some points by bad leverage points
    X_BLP <- ifelse(eps < pout & pout <= 2*eps, rmvnorm(n, mean = rep(30,p), sigma = diag(p)), X_GLP)
    Y_BLP <- ifelse(eps_vector < poutvalue & poutvalue <= 2*eps_vector, intercept + X_BLP%*%b + rnorm(n, mean = -40, sd = 1), Y_GLP)
    # replace some points by vertical outliers
    Y_VO <- ifelse(2*eps_vector < poutvalue & poutvalue <= 3*eps_vector, intercept + X_BLP%*%b + rnorm(n, mean = 10, sd = 1), Y_BLP)
    Y <- Y_VO
    X <- X_BLP
    df4 <- data.frame(Y,X)
    # compute estimators
    lm <- lm(Y ~ ., df)
    beta_ols <- lm$coefficients
    lts <- ltsReg(X,Y,alpha=alpha)
    beta_lts <- lts$coefficients
    plugin <- lmDetMCD(X,Y,alpha=alpha)
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
  
  average_sampleContam <- apply(sampleContam,1,mean)
  
  results <- list(
    "Clean" = average_sampleClean,
    "GLP" = average_sampleGLP,
    "BLP" = average_sampleBLP,
    "VO" = average_sampleVO,
    "allOutliers" = average_sampleContam,
    "configuration" = c(n, epsilon, R, alpha, p)
  )
  return(results)
}

# Perform the simulations for different configurations
sim_1 <- simulateOutliers(100,0.1,200,0.75,1)
sim_2 <- simulateOutliers(100,0.3,200,0.75,1)
sim_3 <- simulateOutliers(100,0.5,200,0.75,1)
sim_4 <- simulateOutliers(100,0.1,200,0.75,5)
sim_5 <- simulateOutliers(100,0.3,200,0.75,5)
sim_6 <- simulateOutliers(100,0.5,200,0.75,5)
sim_7 <- simulateOutliers(100,0.1,200,0.75,40)
sim_8 <- simulateOutliers(100,0.3,200,0.75,40)
sim_9 <- simulateOutliers(100,0.5,200,0.75,40)
sim_10 <- simulateOutliers(100,0.1,200,0.5,1)
sim_11 <- simulateOutliers(100,0.3,200,0.5,1)
sim_12 <- simulateOutliers(100,0.5,200,0.5,1)
sim_13 <- simulateOutliers(100,0.1,200,0.5,5)
sim_14 <- simulateOutliers(100,0.3,200,0.5,5)
sim_15 <- simulateOutliers(100,0.5,200,0.5,5)
sim_16 <- simulateOutliers(100,0.1,200,0.5,40)
sim_17 <- simulateOutliers(100,0.3,200,0.5,40)
sim_18 <- simulateOutliers(100,0.5,200,0.5,40)

#construct dataframes (per type of outlier)
Clean <-
  rbind(
    sim_1$Clean[(length(sim_10$Clean) - 5):length(sim_10$Clean)],
    sim_4$Clean[(length(sim_10$Clean) - 5):length(sim_10$Clean)],
    sim_7$Clean[(length(sim_10$Clean) - 5):length(sim_10$Clean)],
    sim_10$Clean[(length(sim_10$Clean) - 5):length(sim_10$Clean)],
    sim_13$Clean[(length(sim_10$Clean) - 5):length(sim_10$Clean)],
    sim_16$Clean[(length(sim_10$Clean) - 5):length(sim_10$Clean)]
  )
rownames(Clean) <- c("alpha=0.75 - p=1","alpha=0.75 - p=5","alpha=0.75 - p=40","alpha=0.5 - p=1","alpha=0.5 - p=5","alpha=0.5 - p=40")
colnames(Clean) <- c("RMSE - OLS", "RMSE - LTS", "RMSE - PR", "RMSPE - OLS", "RMSPE - LTS", "RMSPE - PR")

GLP <- rbind(
  sim_1$GLP[(length(sim_1$GLP) - 5):length(sim_1$GLP)],
  sim_2$GLP[(length(sim_2$GLP) - 5):length(sim_2$GLP)],
  sim_3$GLP[(length(sim_3$GLP) - 5):length(sim_3$GLP)],
  sim_4$GLP[(length(sim_4$GLP) - 5):length(sim_4$GLP)],
  sim_5$GLP[(length(sim_5$GLP) - 5):length(sim_5$GLP)],
  sim_6$GLP[(length(sim_6$GLP) - 5):length(sim_6$GLP)],
  sim_7$GLP[(length(sim_7$GLP) - 5):length(sim_7$GLP)],
  sim_8$GLP[(length(sim_8$GLP) - 5):length(sim_8$GLP)],
  sim_9$GLP[(length(sim_9$GLP) - 5):length(sim_9$GLP)],
  sim_10$GLP[(length(sim_10$GLP) - 5):length(sim_10$GLP)],
  sim_11$GLP[(length(sim_11$GLP) - 5):length(sim_11$GLP)],
  sim_12$GLP[(length(sim_12$GLP) - 5):length(sim_12$GLP)],
  sim_13$GLP[(length(sim_13$GLP) - 5):length(sim_13$GLP)],
  sim_14$GLP[(length(sim_14$GLP) - 5):length(sim_14$GLP)],
  sim_15$GLP[(length(sim_15$GLP) - 5):length(sim_15$GLP)],
  sim_16$GLP[(length(sim_16$GLP) - 5):length(sim_16$GLP)],
  sim_17$GLP[(length(sim_17$GLP) - 5):length(sim_17$GLP)],
  sim_18$GLP[(length(sim_18$GLP) - 5):length(sim_18$GLP)]
)
rownames(GLP) <- c("0.75, 1, 0.1","0.75, 1, 0.3", "0.75, 1, 0.5","0.75, 5, 0.1","0.75, 5, 0.3", "0.75, 5, 0.5","0.75, 40, 0.1","0.75, 40, 0.3", "0.75, 40, 0.5",
                   "0.5, 1, 0.1","0.5, 1, 0.3", "0.5, 1, 0.5","0.5, 5, 0.1","0.5, 5, 0.3", "0.5, 5, 0.5","0.5, 40, 0.1","0.5, 40, 0.3", "0.5, 40, 0.5")
colnames(GLP) <- c("RMSE - OLS", "RMSE - LTS", "RMSE - PR", "RMSPE - OLS", "RMSPE - LTS", "RMSPE - PR")

BLP <- rbind(
  sim_1$BLP[(length(sim_1$BLP) - 5):length(sim_1$BLP)],
  sim_2$BLP[(length(sim_2$BLP) - 5):length(sim_2$BLP)],
  sim_3$BLP[(length(sim_3$BLP) - 5):length(sim_3$BLP)],
  sim_4$BLP[(length(sim_4$BLP) - 5):length(sim_4$BLP)],
  sim_5$BLP[(length(sim_5$BLP) - 5):length(sim_5$BLP)],
  sim_6$BLP[(length(sim_6$BLP) - 5):length(sim_6$BLP)],
  sim_7$BLP[(length(sim_7$BLP) - 5):length(sim_7$BLP)],
  sim_8$BLP[(length(sim_8$BLP) - 5):length(sim_8$BLP)],
  sim_9$BLP[(length(sim_9$BLP) - 5):length(sim_9$BLP)],
  sim_10$BLP[(length(sim_10$BLP) - 5):length(sim_10$BLP)],
  sim_11$BLP[(length(sim_11$BLP) - 5):length(sim_11$BLP)],
  sim_12$BLP[(length(sim_12$BLP) - 5):length(sim_12$BLP)],
  sim_13$BLP[(length(sim_13$BLP) - 5):length(sim_13$BLP)],
  sim_14$BLP[(length(sim_14$BLP) - 5):length(sim_14$BLP)],
  sim_15$BLP[(length(sim_15$BLP) - 5):length(sim_15$BLP)],
  sim_16$BLP[(length(sim_16$BLP) - 5):length(sim_16$BLP)],
  sim_17$BLP[(length(sim_17$BLP) - 5):length(sim_17$BLP)],
  sim_18$BLP[(length(sim_18$BLP) - 5):length(sim_18$BLP)]
)
rownames(BLP) <- c("0.75, 1, 0.1","0.75, 1, 0.3", "0.75, 1, 0.5","0.75, 5, 0.1","0.75, 5, 0.3", "0.75, 5, 0.5","0.75, 40, 0.1","0.75, 40, 0.3", "0.75, 40, 0.5",
                   "0.5, 1, 0.1","0.5, 1, 0.3", "0.5, 1, 0.5","0.5, 5, 0.1","0.5, 5, 0.3", "0.5, 5, 0.5","0.5, 40, 0.1","0.5, 40, 0.3", "0.5, 40, 0.5")
colnames(BLP) <- c("RMSE - OLS", "RMSE - LTS", "RMSE - PR", "RMSPE - OLS", "RMSPE - LTS", "RMSPE - PR")

VO <- rbind(
  sim_1$VO[(length(sim_1$VO) - 5):length(sim_1$VO)],
  sim_2$VO[(length(sim_2$VO) - 5):length(sim_2$VO)],
  sim_3$VO[(length(sim_3$VO) - 5):length(sim_3$VO)],
  sim_4$VO[(length(sim_4$VO) - 5):length(sim_4$VO)],
  sim_5$VO[(length(sim_5$VO) - 5):length(sim_5$VO)],
  sim_6$VO[(length(sim_6$VO) - 5):length(sim_6$VO)],
  sim_7$VO[(length(sim_7$VO) - 5):length(sim_7$VO)],
  sim_8$VO[(length(sim_8$VO) - 5):length(sim_8$VO)],
  sim_9$VO[(length(sim_9$VO) - 5):length(sim_9$VO)],
  sim_10$VO[(length(sim_10$VO) - 5):length(sim_10$VO)],
  sim_11$VO[(length(sim_11$VO) - 5):length(sim_11$VO)],
  sim_12$VO[(length(sim_12$VO) - 5):length(sim_12$VO)],
  sim_13$VO[(length(sim_13$VO) - 5):length(sim_13$VO)],
  sim_14$VO[(length(sim_14$VO) - 5):length(sim_14$VO)],
  sim_15$VO[(length(sim_15$VO) - 5):length(sim_15$VO)],
  sim_16$VO[(length(sim_16$VO) - 5):length(sim_16$VO)],
  sim_17$VO[(length(sim_17$VO) - 5):length(sim_17$VO)],
  sim_18$VO[(length(sim_18$VO) - 5):length(sim_18$VO)]
)
rownames(VO) <- c("0.75, 1, 0.1","0.75, 1, 0.3", "0.75, 1, 0.5","0.75, 5, 0.1","0.75, 5, 0.3", "0.75, 5, 0.5","0.75, 40, 0.1","0.75, 40, 0.3", "0.75, 40, 0.5",
                   "0.5, 1, 0.1","0.5, 1, 0.3", "0.5, 1, 0.5","0.5, 5, 0.1","0.5, 5, 0.3", "0.5, 5, 0.5","0.5, 40, 0.1","0.5, 40, 0.3", "0.5, 40, 0.5")
colnames(VO) <- c("RMSE - OLS", "RMSE - LTS", "RMSE - PR", "RMSPE - OLS", "RMSPE - LTS", "RMSPE - PR")

allOutliers <-
  rbind(
    sim_1$allOutliers[(length(sim_1$allOutliers) - 5):length(sim_1$allOutliers)],
    sim_2$allOutliers[(length(sim_2$allOutliers) - 5):length(sim_2$allOutliers)],
    sim_3$allOutliers[(length(sim_3$allOutliers) - 5):length(sim_3$allOutliers)],
    sim_4$allOutliers[(length(sim_4$allOutliers) - 5):length(sim_4$allOutliers)],
    sim_5$allOutliers[(length(sim_5$allOutliers) - 5):length(sim_5$allOutliers)],
    sim_6$allOutliers[(length(sim_6$allOutliers) - 5):length(sim_6$allOutliers)],
    sim_7$allOutliers[(length(sim_7$allOutliers) - 5):length(sim_7$allOutliers)],
    sim_8$allOutliers[(length(sim_8$allOutliers) - 5):length(sim_8$allOutliers)],
    sim_9$allOutliers[(length(sim_9$allOutliers) - 5):length(sim_9$allOutliers)],
    sim_10$allOutliers[(length(sim_10$allOutliers) - 5):length(sim_10$allOutliers)],
    sim_11$allOutliers[(length(sim_11$allOutliers) - 5):length(sim_11$allOutliers)],
    sim_12$allOutliers[(length(sim_12$allOutliers) - 5):length(sim_12$allOutliers)],
    sim_13$allOutliers[(length(sim_13$allOutliers) - 5):length(sim_13$allOutliers)],
    sim_14$allOutliers[(length(sim_14$allOutliers) - 5):length(sim_14$allOutliers)],
    sim_15$allOutliers[(length(sim_15$allOutliers) - 5):length(sim_15$allOutliers)],
    sim_16$allOutliers[(length(sim_16$allOutliers) - 5):length(sim_16$allOutliers)],
    sim_17$allOutliers[(length(sim_17$allOutliers) - 5):length(sim_17$allOutliers)],
    sim_18$allOutliers[(length(sim_18$allOutliers) - 5):length(sim_18$allOutliers)]
  )
rownames(allOutliers) <- c("0.75, 1, 0.1","0.75, 1, 0.3", "0.75, 1, 0.5","0.75, 5, 0.1","0.75, 5, 0.3", "0.75, 5, 0.5","0.75, 40, 0.1","0.75, 40, 0.3", "0.75, 40, 0.5",
                   "0.5, 1, 0.1","0.5, 1, 0.3", "0.5, 1, 0.5","0.5, 5, 0.1","0.5, 5, 0.3", "0.5, 5, 0.5","0.5, 40, 0.1","0.5, 40, 0.3", "0.5, 40, 0.5")
colnames(allOutliers) <- c("RMSE - OLS", "RMSE - LTS", "RMSE - PR", "RMSPE - OLS", "RMSPE - LTS", "RMSPE - PR")

#### plots ####
#1 - GLP
# cols1 <- c("OLS"="blue","Plug-in"="red","LTS"="blue")
# data <- df1
# plugin <- lmDetMCD(data[,2],data[,1],alpha=0.5)
# lts <- ltsReg(data[,2],data[,1],alpha=0.5)
# ols <- lm(data[,1]~data[,2])
# predicted_data_plugin <- data.frame(X=data[,2],pred_Y = plugin$fitted.values)
# predicted_data_lts <- data.frame(X=data[,2],pred_Y=lts$fitted.values)
# predicted_data_ols <- data.frame(X=data[,2],pred_Y = ols$fitted.values)
# p1 <- ggplot(NULL) +
#   geom_point(data = data, aes(x=X, y=Y)) +
#   geom_line(data = predicted_data_plugin,aes(x=X,y=pred_Y,color="Plug-in")) +
#   geom_line(data = predicted_data_ols,aes(x=X,y=pred_Y,color="OLS")) +
#   geom_line(data = predicted_data_lts,aes(x=X,y=pred_Y,color="LTS")) +
#   ggtitle("Good Leverage Points, a=0.5") +
#   labs( color ="Legend",
#         x ="X",
#         y = "Y") + 
#   scale_colour_manual(values = cols1) +
#   theme(legend.position = c(0.92,0.92),legend.key = element_blank())
# p1

#2 - BLP
# cols2 <- c("OLS"="green","Plug-in"="blue","LTS"="blue")
# data <- df2
# plugin <- lmDetMCD(data[,2],data[,1],alpha=0.5)
# lts <- ltsReg(data[,2],data[,1],alpha=0.5)
# ols <- lm(data[,1]~data[,2])
# predicted_data_plugin <- data.frame(X=data[,2],pred_Y = plugin$fitted.values)
# predicted_data_lts <- data.frame(X=data[,2],pred_Y=lts$fitted.values)
# predicted_data_ols <- data.frame(X=data[,2],pred_Y = ols$fitted.values)
# p2 <- ggplot(NULL) +
#   geom_point(data = data, aes(x=X, y=Y)) +
#   geom_line(data = predicted_data_plugin,aes(x=X,y=pred_Y,color="Plug-in")) +
#   geom_line(data = predicted_data_ols,aes(x=X,y=pred_Y,color="OLS")) +
#   geom_line(data = predicted_data_lts,aes(x=X,y=pred_Y,color="LTS")) +
#   ggtitle("Bad Leverage Points, a=0.5") +
#   labs( color ="Legend",
#         x ="X",
#         y = "Y") + 
#   scale_colour_manual(values = cols2) +
#   theme(legend.position = c(0.92,0.92),legend.key = element_blank())
# p2

#3 - VO
# cols3 <- c("OLS"="green","Plug-in0.5"="blue","LTS0.5"="blue","Plug-in0.75"="orange","LTS0.75"="black")
# data <- df3
# plugin <- lmDetMCD(data[,2],data[,1],alpha=0.5)
# lts <- ltsReg(data[,2],data[,1],alpha=0.5)
# plugin2 <- lmDetMCD(data[,2],data[,1],alpha=0.75)
# lts2 <- ltsReg(data[,2],data[,1],alpha=0.75)
# ols <- lm(data[,1]~data[,2])
# predicted_data_plugin <- data.frame(X=data[,2],pred_Y = plugin$fitted.values)
# predicted_data_lts <- data.frame(X=data[,2],pred_Y=lts$fitted.values)
# predicted_data_plugin2 <- data.frame(X=data[,2],pred_Y = plugin2$fitted.values)
# predicted_data_lts2 <- data.frame(X=data[,2],pred_Y=lts2$fitted.values)
# predicted_data_ols <- data.frame(X=data[,2],pred_Y = ols$fitted.values)
# p3 <- ggplot(data, aes(x=X, y=Y)) +
#   geom_point() +
#   geom_line(data = predicted_data_plugin,aes(x=X,y=pred_Y,color="Plug-in0.5")) +
#   geom_line(data = predicted_data_ols,aes(x=X,y=pred_Y,color="OLS")) +
#   geom_line(data = predicted_data_lts,aes(x=X,y=pred_Y,color="LTS0.5")) +
#   geom_line(data = predicted_data_plugin2,aes(x=X,y=pred_Y,color="Plug-in0.75")) +
#   geom_line(data = predicted_data_lts2,aes(x=X,y=pred_Y,color="LTS0.75")) +
#   ggtitle("Vertical Outliers") +
#   labs( color ="Legend",
#         x ="X",
#         y = "Y") + 
#   scale_colour_manual(values = cols3) +
#   theme(legend.position = c(0.92,0.92),legend.key = element_blank())
# p3

#4 - ALL
# cols4 <- c("OLS"="green","Plug-in"="red","LTS"="blue")
# data <- df4
# plugin <- lmDetMCD(data[,2],data[,1],alpha=0.75)
# lts <- ltsReg(data[,2],data[,1],alpha=0.75)
# ols <- lm(data[,1]~data[,2])
# predicted_data_plugin <- data.frame(X=data[,2],pred_Y = plugin$fitted.values)
# predicted_data_lts <- data.frame(X=data[,2],pred_Y=lts$fitted.values)
# predicted_data_ols <- data.frame(X=data[,2],pred_Y = ols$fitted.values)
# p4 <- ggplot(data, aes(x=X, y=Y)) +
#   geom_point() +
#   geom_line(data = predicted_data_plugin,aes(x=X,y=pred_Y,color="Plug-in")) +
#   geom_line(data = predicted_data_ols,aes(x=X,y=pred_Y,color="OLS")) +
#   geom_line(data = predicted_data_lts,aes(x=X,y=pred_Y,color="LTS")) +
#   ggtitle("All types of outliers, a=0.75") +
#   labs( color ="Legend",
#         x ="X",
#         y = "Y") + 
#   scale_colour_manual(values = cols4) +
#   theme(legend.position = c(0.95,0.92),legend.key = element_blank())
# p4

#plot together in 1 window
# grid.arrange(p1,p2,p3,p4,nrow=2)





