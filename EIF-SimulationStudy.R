# --------------------------------------------------------------------
# Author:
# Group 14: Ana Tasia Bueno de Mesquita (424406), Rosa ten Kate (413775), Frederique Ram (431623), Tiemen de Jong (430839)
# --------------------------------------------------------------------
#The library 
library(mvtnorm)

#Calculate the regression coefficients based on the covariance matrix
Regression<- function(S,mu)
{
  B <-  solve(S[-1,-1])%*%S[-1,1]
  a <- mu[1] - crossprod(mu[-1], B) 
  Result <- list("int."=a,"slope"=B)
  return(Result)
}

#Calculate the Root Mean Squared (Prediction) Error
RMSE <- function (true, estimated) {
  MSE <- mean((true-estimated)^2)
  RMSE <- sqrt(MSE)
  return (RMSE)
}

#Emperical influence function in the x-direction 
# Input: 
# x ....... well-behaved data set 
# R ....... range of x-values which will be replaced
# type .... type of estimation method "OLS", "LTS", "PR"

# Output: 
# Slope ... values of the EIF for the slope coefficient
# Int ..... values of the EIF for the intercept coefficient
# R ....... range of x-values which are replaced
# K ....... ID of the observation which is randomly replaced
EIF <- function(x,R,type)
{ X.r <- x
n <- nrow(x)
p <- ncol(x)
#random observation i to be replaced
K <- sample(1:n,1) 
#initialize vectors 
B.r <- 0
a.r <- 0
EIF.B <- 0 
EIF.A <-0
if (type == "OLS"){
  #### intial values slope and intercept OLS  ####
  mu <- apply(x,2,mean)
  B <- Regression(cov(x),mu)$slope
  a <- Regression(cov(x),mu)$int.
  #randomly replace one observation 
  for (i in 1:length(R))  {
    X.r <- x
    X.r[K,2] <- R[i] #replace X
    # OLS for the replaced data set
    mu.r <- apply(X.r,2,mean)
    B.r[i] <-  Regression(cov(X.r),mu.r)$slope
    a.r[i] <-  Regression(cov(X.r),mu)$int
    EIF.B[i]<-  n*(B.r[i]-B)
    EIF.A[i]<-  n*(a.r[i]-a)
  }}
if (type == "LTS"){#### intial values slope and intercept LTS ####
  lts <- ltsReg(x[,1]~x[,2])
  B<- lts$coefficients[[2]]
  a<- lts$coefficients[[1]]
  #### Replace a random observation with a values in R####  
  for (i in 1:length(R)) {
    X.r <- x
    X.r[K,2] <- R[i]
    # LTS for the replaced dataset
    B.r[i]<- ltsReg(X.r[,1]~X.r[,2])$coefficients[[2]]
    a.r[i]<- ltsReg(X.r[,1]~X.r[,2])$coefficients[[1]]
    EIF.B[i]<-  n*(B.r[i]-B)
    EIF.A[i]<-  n*(a.r[i]-a)}
}
if (type == "PLUG")
{ alpha <- as.double(readline(prompt = "Enter alpha:"))
#### intial values slope and intercept LTS ####
MCD <- covDetMCD(x,alpha)
S <- MCD$cov
C <- MCD$center
b <- solve(S[2:p,2:p])%*%S[1,2:p]
a <- C[1] - crossprod(C[-1], b)
for (i in 1:length(R)) {
  X.r <- x
  X.r[K,2] <- R[i]
  #CovdetMCD for the replaced dataset
  MCD <- covDetMCD(X.r,alpha)
  S <- MCD$cov
  C <- MCD$center
  B.r[i]<- solve(S[2:p,2:p])%*%S[1,2:p]
  a.r[i]<- C[1] - crossprod(C[-1], b)
  EIF.B[i]<-  n*(B.r[i]- b)
  EIF.A[i]<-  n*(a.r[i]- a )
}

}
resutls <- list("slope" = EIF.B,"int." = EIF.A,"Replace" = R, "obs." = K)
return(resutls)
}
# calculate the EIF in the x-direction for every type of estimation method
EIF.ols <- as.data.frame(EIF(x, -5:5, "OLS" ))
EIF.lts <-as.data.frame(EIF(x, -5:5, "LTS"))
EIF.plug <- as.data.frame(EIF(x, -5:5, "PLUG"))


#Emperical influence function in the (x,y)-direction 
# Input: 
# x ....... well-behaved data set 
# R_x ....... range of x-values which will be replaced
# R_y ....... range of x-values which will be replaced
# type .... type of estimation method "OLS", "LTS", "PR"

# Output: 
# Slope ... values of the EIF for the slope coefficient
# Int ..... values of the EIF for the intercept coefficient
EIF.square <- function(x,R_y,R_x,type)
{
  n <- nrow(x)
  p <- ncol(x)
  #random observation i to be replaced
  #K <- 50 #sample(1:n,1) 
  #initialize vectors 
  B.r <- matrix(0, nrow = length(R_y),ncol=length(R_x))
  a.r <- matrix(0, nrow = length(R_y),ncol=length(R_x) )
  EIF.B <- matrix(0, nrow = length(R_y),ncol=length(R_x) )
  EIF.A <-matrix(0, nrow = length(R_y),ncol=length(R_x) )
  if (type == "OLS"){
    #### intial values slope and intercept OLS  ####
    mu <- apply(x,2,mean)
    B <- Regression(cov(x),mu)$slope
    a <- Regression(cov(x),mu)$int.
    #randomly replace one observation 
    for (i in 1:length(R_y)) {
      for(j in 1:length(R_x))
      {
        X.r <- rbind(x,c(R_y[i],R_x[j]))
        #X.r[,1] <- R[i]
        #X.r[K,2] <- R[j]
        # OLS for the replaced data set
        mu.r <- apply(X.r,2,mean)
        B.r[i,j] <-  Regression(cov(X.r),mu.r)$slope
        a.r[i,j] <-  Regression(cov(X.r),mu.r)$int
        EIF.B[i,j]<-  (n+1)*(B.r[i,j]-B)
        EIF.A[i,j]<-  (n+1)*(a.r[i,j]-a)
      }}
  }  
  if (type == "LTS"){#### intial values slope and intercept LTS ####
    lts <- ltsReg(x[,1]~x[,2])
    B<- lts$coefficients[[2]]
    a<- lts$coefficients[[1]]
    #### Replace a random observation with a values in R####  
    for (i in 1:length(R_y)) {
      for(j in 1:length(R_x)){
        X.r <- rbind(x,c(R_y[i],R_x[j]))
        # LTS for the replaced dataset
        B.r[i,j]<- ltsReg(X.r[,1]~X.r[,2])$coefficients[[2]]
        a.r[i,j]<- ltsReg(X.r[,1]~X.r[,2])$coefficients[[1]]
        EIF.B[i,j]<-  (n+1)*(B.r[i,j]-B)
        EIF.A[i,j]<-  (n+1)*(a.r[i,j]-a)}
    }}
  
  if (type == "PLUG")
  { alpha <- as.double(readline(prompt = "Enter alpha:"))
  #### intial values slope and intercept LTS ####
  MCD <- covDetMCD(x,alpha)
  S <- MCD$cov
  C <- MCD$center
  b <- solve(S[2:p,2:p])%*%S[1,2:p]
  a <- C[1] - crossprod(C[-1], b)
  for (i in 1:length(R_y)) {
    for(j in 1:length(R_x)){
      X.r <- rbind(x,c(R_y[i],R_x[j]))
      #CovdetMCD for the replaced dataset
      MCD <- covDetMCD(X.r,alpha)
      S <- MCD$cov
      C <- MCD$center
      B.r[i,j]<- solve(S[2:p,2:p])%*%S[1,2:p]
      a.r[i,j]<- C[1] - crossprod(C[-1], b)
      EIF.B[i,j]<-  (n+1)*(B.r[i,j]- b)
      EIF.A[i,j]<-  (n+1)*(a.r[i,j]- a )
    }}
  }
  
  resutls <- list("slope" = EIF.B,"int." = EIF.A)
  return(resutls)
}

# calculate the EIF in the (x,y)-direction for every type of estimation method
EIF.ols <- EIF.square(data,6:20,6:38, "OLS" )
EIF.lts <- EIF.square(data, 6:20,6:38, "LTS")
EIF.plug <- EIF.square(data,6:20,6:38, "PLUG")


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
# Please instal to run the code: 
#library(ggplot2) 
#library(plotly)

# Slope x-direciton  
#df <- as.data.frame(x)
#Slope <- ggplot(df, aes(x=V1, y=V2))  + 
#  geom_point(color="darkgrey") +
#  geom_line(data = EIF.ols,aes(x=Replace,y=slope,color = "OLS")) +
# geom_line(data = EIF.lts,aes(x=Replace,y=slope, color = "LTS")) +
#  geom_line(data = EIF.plug,aes(x=Replace,y=slope, color = "Plug-in")) + 
#  geom_hline(yintercept = 0, linetype="dashed", color = "black") + xlab("X") + ylab("Emperical influence function")+
#  geom_point(data=df[36,], color="black") + 
#  labs( color ="",
#        x ="X",
#        y = "Emperical influence fucntion") + 
#  scale_colour_manual(values = cols) +
#  theme(legend.position = c(0.92,0.92),legend.key = element_blank())
#Intercept x-direciton
#Int <- ggplot(df, aes(x=V1, y=V2))  + 
#  geom_point(color="darkgrey") +
# geom_line(data = EIF.ols,aes(x=Replace,y=int.,color = "OLS")) +
#  geom_line(data = EIF.lts,aes(x=Replace,y=int., color = "LTS")) +
#  geom_line(data = EIF.plug,aes(x=Replace,y=int., color = "Plug-in")) + 
#  geom_hline(yintercept = 0, linetype="dashed", color = "black") + xlab("X") + ylab("Emperical influence function")+
#  geom_point(data=df[36,], color="black") + 
#  labs( color ="",
#       x ="X",
#        y = "Emperical influence fucntion") + 
#  scale_colour_manual(values = cols) +
#  theme(legend.position = c(0.92,0.92),legend.key = element_blank())

#3D plots (x,y)-direction 
#slope
#ols <- plot_ly(x =  6:38, y =6:20, z = ~EIF.ols$slope) %>% add_surface() %>%
#  layout(title ="Slope - OLS",
#         scene = list(
#           xaxis = list(title = "x"),
#           yaxis = list(title = "y"),
#           zaxis = list(title = "EIF")
#         ))
#lts<- plot_ly(x = 6:38, y =6:20, z = ~EIF.lts$slope) %>% add_surface()%>%
#  layout(title ="Slope - LTS",
#         scene = list(
#           xaxis = list(title = "x"),
#           yaxis = list(title = "y"),
#           zaxis = list(title = "EIF")
#         ))
#
#plug<- plot_ly(x = 6:38, y =6:20, z =EIF.plug$slope) %>% add_surface()%>%
#  layout(title ="Slope - Plug-in",
#         scene = list(
#           xaxis = list(title = "x"),
#           yaxis = list(title = "y"),
#           zaxis = list(title = "EIF")
#         ))

#intercept
#ols.int <- plot_ly(x =  6:38, y =6:20, z = ~EIF.ols$int.) %>% add_surface() %>%
#  layout(title ="Intercept- OLS",
#         scene = list(
#           xaxis = list(title = "x"),
#           yaxis = list(title = "y"),
#           zaxis = list(title = "EIF")
#         ))
#lts.int<- plot_ly(x = 6:38, y =6:20, z = ~EIF.lts$int.) %>% add_surface()%>%
#  layout(title ="Intercept- LTS",
#         scene = list(
#           xaxis = list(title = "x"),
#           yaxis = list(title = "y"),
#           zaxis = list(title = "EIF")
#         ))

#plug.int<- plot_ly(x = 6:38 ,y=6:20, z =~EIF.plug$int.) %>% add_surface()%>%
#  layout(title ="Intercept - Plug-in",
#         scene = list(
#           xaxis = list(title = "x"),
#           yaxis = list(title = "y"),
#           zaxis = list(title = "EIF")
#         ))


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
