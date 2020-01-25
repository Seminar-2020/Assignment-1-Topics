# --------------------------------------------------------------------
# Author:
# Group 14: Tasia Bueno de Mesquita (), Rosa ten Kate (), Frederique Ram (431623), Tiemen de Jong ()
# --------------------------------------------------------------------

## Use this code skeleton to implement the deterministic MCD and the plug-in
## robust regression estimator.  Please submit your implementation together
## with your report.

## IMPORTANT: Please do not change any function names and make sure that you
##            use the correct input and output as specified for each function.
##            This will simplify grading because each group's code will be
##            similarly structured.  However, feel free to add other arguments
##            to the function definitions as needed.



## Functions for initial estimators

# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent
corHT <- function(z) {
  # *enter your code here*
}

# spearman correlation matrix
corSpearman <- function(z) {
  # *enter your code here*
}

# correlation matrix based on normal scores of the ranks
corNSR <- function(z) {
  # *enter your code here*
}

# modified spatial sign covariance matrix
covMSS <- function(z) {
  # *enter your code here*
}

# covariance matrix based on first step of BACON
covBACON1 <- function(z) {
  # *enter your code here*
}

# raw OGK estimator of the covariance matrix with median and Qn
rawCovOGK <- function(z) {
  # *enter your code here*
  # Hint: have a look at function covOGK() in package robustbase
}



## Main function for deterministic MCD algorithm

# Input:
# x ....... data matrix
# alpha ... proportion of observations to be used for subset size
# anything else you need

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# any other output you want to return

covDetMCD <- function(x, alpha, ...) {
  n <- nrow(x)
  h <- h.alpha.n(alpha,n,ncol(x))
  
  # standardize data
  rob_std <- function(col) { #not sure if Qn() is allowed, estimate for variation in a vector
    std <- (col-median(col))/Qn(col)
    return(std)
  }
  z <- as.data.frame(apply(x, 2, "rob_std"))
  
  # obtain 6 initial estimates
  S1 <- corHT(z)
  S2 <- corSpearman(z)
  S3 <- corNSR(z)
  S4 <- covMSS(z)
  S5 <- covBACON1(z)
  S6 <- rawCovOGK(z)
  Sk <- list(S1,S2,S3,S4,S5,S6)
  
  proper_ev <- function(z,S) { #z cannot be a dataframe but has to be a matrix!
    E <- eigen(S)$vectors
    B <- z%*%E
    Qn2 <- apply(z,2,Qn)^2
    L <- diag(Qn2)
    Sigma_hat <- E%*%L%*%t(E)
    mu_hat <- Sigma_hat^(1/2)%*%apply(z%*%Sigma_hat^(-1/2), 2, median) #not 100% sure dat dit de goede median is, returnt nu NA maar testen als een vd functies geschreven is
    result <- list("center"= mu_hat, "scatter" = Sigma_hat)
    return(result)
  }
  
  results_mu <- list()
  results_Sigma <- list()
  for (i in 1:6) {
    result <- proper_ev(z,Sk[[i]])
    results_mu[[i]] <- result[[1]]
    results_Sigma[[i]] <- result[[2]] #xoxox
  }
  
  algorithm <- function(z, mu_hat, Sigma_hat) {
    z$distances <- mahalanobis(z, mu_hat, Sigma_hat)^(1/2)
    h0 <- ceiling(nrow(z)/2)
    h <- h.alpha.n(alpha,nrow(z),ncol(z))
    z$rank <- rank(z$distances, ties.method="random")
    initial_subset <- z[z$rank<=h0,1:(ncol(z)-2)] #pick h0 smallest distances as initial subset and delete columns for distance & rank
    T_H0 <-colSums(initial_subset)/h0
    T_H0_rep <- matrix(T_H0,nrow=h,ncol=ncol(initial_subset),byrow=TRUE)
    S_H0 <- (initial_subset-T_H0_rep)*t(initial_subset-T_H_rep)/h0
    for (i in 1:1000) {
      # calculate mahalanobis distance for all data points given estimated mean and cov
      z$distances <- mahalanobis(z, T_H0, S_H0)^(1/2)
      # order all mahalanobis distances
      z$rank <- rank(z$distances, ties.method="random")
      # pick h smallest, i.e. set weights of h smallest to 1
      z$weights <- ifelse(z$distances<=h, 1, 0)
      # calculate new mean and cov based on this subset
      T_H <-colSums(z[z$weights==1,1:(ncol(z)-3)])/h #don't use distances, rank & weights to calculate
      T_H_rep <- matrix(T_H,nrow=h,ncol=ncol(z)-3,byrow=TRUE)
      S_H <- (z[z$weights==1,1:(ncol(z)-3)]-T_H_rep)*t(z[z$weights==1,1:(ncol(z)-3)]-T_H_rep)/h
      if (det(S_H)==det(S_H0)) { # stopping criterium
        break
      } 
      T_H0 <- T_H
      S_H0 <- S_H
      # repeat
    }
    result <- list("raw.center"=T_H0, "raw.cov"=S_H0,"weights"=z$weights)
    return(result)
  }
  
  results_raw_center <- list()
  results_raw_cov <- list()
  results_weights <- list()
  
  for (i in 1:6) {
    result <- algorithm(z,results_mu[[i]],results_Sigma[[i]])
    results_raw_center[[i]] <- result[[1]]
    results_raw_cov[[i]] <- result[[2]]
    results_weights[[i]] <- result[[3]]
  }
  
  best_det <- function(benchmark){
    best_det <- 1
    for (i in 2:6) {
      if (det(results_raw_cov[[i]])<benchmark) {
        best_det <- i
        benchmark <- det(results_raw_cov[[i]])
      }
    }
    return(c(i,benchmark))
  }
  det_1 <- det(results_raw_cov[[1]])
  ind_det <- best_det(det_1) 
  raw.cov <- results_raw_cov[[ind_det[1]]]
  raw.center <- results_raw_center[[ind_det[1]]]
  #dit werkt nog niet helemaal
  best$weights <- results_weights[[ind_det[1]]]; #vector with ones and zeros
  best$indices <- seq(1:nrow(z))
  best <- best[weights==1,2]
  
  # Reweighting step to transform and obtain final estimate
  Q <- qchisq(1-0.025, ncol(z), ncp = 0, lower.tail = TRUE, log.p = FALSE)
  z$weights_r <- ifelse(mahalanobis(z, raw.center, raw.cov)<=Q, 1, 0)
  weights <- z$weights_r
  center <-colSums(z[z$weights_r==1,1:(ncol(z)-3)])/sum(z$weights_r) #don't use distances, rank & weights to calculate
  center_rep <- matrix(center,nrow=sum(z$weights_r),ncol=ncol(z)-3,byrow=TRUE)
  cov <- (z[z$weights_r==1,1:(ncol(z)-3)]-center_rep)*t(z[z$weights_r==1,1:(ncol(z)-3)]-center_rep)/sum(z$weights_r)
  
  # Please note that the subset sizes for the MCD are not simply fractions of 
  # the number of observations in the data set, as discussed in the lectures.
  # You can use function h.alpha.n() from package robustbase to compute the 
  # subset size.
}



## Function for regression based on the deterministic MCD

# Input:
# x ........ matrix of explanatory variables
# y ........ response variable
# alpha .... proportion of observations to be used for the subset size in the 
#            MCD estimator
# anything else you need

# Output
# A list with the following components:
# coefficients .... estimated regression coefficients based on the reweighted
#                   deterministic MCD covariance matrix
# fitted.values ... fitted values for all observations in the data
# residuals ....... residuals for all observations in the data
# MCD ............. R object for the deterministic MCD (entire output from
#                   function covDetMCD())
# any other output you want to return

lmDetMCD <- function(x, y, alpha, ...) {
  # *enter your code here*
}
