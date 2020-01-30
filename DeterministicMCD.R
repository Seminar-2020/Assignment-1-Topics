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

library(robustbase)

#euclidean distance function 
euc.dist<- function(x1) sqrt(sum((x1) ^ 2))

#### Functions for initial estimators ####
# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent 
corHT <- function(z) {
  y  = tanh(z)
  s.1 = cor(y, method = "pearson")
  return(s.1)
}

# spearman correlation matrix s.2
corSpearman <- function(z) {
  r <- apply(z, 2, rank)
  s.2 = cor(r, method = "pearson")
  return(s.2)
}

# correlation matrix based on normal scores of the ranks s.3
corNSR <- function(z) {
  #calculate the rank of each colum
  r <- apply(z, 2, rank)
  #normalize the ranks 
  T = qnorm(((r-c(1/3,1/3))/(nrow(z) + 1/3))) #inv of normal distribution
  #calculate the correlation 
  s.3 = cor(T, method = "pearson")
  return(s.3)
}

# modified spatial sign covariance matrix s.4
covMSS <- function(z) {
  z <- as.matrix(z)
  d  <-  apply(z,1, euc.dist)
  k  <- z/d
  k[is.na(k)] <- 0  #comment, paper sign, if function
  s.4 <- 1/nrow(z)*crossprod(k)
  return(s.4)
} 

# covariance matrix based on first step of BACON #Juist variables, beter schrijven s.5
covBACON <- function(z) {
  #Based on the norm 
  dist <- apply(z, 1, euc.dist)
  zz <- cbind(z,dist)
  zz <- zz[order(dist),]
  z.s <- zz[1:ceiling(nrow(z)/2), 1:ncol(z)]
  s.5 <- cov(z.s)
  return (s.5)
} 

# raw OGK estimator of the covariance matrix with median and Qn s.6
rawCovOGK <- function(z) {
  raw<-covOGK(z,n.iter=2,sigmamu = s_Qn, rcov = covGK, weight.fn = hard.rejection)
  s.6 <- raw$cov
  return(s.6)
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
  p <- ncol(x)
  h <- h.alpha.n(alpha,n,p)
  
  rob_paper <- function(x) {
    Var <- apply(x, 2, mad)
    medians <- apply(x,2,median)
    A <- diag(Var^(-1))
    ones <- matrix(1, nrow=n, ncol=1)
    v <- medians/Var
    z <- as.matrix(x)%*%A-ones%*%v
    # evt nog kolomnamen herstellen
    return(z)
  }
  z <- as.data.frame(rob_paper(x))
  
  # obtain 6 initial estimates
  S1 <- corHT(z)
  S2 <- corSpearman(z)
  S3 <- corNSR(z)
  S4 <- covMSS(z)
  S5 <- covBACON(z)
  S6 <- rawCovOGK(z)
  Sk <- list(S1,S2,S3,S4,S5,S6)
  
  proper_ev <- function(z,S) { #z cannot be a dataframe but has to be a matrix!
    E <- eigen(S)$vectors
    B <- as.matrix(z)%*%E
    V <- apply(B,2,mad)^2
    L <- diag(V)
    Sigma_hat <- E%*%L%*%t(E)
    mu_hat <- chol(Sigma_hat)%*%apply(as.matrix(z)%*%solve(chol(Sigma_hat)), 2, median)
    result <- list("center"= mu_hat, "scatter" = Sigma_hat)
    return(result)
  }
  
  results_mu <- list()
  results_Sigma <- list()
  for (i in 1:6) {
    result <- proper_ev(z,Sk[[i]])
    results_mu[[i]] <- result[[1]]
    results_Sigma[[i]] <- result[[2]]
  }
#iteration_medians <- list()
#iteration_sigmas <- list()
algorithm <- function(z, mu_hat, Sigma_hat, alpha) {
    z$distances <- mahalanobis(as.matrix(z), mu_hat, Sigma_hat)^(1/2)
    h0 <- ceiling(nrow(z)/2)
    h <- h.alpha.n(alpha,n,p) #don't take distance column into account
    z$rank <- rank(z$distances, ties.method="random")
    initial_subset <- z[z$rank<=h0,1:p] #pick h0 smallest distances as initial subset and delete columns for distance & rank
    #T_H0 <- colSums(initial_subset)/h0
    T_H0 <- colMeans(initial_subset)
    #iteration_medians[[1]] <- T_H0
    #T_H0_rep <- matrix(T_H0,nrow=h,ncol=ncol(initial_subset),byrow=TRUE)
    #S_H0 <- crossprod(as.matrix(initial_subset-T_H0_rep))/h0
    S_H0 <- cov(initial_subset)
    print(T_H0)
    print(S_H0)
    #iteration_sigmas[[1]] <- S_H0
    for (i in 1:1000) {
      # calculate mahalanobis distance for all data points given estimated mean and cov
      z$distances <- mahalanobis(as.matrix(z[,1:p]), T_H0, S_H0)^(1/2)
      # order all mahalanobis distances
      z$rank <- rank(z$distances, ties.method="random")
      # pick h smallest, i.e. set "select' of h smallest ranks to 1
      z$select <- ifelse(z$rank<=h, 1, 0)
      # calculate new mean and cov based on this subset
      #T_H <-colSums(z[z$select==1,1:p])/h #don't use distances, rank & weights to calculate
      #T_H_rep <- matrix(T_H,nrow=h,ncol=p,byrow=TRUE)
      #S_H <- crossprod(as.matrix(z[z$select==1,1:p]-T_H_rep))/h
      T_H <- colMeans(z[z$select==1,1:p])
      print(T_H)
      #iteration_medians[[i+1]] <- T_H
      S_H <- cov(z[z$select==1,1:p])
      print(S_H)
      #iteration_sigmas[[i+1]] <- S_H
      if (det(S_H)==det(S_H0)) { # stopping criterium
        break
      } 
      T_H0 <- T_H
      S_H0 <- S_H
      # repeat
    }
    result <- list("raw.center"=T_H0, "raw.cov"=S_H0,"selection"=z$select)
    #result <- list("raw.center"=iteration_medians, "raw.cov"=iteration_sigmas,"selection"=z$select)
    return(result)
  }
  
  results_raw_center <- list()
  results_raw_cov <- list()
  results_weights <- list()
  
  for (i in 1:6) {
    result <- algorithm(z,results_mu[[i]],results_Sigma[[i]],alpha)
    results_raw_center[[i]] <- result[[1]]
    results_raw_cov[[i]] <- result[[2]]
    results_weights[[i]] <- result[[3]]
  }
  
  best_start <- 1
  best_det <- function(benchmark){
    for (i in 2:6) {
      if (det(results_raw_cov[[i]])<benchmark) {
        best_start <- i
        benchmark <- det(results_raw_cov[[i]])
      }
    }
    result <- list(best_start, benchmark)
    return(result)
  }
  det_1 <- det(results_raw_cov[[1]])
  ind_det <- best_det(det_1) 
  raw.center <- results_raw_center[[ind_det[[1]]]] #pick one with smallest determinant
  raw.cov <- results_raw_cov[[ind_det[[1]]]] #pick one with smallest determinant
  
  z$weights <- results_weights[[ind_det[[1]]]]; #vector with ones and zeros
  z$indices <- seq(1:nrow(z))
  best <- z[z$weights==1,4]
  
  # Reweighting step to transform and obtain final estimate
  Q <- sqrt(qchisq(1-0.025, p, ncp = 0, lower.tail = TRUE, log.p = FALSE))
  z$weights_r <- ifelse(sqrt(mahalanobis(as.matrix(z[,1:p]), raw.center, raw.cov))<=Q, 1, 0)
  weights <- z$weights_r
  center <-colMeans(z[z$weights_r==1,1:p]) 
  cov <- cov(as.matrix(z[z$weights_r==1,1:p]))
  
  #Fisher correction factors
  fisher <- function(frac, p) {
    chisq <- qchisq(frac, p)
    c_alpha <- frac/pgamma(chisq/2, p/2+1,1)
    return(c_alpha)
  }
  fisher_cor_raw <- fisher(h/n, p)
  fisher_cor <- fisher(sum(z$weights_r)/n,p)
  
  Var <- apply(x, 2, mad)
  medians <- apply(x,2,median)
  A <- diag(Var^(-1))
  v <- -medians/Var
  
  #transformed results (in terms of original data)
  center.x <- (center-v)%*%solve(A)
  cov.x <- solve(A)%*%(fisher_cor*cov)%*%solve(A)
  raw.center.x <- (raw.center-v)%*%solve(A)
  raw.cov.x <- solve(A)%*%(fisher_cor_raw*raw.cov)%*%solve(A) 
  
  results <- list("rwgt.center"=center, "rwgt.cov" =cov, "weights" = weights, "raw.center" = raw.center, "raw.cov"=raw.cov, "best.raw"=best, "center.x"=as.numeric(center.x),"cov.x"= cov.x,"raw.center.x"= as.numeric(raw.center.x),"raw.cov.x"= raw.cov.x)#, iteration_medians, iteration_sigmas)
  return(results)
}

E_center <- colMeans(Eredivisie28)
E_cov <- cov(Eredivisie28)
E_radius <- sqrt(qchisq(0.975, df = p))

library(car)
E_ellipse <- data.frame(ellipse(center=E_center,
                                shape=E_cov,
                                radius=E_radius,
                                segments=100,
                                draw=FALSE))

ellipse_mcd_raw <- data.frame(ellipse(center=results_own$raw.center.x,
                                      shape = results_own$raw.cov.x,
                                      radius=sqrt(qchisq(0.975, df = p)),
                                      segments=100,
                                      draw=FALSE))

ellipse_mcd <- data.frame(ellipse(center=results_own$center.x,
                                      shape = results_own$cov.x,
                                      radius=sqrt(qchisq(0.975, df = p)),
                                      segments=100,
                                      draw=FALSE))

colnames(ellipse_mcd_raw) <- colnames(ellipse_mcd) <- colnames(E_ellipse) <- colnames(Eredivisie28)

ggplot(Eredivisie28, aes(Age, MarketValue)) +
  geom_point() +
  geom_polygon(data=ellipse_mcd_raw, color="red", fill="red", alpha=0.3) +
  geom_polygon(data=ellipse_mcd, color="yellow", fill="yellow", alpha=0.3) 
#+ geom_polygon(data=E_ellipse, color="yellow", fill="yellow", alpha=0.3)

# compare to
results <- covMcd(log(Eredivisie28), alpha=0.75, nsamp="deterministic", use.correction=FALSE)

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
