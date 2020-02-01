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

# Output
# A list with the following components:
# center ....... mean of the reweighted estimator
# cov .......... covariance matrix of the reweighted estimator
# weights ...... binary weights of the observations used for the reweighted
#                estimator (0 if outlier, 1 if used for estimation)
# weights_final. binary vector indicating if the observations falls within 
#                the confidence ellips of the reweighted estimators 
# raw.center ... mean of the raw estimator
# raw.cov ...... covariance matrix of the raw estimator
# best ......... indices of the observations in the best subset found by the
#                raw estimator
# raw.n.outliers number of outliers according to the raw estimator
# n.outliers ... number of outliers according to the reweighted estimator

covDetMCD <- function(x, alpha) {
  n <- nrow(x)
  p <- ncol(x)
  h <- h.alpha.n(alpha,n,p)
  z <- as.data.frame(rob_standard(x))
  
  #obtain 6 initial estimates and place them in a list
  Sk <- list("S1" = corHT(z),
             "S2" = corSpearman(z),
             "S3" = corNSR(z),
             "S4" = covMSS(z),
             "S5" = covBACON(z),
             "S6" = rawCovOGK(z))
  
  results_mu <- list()
  results_Sigma <- list()
  for (i in 1:6) {
    result <- proper_ev(z,Sk[[i]])
    results_mu[[i]] <- result[[1]] #list with 6 initial center estimates
    results_Sigma[[i]] <- result[[2]] #list with 6 initial scatter estimates 
  }
  
  results_raw_center <- list()
  results_raw_cov <- list()
  results_selection <- list()
  results_distances <- list()
  
  for (i in 1:6) {
    result <- algorithm(z,results_mu[[i]],results_Sigma[[i]],alpha)
    results_raw_center[[i]] <- result[[1]]
    results_raw_cov[[i]] <- result[[2]]
    results_selection[[i]] <- result[[3]]
  }
  
  ind_det <- best_det(det(results_raw_cov[[1]])) #use the determinant of estimate 1 as an initial benchmark 
  raw.center.z <- results_raw_center[[ind_det[[1]]]] 
  fisher_cor_raw <- fisher(h/n, p)
  raw.cov.z <- fisher_cor_raw*results_raw_cov[[ind_det[[1]]]]
  z$selected <- results_selection[[ind_det[[1]]]]
  z$indices <- seq(1:nrow(z))
  best <- z[z$selected==1,4]
  
  #reweighting step
  Q <- sqrt(qchisq(0.975, p))
  obs_out_raw <- sum(sqrt(mahalanobis(as.matrix(z[,1:p]), raw.center.z, raw.cov.z))>Q)
  z$weights_r <- ifelse(sqrt(mahalanobis(as.matrix(z[,1:p]), raw.center.z, raw.cov.z))<=Q, 1, 0)
  weights <- z$weights_r
  center.z <-colMeans(z[z$weights_r==1,1:p]) 
  fisher_cor <- fisher(sum(z$weights_r)/n,p)
  cov.z <- fisher_cor*cov(as.matrix(z[z$weights_r==1,1:p]))
  obs_out_reweighted <- sum(sqrt(mahalanobis(as.matrix(z[,1:p]), center.z, cov.z))>Q)
  weights_final <- z$weights_final <- ifelse(sqrt(mahalanobis(as.matrix(z[,1:p]), center.z, cov.z))<=Q,1,0)
  
  #reverse standardization, transform back to the units
  Var <- apply(x, 2, mad)
  medians <- apply(x,2,median)
  A <- diag(Var^(-1))
  v <- -medians/Var
  center <- (center.z-v)%*%solve(A)
  cov <- solve(A)%*%cov.z%*%solve(A)
  raw.center <- (raw.center.z-v)%*%solve(A)
  raw.cov <- solve(A)%*%raw.cov.z%*%solve(A) 
  
  #gather all important results
  results <- list("center"=as.numeric(center), 
                  "cov" =cov, 
                  "weights" = weights,
                  "weights_final"=weights_final,
                  "raw.center" = as.numeric(raw.center), 
                  "raw.cov"=raw.cov, 
                  "best"=best, 
                  "raw.n.outliers"=obs_out_raw,
                  "n.outliers"=obs_out_reweighted)
  return(results)
  
  # standardizes a dataframe by subtracting the median and dividing by the MAD (instead of Qn proposed in paper)
  # returns the standardized dataframe
  rob_standard <- function(x) {
    Var <- apply(x, 2, mad)
    medians <- apply(x,2,median)
    A <- diag(Var^(-1))
    ones <- matrix(1, nrow=n, ncol=1)
    v <- medians/Var
    z <- as.matrix(x)%*%A-ones%*%v
    return(z)
  }
  # transform S into initial scatter estimate with accurate eigenvalues
  # obtain an initial estimate for the center
  # returns a list with covariance and center estimate for z
  proper_ev <- function(z,S) {
    E <- eigen(S)$vectors
    B <- as.matrix(z)%*%E
    V <- apply(B,2,mad)^2
    L <- diag(V)
    Sigma_hat <- E%*%L%*%t(E)
    mu_hat <- chol(Sigma_hat)%*%apply(as.matrix(z)%*%solve(chol(Sigma_hat)), 2, median)
    result <- list("center"= mu_hat, "scatter" = Sigma_hat)
    return(result)
  }
  #find the optimal subset of h observations given the initial center and scatter estimates
  #return the estimated center, scatter and which h variables are selected in the subset 
  algorithm <- function(z, mu_hat, Sigma_hat, alpha) {
    #compute statistical distances and select h0 smallest
    z$distances <- mahalanobis(as.matrix(z), mu_hat, Sigma_hat)^(1/2)
    h0 <- ceiling(nrow(z)/2)
    h <- h.alpha.n(alpha,n,p) 
    z$rank <- rank(z$distances, ties.method="random") #method set to "random" to avoid similar ranks
    initial_subset <- z[z$rank<=h0,1:p] 
    #compute new center and scatter estimate based on the concentrated subset
    T_H0 <- colMeans(initial_subset)
    S_H0 <- cov(initial_subset)
    #compute distances based on this subset and iterate this procedure until convergence
    for (i in 1:1000) {
      z$distances <- mahalanobis(as.matrix(z[,1:p]), T_H0, S_H0)^(1/2)
      z$rank <- rank(z$distances, ties.method="random")
      z$select <- ifelse(z$rank<=h, 1, 0) #pick h smallest, i.e. set "select' of h smallest ranks to 1
      T_H <- colMeans(z[z$select==1,1:p])
      S_H <- cov(z[z$select==1,1:p])
      if (det(S_H)==det(S_H0)) { # stopping criterion
        break
      } 
      #if not yet converged: update estimates and iterate again
      T_H0 <- T_H
      S_H0 <- S_H
    }
    result <- list("raw.center"=T_H0, "raw.cov"=S_H0,"selection"=z$select)
    return(result)
  }
  #identify which of the 6 estimators has the lowest determinant
  #return a list with the index of the chosen initial estimator and the corresponding determinant value
  best_det <- function(benchmark){
    best_start <- 1
    for (i in 2:6) {
      if (det(results_raw_cov[[i]])<benchmark) {
        best_start <- i
        benchmark <- det(results_raw_cov[[i]])
      }
    }
    result <- list(best_start, benchmark)
    return(result)
  }
  #compute the fisher correction factors based on the fraction of observations used and the number of variables
  fisher <- function(frac, p) {
    chisq <- qchisq(frac, p)
    c_alpha <- frac/pgamma(chisq/2, p/2+1,1)
    return(c_alpha)
  }
}

#this function uses the output of covDetMCD to construct a plot of the data and confidence ellipses
plot_ellipses <- function(data,mcd_obj) {
  require(car)
  require(ggplot2)
  ellipse_mcd_raw <- data.frame(ellipse(center=mcd_obj$raw.center,
                                        shape = mcd_obj$raw.cov,
                                        radius=sqrt(qchisq(0.975, df = p)),
                                        segments=100,
                                        draw=FALSE))
  
  ellipse_mcd <- data.frame(ellipse(center=mcd_obj$center,
                                    shape = mcd_obj$cov,
                                    radius=sqrt(qchisq(0.975, df = p)),
                                    segments=100,
                                    draw=FALSE))
  
  colnames(ellipse_mcd_raw) <- colnames(ellipse_mcd) <- colnames(data)
  
  ggplot(data, aes(Age, MarketValue)) +
    geom_point() +
    geom_polygon(color="red", fill="red", alpha=0.3, data=ellipse_mcd_raw) +
    geom_polygon(color="yellow", fill="yellow", alpha=0.3, data=ellipse_mcd)
}
plot_ellipses(as.data.frame(log(Eredivisie28)), covDetMCD(as.data.frame(log(Eredivisie28)),alpha=0.75))

covDetMCD(as.data.frame(log(Eredivisie28)),alpha=0.75)

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
