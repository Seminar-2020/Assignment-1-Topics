#### Functions for initial estimators ####
# Input: the standardized data matrix z
# Output: the estimated covariance or correlation matrix
# Please do not use any other input or output for the initial estimators

# correlation matrix based on the hyperbolic tangent s.1 
corHT <- function(z) {
  y  = atan(z)
  s.1 = cor(y, method = "pearson")
  return(s.1)
}
s.1 <-corHT(df.s)
  
# spearman correlation matrix s.2
corSpearman <- function(z) {
  r <- apply(df.s, 2, rank)
  s.2 = cor(r, method = "pearson")
  return(s.2)
}
s.2 <-corSpearman(df.s)

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
s.3 <- corNSR(df.s)

#euclidean distance function 
euc.dist<- function(x1) sqrt(sum((x1) ^ 2))

# modified spatial sign covariance matrix s.4
covMSS <- function(z) {
     d  <-  apply(z,1, euc.dist)
     k  <- df.s/d
     k[is.na(k)] <- 0  #comment, paper sign, if function
     s.4 <- 1/nrow(z)*crossprod(k)
     return(s.4)
} 
s.4 <- covMSS(df.s)
# covariance matrix based on first step of BACON #Juist variables, beter schrijven s.5
covBACON <- function(z) {
  #Based on the norm 
  dist <- apply(z, 1, euc.dist)
  zz <- cbind(z,dist)
  zz <- zz[order(dist),]
  z.s <- zz[1:ceiling(nrow(df.s)/2), 1:ncol(z)]
  s.5 <- cov(z.s)
  return (s.6)
} 
s.5 <- covBACON(df.s)

# raw OGK estimator of the covariance matrix with median and Qn s.6
rawCovOGK <- function(z) {
  raw<-covOGK(z, sigmamu = s_Qn, rcov = covGK, weight.fn = hard.rejection)
  s.6 <- raw$cov
  return(s.6)
}  #documentation 
s.6 <- rawCovOGK(df.s)


