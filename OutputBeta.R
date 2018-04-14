authors <- function() {
  c("Annie Yang", "Fuyu Zou", "Lisha Wang", "Weiqi Pan")
}

## Install packages 
install.packages("glmnet")
library(glmnet)

## Import dataset 
X <- read.csv("X.csv",header = FALSE)
X <- as.matrix(X)

y <- read.csv("y.csv",header = FALSE)
y <- as.matrix(y)

n <- nrow(X)

## soft threshold ## 
softThresh <- function(x, lambda) {
  sign(x)*pmax(0, abs(x) - lambda)
}

###### Use ADMM algorithm and write function admmEN to estimate beta ######

admmEN <- function(X, y, tau, alpha=0.95, maxit = 10000, tol=1e-7) {# tau = n*lambda
  XX <- t(X) %*% X
  Xy <- t(X) %*% y
  
  p <- ncol(X) # number of variables
  L <- rep(0, p) # initialize Lagrange multiplier
  rho <- 4 
  
  z0 <- z <- beta0 <- beta <- rep(0, p)
  
  ## Compared to ADMMLasso, add tau*(1-alpha)*diag(rep(1, p)) because of the L2 regularization
  Sinv <- solve(XX + tau*(1-alpha)*diag(rep(1, p)) + rho*diag(rep(1, p)))
  
  for (it in 1:maxit) {
    ## update beta
    ## beta <- solve(XX + tau*(1-alpha)*diag(rep(1, p)) + rho*diag(rep(1, p)) ) %*% (Xy + rho * z - lambda)
    beta <- Sinv %*% (Xy + rho * z - L)
    
    ## update z
    z <- softThresh(beta + L/rho, tau*alpha/rho) # Compared to tau/rho lasso, here we multiply it by alpha  
    
    ## update Lagrange multiplier
    L <- L + rho* (beta - z) 
    
    change <- max(  c( base::norm(beta - beta0, "F"),
                       base::norm(z - z0, "F") ) )
    if (change < tol || it > maxit) { # Check convergence
      break
    }
    beta0 <-  beta
    z0 <-  z
    
  }
  z
}

###### Write a function to output beta as csv file ######

OutputBeta <- function(X, y, alpha, lambda){
  # Create null matrix
  beta <- matrix(NA,nrow = ncol(X), ncol = length(lambda))
  for( i in 1:length(lambda)){
    tau <- length(y)*lambda[i]
    beta[,i] <- admmEN(X, y, alpha, tau=tau)
    names <- paste("Lambda","=",lambda)
  }
  # Name columns
  colnames(beta) <- names
  return(beta)
}

# Create lambda sequence
Lseq <- c(0.01,0.1,1,10)


############################## Q1 ##############################
## Set alpha = 0.95
beta1 <- OutputBeta(X,y,alpha = 0.95, Lseq)
write.csv(beta1, "b1.csv", row.names = FALSE)


############################## Q2 ##############################
## Set alpha = 1
beta2 <- OutputBeta(X,y,alpha = 1, Lseq)
write.csv(beta1, "b2.csv", row.names = FALSE)


############################## DISCUSSION ##############################
## we also tried another three methods to test the cofficients
## all the methods return same results with same alpha and lambda
## chose admm becasue admm is much faster than cd algorithm  

## ADMM2 - estimate beta by construct an artificial data set (ynew, Xnew)
admmEN2 <- function(X, y, tau, alpha, maxit = 10000, tol=1e-7) {
  
  lambda2 <- tau*(1-alpha)
  
  n <- ncol(X)
  
  I <- sqrt(lambda2)*diag(rep(1,n))
  
  Xnew <- as.matrix(rbind(X,I)/sqrt(1+lambda2))
  
  ynew <- rbind(y,as.matrix(rep(0,n)))
  
  XX <- t(Xnew) %*% Xnew
  Xy <- t(Xnew) %*% ynew
  
  p <- ncol(Xnew)
  L <- rep(0, p)
  rho <- 4
  
  z0 <- z <- beta0 <- beta <- rep(0, p)
  Sinv <- solve(XX + rho*diag(rep(1, p)))
  
  for (it in 1:maxit) {
    ## update beta
    ## beta <- solve(XX + rho*diag(rep(1, p)) ) %*% (Xy + rho * z - lambda)
    beta <- Sinv %*% (Xy + rho * z - L)
    
    ## update z
    z <- softThresh(beta + L/rho, tau*alpha/(sqrt(1+lambda2)*rho))
    
    ## update Lagrange multiplier
    L <- L + rho* (beta - z) 
    
    change <- max(  c( base::norm(beta - beta0, "F"),
                       base::norm(z - z0, "F") ) )
    if (change < tol || it > maxit) {
      break
    }
    beta0 <-  beta
    z0 <-  z
    
  }
  return (z/sqrt(1+lambda2))
}

## CD1
EN.cd1 <- function(X,y,beta,alpha, tau ,tol=1e-7,maxiter=100000,quiet=FALSE){# tau = n*lambda
  
  beta <- as.matrix(beta); X <- as.matrix(X)
  betalist <- list(length=(maxiter+1))
  betalist[[1]] <- beta
  
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      
      r <- y - X[,-k]%*%beta[-k]
      
      beta[k] <- (1/(base::norm(as.matrix(X[,k]),"F")^2 + (1-alpha)*tau))*softThresh(t(r)%*%X[,k],alpha*tau)
    }
    
    betalist[[(j+1)]] <- beta
    
    if (norm(betalist[[j]] - beta,"F") < tol) { break }
  } 
  
  return (beta) 
}

## CD2 - estimate beta by construct an artificial data set (ynew, Xnew)
EN.cd2 <- function(X,y,beta,alpha, tau ,tol=1e-7,maxiter=100000,quiet=FALSE){
  
  lambda2 <- tau*(1-alpha)
  
  I <- sqrt(lambda2)*diag(rep(1,length(beta)))
  
  Xnew <- as.matrix(rbind(X,I)/sqrt(1+lambda2))
  
  ynew <- rbind(y,as.matrix(rep(0,length(beta))))
  
  beta <- as.matrix(beta)
  
  obj <- numeric(length=(maxiter+1))
  betalist <- list(length=(maxiter+1))
  betalist[[1]] <- beta
  
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      
      r <- ynew - Xnew[,-k]%*%beta[-k]
      
      beta[k] <- (1/norm(as.matrix(Xnew[,k]),"F")^2)*softThresh(t(r)%*%Xnew[,k],tau*alpha/sqrt(1+lambda2))
    }
    
    betalist[[(j+1)]] <- beta
    
    if (norm(betalist[[j]] - beta,"F") < tol) { break }
  } 
  
  return (beta/sqrt(1+lambda2)) 
}


## We could then use two methods to test the difference between the estimated coefficients
## First, we could use norm to test the difference between two coefficient vectors, and comparing the result with glmnet package
## For example, set alpha = 0.95 and lambda = 10
re.admm <- admmEN(X, y, alpha = 0.95, tau = 10*n)
## use glmnet package to test the result
## change default argument in glmnet to solve problems here: Set standardize and intercept equal FALSE 
## so we could plug in the original dataset directly
re.EN <- as.matrix(coef(glmnet(X, y, lambda = 10, alpha = 0.95, standardize  = F, intercept = F))[-1])
## use norm to test the difference between beta vectors
base::norm(re.admm-re.EN,"F")

## Second, we could the estimated objective to test the estimation accuracy following the function that Dr.Luo shows on the hw4
## Multiply the objective function by n. And set tau = n*lambda. We write the following obj function:

obj <- function(beta, X, y, tau,alpha) {
  1/2*base::norm(y - X%*%beta, "F")^2 + 0.5*tau*(1-alpha)*base::norm(beta,"F")^2 + tau*alpha*sum(abs(beta))
}

## In this case, when alpha is 1, which is the lasso regression, the results are the same with glmnet package
## You can run the following codes to test the lasso situation:

# re.admm_l <- admmEN(X, y, alpha = 1, tau = 1*n)
# re.EN_l <- as.matrix(coef(glmnet(X, y, lambda = 1,alpha = 1, standardize  = F, intercept = F))[-1])

# obj_admm <- obj(re.admm_l, X, y, tau = 1*n, alpha = 1)
# obj_glmnet <- obj(re.EN_l, X, y, tau = 1*n, alpha = 1)

## You will get the same result as 1300.847.


## Also, when lambda is 0.01 and alpha is 0.95, the results are also very close.
## Although not perfectly match.
## Conclusion:the results returned from obj function using admm estimated beta 
## are all smaller than or equal to the results using beta from glmnet function.
## And the difference is small.


