soft_thresholding <- function(x,a){
  ## This could be done more efficiently using vector multiplication
  ## See the forumula in slides
  ##  sign(x)*pmax(abs(x) - a, 0)
  result <- numeric(length(x))
  result[which(x > a)] <- x[which(x > a)] - a
  result[which(x < -a)] <- x[which(x < -a)] + a
  return(result)
}



lasso_kkt_check <- function(X,y,beta,lambda, tol=1e-3){
  ## check convergence 
  beta <- as.matrix(beta); X <- as.matrix(X)
  ## Assuming no intercepts 
  G <- t(X)%*%(y-X%*%beta)/length(y)
  ix <- which(beta == 0 )
  iy <- which(beta != 0)
  if (any(abs(G[ix]) > (lambda + tol) )) { return(pass=0) }
  if (any(abs( G[iy] - lambda*sign(beta[iy] ) ) > tol)) { return(pass=0) }  
  return(pass=1)
}



##### lasso.en ####

en <- function(X,y,beta,lambda,alpha, tol=1e-6,maxiter=1000,quiet=FALSE){
  # note that the LS part  in this function is the one in slides divided by length(y) = n 
  ## Or equivalently  lambda here = n * lambda in class
  beta <- as.matrix(beta); X <- as.matrix(X)
  obj <- numeric(length=(maxiter+1))
  betalist <- list(length=(maxiter+1))
  betalist[[1]] <- beta
  
  for (j in 1:maxiter){
    for (k in 1:length(beta)){
      r <- y - X[,-k]%*%beta[-k]
      beta[k] <- (1/norm(as.matrix(X[,k]),"F")^2)*soft_thresholding(t(r)%*%X[,k],length(y)*lambda*alpha) / (1+lambda*(1-alpha))
    }
    betalist[[(j+1)]] <- beta
    obj[j] <- (1/2)*(1/length(y))*norm(y - X%*%beta,"F")^2 + lambda*(alpha*sum(abs(beta)) + (1/2) * (1 - alpha) * sum(beta^2))
    if (norm(betalist[[j]] - beta,"F") < tol) { break }
  } 
  check <- lasso_kkt_check(X,y,beta,lambda) 
  
  if (quiet==FALSE){
    if (check==1) {
      cat(noquote("Minimum obtained.\n"))
    }
    else { cat(noquote("Minimum not obtained.\n")) } 
  }
  return(list(obj=obj[1:j],beta=beta)) 
}

#### model ####
X <- read.csv("X.csv", header = FALSE)
Y <- read.csv("y.csv", header = FALSE)
colnames(Y) <- ("y")
dataset <- cbind(X,Y)

x<-model.matrix(y~ ., data = dataset)[,-1]
y <- dataset$y


#### model test ####
## could try, will delete in the final version
#re.en <- en(x,y,rep(0,400),lambda = 1,alpha = 0.95)
# re.en.beta <- re.en$beta[which(re.en$beta != "0")]
#fit <- glmnet(x, y, lambda = 1, alpha = 0.95, standardize  = F, intercept = F)
#re.glmnet <- as.numeric( coef(fit) )[-1]
# re.glmnet.beta <- re.glmnet[which(re.glmnet != "0")]
# lambda.test <- c(0.01,0.1,1,10)

################################# Q1 #################################
lambda.test <- c(0.01,0.1,1,10)
b1 <- matrix(NA,nrow=400,ncol=4)
system.time( for (i in 1:4){
  re.en <- en(x,y,rep(0,400),lambda = lambda.test[i],alpha = 0.95)
  b1[,i] <- re.en$beta
})

################################# Q2 #################################
b2 <- matrix(NA,nrow=400,ncol=4)
system.time( for (i in 1:4){
  re.en <- en(x,y,rep(0,400),lambda = lambda.test[i],alpha = 1)
  b2[,i] <- re.en$beta
})

#write.csv(b1,"b1.csv")
#write.csv(b2,"b2.csv")