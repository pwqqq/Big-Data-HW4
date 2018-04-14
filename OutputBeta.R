############### USE ADMM ALGORITHM ESTIMATE BETA ###############

# Write a function to output beta as csv file
OutputBeta <- function(X, y, alpha, lambda){
  # Create null matrix
  beta <- matrix(,nrow = ncol(X), ncol = length(lambda))
  
  for( i in 1:length(lambda)){
    
    tau <- length(y)*lambda[i]
    
    beta[,i] <- admmEN(X, y, alpha, tau=tau)
    
    names <- paste("Lambda","=",lambda)
  }
  # name columns
  colnames(beta) <- names
  
  return(beta)
  
}

# Create lambda sequence
Lseq <- c(0.01,0.1,1,10)

# alpha = 0.95
beta1 <- OutputBeta(X,y,alpha = 0.95, Lseq)

write.csv(beta1, "b1.csv", row.names = FALSE)


# alpha = 1 - Lasso
beta2 <- OutputBeta(X,y,alpha = 1, Lseq)

write.csv(beta2, "b2.csv", row.names = FALSE)