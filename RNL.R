
### Some helper function
rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

#### Adapted QIS function taken from 
#### https://github.com/MikeWolf007/covShrinkage/blob/main/qis.R
#### Only the output and the name of the function was changed
myQIS <- function(Y, k = 1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  sample <- (t(Y) %*% Y) / n    # sample covariance matrix    
  spectral <- eigen(sample)    # spectral decompositon
  lambda <- spectral$values[p:1]    # sort eigenvalues in ascending order
  u <- spectral$vectors[,p:1]    # eigenvectors follow their eigenvalues
  h <- min(c^2, 1/c^2)^0.35 / p^0.35    # smoothing parameter
  invlambda <- 1 / lambda[max(1, p-n+1):p]    # inverse of non-null eigenvalues   
  Lj <- rep.row(invlambda, min(p, n))    # like 1 / lambda_j
  Lj.i <- Lj - t(Lj)    # like (1 / lambda_j) - (1 / lambda_i)
  theta <- rowMeans(Lj * Lj.i / (Lj.i^2 + h^2 * Lj^2))    # smoothed Stein shrinker
  Htheta <- rowMeans(Lj * (h * Lj) / (Lj.i^2 + h^2 * Lj^2)) # its conjugate
  Atheta2 <- theta^2 + Htheta^2    # its squared amplitude
  if (p <= n)    # case where sample covariance matrix is not singular
    delta <- 1 / ((1 - c)^2 * invlambda + 2 * c * (1 - c) * invlambda * theta +
                    c^2 * invlambda * Atheta2)           # optimally shrunk eigenvalues
  else {    # case where sample covariance matrix is singular
    delta0 <- 1 / ((c - 1) * mean(invlambda))     # shrinkage of null eigenvalues
    delta <- c(rep(delta0, p - n), 1 / (invlambda * Atheta2));
  }
  deltaQIS <- delta * (sum(lambda) / sum(delta))    # preserve trace
  sigmahat <- u %*% diag(deltaQIS) %*% t(u)    #reconstruct covariance matrix
  
  return(list(Sig=sigmahat, Lambda=deltaQIS, U=u))
  
}

# Helper function used in RNL, where eig_fix is a vector of eigenvalues.
# Note: The first p-n elements in eig_fix necessarily correspond to the zero eigenvalues.
sortLambda <- function(eig_fix, n, p) {
  eig_fix <- sort(eig_fix)
  if (p > n) {
    # First eigenvalues = smallest eigenvalues (in this case corresponding
    # to the sample eigenvalues with zeros)
    eig_fix_0 <- eig_fix[1:(p - n + 1)]
    un <- unique(eig_fix_0)
    if (length(un) > 1) {
      nr <- rep(0, length(un))
      for (i in 1:length(un)) {
        nr[i] <- sum(eig_fix_0 == un[i])
      }
      ind <- which.max(nr)
      eig_fix_0 <- rep(un[ind], length(eig_fix_0))
    }
    
    eig_fix_1 <- sort(eig_fix[(p - n + 2):length(eig_fix)])
    eig_fix <- c(eig_fix_0, eig_fix_1)
  }
  return(eig_fix)
}

### helper function for RNL, updating the eigenvectors.
VIteration <- function(Z, Lambda){
  
  tol <- 1e-5
  maxit <- 10
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  V0 <- V <-diag(p)
  Lambdainv <- Lambda^(-1)
  
  i <- 0
  diagnorm <- 1000
  
  lam = rep(1,p)
  H_init = V0%*%(diag(Lambda))%*%t(V0)
  
  lower <- diag((Z%*%solve(H_init)%*%t(Z))/p)
  ratioV <- (t(Z/lower)%*%Z)/n
  ratioV <- (ratioV + t(ratioV))/2
  
  while ((abs(diagnorm) > tol) & (i < maxit)){
    
    i <- i+1
    
    V_old <- V
    ratioV_old <- ratioV
    
    Hinv_old <- V_old%*%diag(Lambdainv)%*%t(V_old)
    H_old <- V_old%*%diag(Lambda)%*%t(V_old)
    
    lower <- diag((Z%*%Hinv_old%*%t(Z))/p)
    ratioV <- (t(Z/lower)%*%Z)/n
    ratioV <- (ratioV + t(ratioV))/2
    
    eigen_Psi<-eigen(ratioV)
    V=eigen_Psi[[2]][,order(eigen_Psi[[1]])]
    
    diagnorm <- norm(t(V_old)%*%ratioV_old%*%V_old%*%diag(Lambdainv)
                     - diag(Lambdainv)%*%t(V)%*%ratioV%*%V, "F")
    
  }
  
  return(V)
  
}

### Robust nonlinear shrinkage, where Y is an nxp data matrix and C = T refers to the R-C-NL approach.
RNL <- function(X, C = T){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  X <- scale(X,scale=F)
  traceShat <- sum(diag(var(X)))
  
  if (C == T){
    
    dsquare <- diag(var(X))
    X <- t(apply(X,1,function(y){y/sqrt(dsquare)}))
    
  } else {
    
    dsquare <- rep(1,p) 
  }
  
  Z <- t(apply(X,1, function(x){x/sqrt(sum(x^2))}))
       
  Lambda<-myQIS(Z)[[2]]
  Lambda<-sortLambda(Lambda,n,p)
  
  if(any(Lambda<1e-10)){
    H0 <- myQIS(X)[[1]]
  } else {
    V <- VIteration(Z=Z,Lambda=Lambda)
    
    H_inv <- V%*%diag(Lambda^(-1))%*%t(V)
    lower <- diag((Z%*%H_inv%*%t(Z))/p)
    X <- Z/sqrt(lower)
    
    H0 <- myQIS(X)[[1]]
  }
  
  H <- diag(sqrt(dsquare))%*%H0%*%diag(sqrt(dsquare))
  
  return(traceShat*H/sum(diag(H)))
  
}








