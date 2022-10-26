rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

#### Adapted QIS function taken from 
#### https://github.com/MikeWolf007/covShrinkage/blob/main/qis.R
#### Only the output and the name of the function was changed
QIS <- function(Y, k = -1) {
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


VIteration <- function(Z, Lambda){
  
  tol<-1e-4
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  V0 <- V <-diag(p)
  Lambdainv<-Lambda^(-1)
  
  
  i=0
  diagnorm=1000
  
  
  while (abs(diagnorm) > tol){
    
    i<-i+1
    
    if (i > 200){
    warning("More than 200 iterations")
    break
    }
    
    V_old=V
    Hinv_old=V_old%*%diag(Lambdainv)%*%t(V_old)
    H_old=V_old%*%diag(Lambda)%*%t(V_old)
    
    lower <- diag( (Z%*%Hinv_old%*%t(Z))/p )
    ratioV <- t(Z/lower)%*%Z
    # Make sure ratioV
    ratioV <- (ratioV + t(ratioV))/2
    
    eigen_Psi<-eigen(ratioV)
    V=eigen_Psi[[2]][,order(eigen_Psi[[1]])]
    
    
    diagnorm=norm(t(V)%*%ratioV%*%V%*%diag(Lambdainv)
                  - diag(Lambdainv)%*%t(V)%*%ratioV%*%V, "F")
    
  
  }
 


  return(V)
  
}


RNL <- function(Y, C=F, cov=T){
  
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  Y<-scale(Y,scale=F)
  traceShat<-sum( diag(var(Y)) )
  
  if (C==T){
    
    dsquare=diag(var(Y))
    Y=Y/sqrt(dsquare)
    
  }else{
    
    dsquare=rep(1,p) 
  }
  
  Z<-t(apply(Y,1, function(x) x/sqrt(sum(x^2)) ))
       
  Lambda<-QIS(Z)[[2]]
  Lambda<-sort(Lambda)
  
  V<-VIteration(Z=Z, Lambda=Lambda)
  
  # 
  H_inv <- V%*%diag(Lambda^(-1))%*%t(V)
  lower <- diag( (Z%*%H_inv%*%t(Z))/p )
  Y<-Z/sqrt(lower)
  
  H0<-QIS(Y)[[1]]
  
  H <- diag( sqrt(dsquare))%*%H0%*%diag( sqrt(dsquare))
  
  if (cov==T){
   # If we are confident the covariance matrix exists
  # (i.e., the data is not too heavy-tailed)
    return(traceShat*H/sum(diag(H)))
   }else {
     # If not, sum( diag(var(Y)) ) does not make sense and so
     # we standardize by p.
     return(p*H/sum(diag(H))) 
    }
  
}








