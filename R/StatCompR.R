#' @title Dataset x for test of regression models
#' @name samplex
#' @description A dataset simulated from sparse singular values and is used to carry out \code{FES} and \code{PLS0}.
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
NULL


#' @title Dataset y for test of regression models
#' @name sampley
#' @description A dataset simulated from multivariate normal distribution and is used to carry out \code{FES} and \code{PLS0}.
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
NULL


#' @title FES method for multiple regression. 
#' @description Factor estimation and selection with the tuning parameter selected algorithm with R. 
#' @param x0 a n*p matrix of explanatory variables
#' @param y0 a n*q matrix of responses
#' @return Return the estimation of coefficient B. 
#' @examples
#' \dontrun{
#' data(samplex)
#' attach(samplex)
#' data(sampley)
#' attach(sampley)
#' B.hat <- FES(samplex[[1]], sampley[[1]])
#' }
#' @export
FES <- function(x0, y0){
  
  B.hat <- function(t, x, y){
    B.OLS <- solve(t(x) %*% x) %*% t(x) %*% y 
    a <- svd(B.OLS)
    D.LS <- a$d
    
    sum.Dii.hat <- min(t, sum(a$d))
    D.LS <- c(D.LS, 0)
    for(k in 1:8){
      k <- 9 - k
      lambda <- (sum(D.LS[1:k]) - sum.Dii.hat)/k
      if(lambda <= D.LS[k] && lambda >= D.LS[k+1]) break
    }
    
    D.hat <- matrix(rep(0, 8*8), nrow = 8, ncol = 8)
    for(i in 1:8) D.hat[i, i] <- max(D.LS[i] - lambda, 0)
    B <- a$u %*% D.hat %*% t(a$v)
    
    K <- matrix(rep(0, 8*8), nrow = 8, ncol = 8)
    for(i in 1:8){
      if(D.hat[i, i] > 0) K <- K + a$v[, i] %*% t(a$v[, i]) / D.hat[i, i]
    }
    
    df <- 8 * sum(diag(x %*% solve(t(x) %*% x + 2 * 20 * lambda * K) %*% t(x)))
    GCV <- sum(diag((y - x %*% B) %*% t(y - x %*% B)))/(8*8 - df)
    
    return(list(B, GCV))
  }
  
  GCVscore <- rep(0, 100)
  t <- seq(from = 5, to = 25, length.out = 100)
  
  for(i in 1:100) GCVscore[i] <- B.hat(t[i], x0, y0)[[2]]
  t0 <- t[which.min(GCVscore)]
  B0 <- B.hat(t0, x0, y0)[[1]]
  
  return(B0)
}


#' @title PLS method for multiple regression. 
#' @description Two-block partial least squares algorithm with R. 
#' @param x0 a n*p matrix of explanatory variables
#' @param y0 a n*q matrix of responses
#' @return Return the estimation of coefficient B. 
#' @importFrom stats coef sd
#' @importFrom pls plsr
#' @examples
#' \dontrun{
#' data(samplex)
#' attach(samplex)
#' data(sampley)
#' attach(sampley)
#' B.hat <- PLS0(samplex[[1]], sampley[[1]])
#' }
#' @export
PLS0 <- function(x0, y0){
  estimate <- plsr(y0 ~ x0, scale=T, validation="CV", segments = 10)
  k <- estimate$validation$ncomp
  
  estimatek <- plsr(y0 ~ x0, ncomp = k, scale=T, validation="CV")
  estimatek.coef <- coef(estimatek)/apply(y0, 2, sd) # regression coefficient
  estimatek.int <- apply(y0, 2, mean) - estimatek.coef[,,1] %*% apply(x0, 2, mean) # intercept
  return(estimatek.coef[,,1])
}


#' @title Model error for N groups. 
#' @description Mean and standard error of model error for N groups. 
#' @param N number of groups of data "x" and "y"
#' @param x a n*p matrix of explanatory variables
#' @param y a n*q matrix of explanatory variables
#' @param B true coefficient matrix for generating samples
#' @param sigma true correlation matrix for generating dataset "x"
#' @param method choice of regression method, 1 for FES and 2 for PLS
#' @return Return the estimation of coefficient B. 
#' @examples
#' \dontrun{
#' data(samplex)
#' attach(samplex)
#' data(sampley)
#' attach(sampley)
#' result <- modelerror(200, samplex, sampley, samplex[[201]], samplex[[202]], method = 1)
#' }
#' @export
modelerror <- function(N = 200, x, y, B, sigma, method){
me <- rep(0, N)
if(method == 1){
  for(m in 1:N){
    B0 <- FES(x[[m]], y[[m]])
    me[m] <- mean(abs(t(B0 - B) %*% sigma %*% (B0 - B)))
  }
}else if(method == 2){
  for(m in 1:N){
    B0 <- PLS0(x[[m]], y[[m]])
    me[m] <- mean(abs(t(B0 - B) %*% sigma %*% (B0 - B)))
  }
}
return(c(mean(me), sd(me)))
}





 