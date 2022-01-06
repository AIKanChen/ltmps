#' @title MPS estimation of the left-truncated Burr distribution
#' @description Apply MPS estimation to the left-truncated Burr distribution.
#' @param data Empirical data that is left-truncated.
#' @param init_value  The initial value of left-truncated Burr distributions parameter alpha, gamma and theta.
#' @param u The left truncated point.
#' @import stats
#' @importFrom actuar rburr
#' @return A \code{list} contains the estimation result of parameters alpha, gamma and theta, standard error, and covariance matrix of estimated parameters.
#' @export
#'
#' @examples
#' data0 = rburr(3000,1,4,0.2)
#' data = data0[data0>5]
#' fit = burr_LTMPS(data = data,init_value = c(1,4,0.2),u=5)
#' fit$par.ests

burr_LTMPS<-function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MPS estimation function
  opt_burr <- function(par){
    k=exp(par[1]);c=exp(par[2]);lambda=exp(par[3])
    est=-sum(log(c(1, 1-(1+(dat/lambda)^c)^(-k))-
                   (1-(1+(c(dat,u)/lambda)^c)^(-k))))+
      (n+1)*log((1+(u/lambda)^c)^(-k))
    return(est)
  }
  fit <- optim(log(init_value),opt_burr,hessian = T)
  #result
  par.ests <- exp(fit$par)
  varcov <- diag(exp(fit$par))%*%solve(fit$hessian)%*%diag(exp(fit$par))
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nllh.final = fit$value)
  names(out$par.ests) <- c("alpha","gamma","theta")
  names(out$par.ses) <- c("alpha","gamma","theta")
  class(out) <- "ltburr"
  return(out)
}
