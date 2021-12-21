#' @title   MLE estimation of the left truncated Log-normal distribution
#' @description Apply the left-truncated MLE estimation to the Log-normal distribution.
#' @param data Empirical data that is left-truncated.
#' @param init_value  The initial value of left-truncated Log-normal distributions parameter mu and sigma.
#' @param u The left-truncated point.
#' @import stats
#' @return a \code{list} contains the estimation result of mu and sigma parameters, standard error and covariance matrix of estimated parameters.
#' @export
#'
#' @examples
#' data0=rlnorm(30000,4,2)
#' data = data0[data0>5]
#' fit = lgnorm_LTMLE(data = data,init_value = c(0.6,2.2),u=5)
#' fit$par.ests


lgnorm_LTMLE<- function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MPS estimation function
  opt_lgnorm <- function(par){
    mu=par[1];sigma=exp(par[2])
    est = -sum(log(dlnorm(dat,mu,sigma)))+
      n*log(1-plnorm(u,mu,sigma))
    return(est)
  }
  fit <- optim(c(init_value[1],log(init_value[2])),opt_lgnorm,hessian = T)

  #result
  par.ests <- c(fit$par[1],exp(fit$par[2]))
  varcov <- diag(c(1,exp(fit$par[2])))%*%solve(fit$hessian)%*%diag(c(1,exp(fit$par[2])))
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nllh.final = fit$value)
  names(out$par.ests) <- c("mu","sigma")
  names(out$par.ses) <- c("mu","sigma")
  class(out) <- "ltlognormal"
  return(out)
}
