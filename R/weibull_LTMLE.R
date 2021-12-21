#' @title  Left truncated MPS estimation of the Weibull distribution
#' @description Apply the left-truncated MPS estimation to the Weibull distribution.
#' @param data Empirical data that is left-truncated.
#' @param init_value  The initial value of left-truncated Weibull distributions parameters shape and scale.
#' @param u The left truncated point.
#' @import stats
#' @return a \code{list} contains the estimation result of parameters shape and scale, standard error, and covariance matrix of estimated parameters.
#' @export
#'
#' @examples
#' data0 = rweibull(3000,0.5,2)
#' data = data0[data0>5]
#' fit = weibull_LTMLE(data = data,init_value = c(0.5,2),u=5)
#' fit$par.ests

weibull_LTMLE<-function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MPS estimation function
  opt_weibull <- function(par){
    shape=exp(par[1]);scale=exp(par[2])
    est = -sum(log(dweibull(dat,shape,scale)))+
      n*log(1-pweibull(u,shape,scale))
    return(est)
  }
  fit <- optim(log(init_value),opt_weibull,hessian = T)

  #result
  par.ests <- exp(fit$par)
  varcov <- diag(exp(fit$par))%*%solve(fit$hessian)%*%diag(exp(fit$par))
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nllh.final = fit$value)
  names(out$par.ests) <- c("shape","scale")
  names(out$par.ses) <- c("shape","scale")
  class(out) <- "ltweibull"
  return(out)
}
