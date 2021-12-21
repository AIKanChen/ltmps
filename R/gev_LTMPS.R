#' @title  Left truncated MPS estimation of the GEV distribution
#' @description Apply the left-truncated MPS estimation to the GEV distribution
#' @param data Empirical data that is left-truncated
#' @param init_value  The initial value of left-truncated Generalised Extreme Value(GEV) distributions parameter xi and sigma
#' @param u The left truncated point
#'
#' @return a \code{list} contains the estimation result of shape and scale parameters, standard error, and covariance matrix of estimated parameters.
#' @export
#' @importFrom stats optim
#' @importFrom evir pgev
#' @importFrom evir dgev
#' @examples
#' data0=rgev(30000,xi=0.5,mu=4,sigma=2)
#' data = data0[data0>5]
#' fit = gev_LTMPS(data = data,init_value = c(0.6,2.2),u=5)
#' fit$par.ests

gev_LTMPS = function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MPS estimation function
  opt_gev<- function(par){
    xi=par[1];sigma=par[2]
    est = -sum(log(c(1,exp(-((1+xi*(dat-sigma/xi)/sigma)^(-1/xi))))-
                     exp(-((1+xi*(c(dat,u)-sigma/xi)/sigma)^(-1/xi)))))+
      (n+1)*log(1-exp(-((1+xi*(u-sigma/xi)/sigma)^(-1/xi))))
    return(est)
  }
  fit = optim(init_value,opt_gev,hessian = T)

  #result
  par.ests <- fit$par
  varcov <- solve(fit$hessian)
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nllh.final = fit$value)
  names(out$par.ests) <- c("xi", "sigma")
  names(out$par.ses) <- c("xi", "sigma")
  class(out) <- "ltgev"
  return(out)
}
