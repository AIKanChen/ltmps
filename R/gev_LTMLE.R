#' @title   MLE estimation of the left truncated GEV distribution
#' @description Apply the left-truncated MLE estimation to the GEV distribution.
#' @param data Empirical data that is left-truncated.
#' @param init_value  The initial value of left-truncated Generalised Extreme Value(GEV) distributions parameter xi and sigma.
#' @param u The left-truncated point.
#' @importFrom stats optim
#' @importFrom evir pgev
#' @importFrom evir dgev
#' @return a \code{list} contains the estimation result of shape and scale parameters, standard error and covariance matrix of estimated parameters.
#' @export
#'
#' @examples
#' data0=rgev(30000,xi=0.5,mu=4,sigma=2)
#' data = data0[data0>5]
#' fit = gev_LTMLE(data = data,init_value = c(0.6,2.2),u=5)
#' fit$par.ests


gev_LTMLE<-function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MLE estimation function
  opt_gev<- function(par){
    xi=par[1];sigma=par[2]
    est = -sum(log(dgev(dat,xi=xi,mu=sigma/xi,sigma=sigma)))+
      n*log(1-pgev(u,xi=xi,mu=sigma/xi,sigma=sigma))
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
