#' @title  Left truncated MPS estimation of the Pareto distribution
#' @description Apply the left-truncated MPS estimation to the Pareto distribution.
#' @param data Empirical data that is left-truncated.
#' @param init_value  The initial value of left-truncated Pareto distributions parameter shape and scale.
#' @param u The left truncated point.
#' @importFrom stats optim
#' @importFrom actuar ppareto
#' @importFrom actuar dpareto
#' @importFrom actuar rpareto
#' @return a \code{list} contains the estimation result of parameters shape and scale, standard error, and covariance matrix of estimated parameters.
#' @export
#'
#' @examples
#' data0=rpareto(30000,0.5,3)
#' data = data0[data0>5]
#' fit = pareto_LTMPS(data = data,init_value = c(0.6,2),u=5)
#' fit$par.ests

pareto_LTMLE = function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MPS estimation function
  opt_pareto <- function(par){
    shape=exp(par[1]);scale=par[2]
    est = -sum(log(dpareto(dat,shape,scale)))+
      n*log(1-ppareto(u,shape,scale))
    return(est)
  }
  fit <- optim(c(log(init_value[1]),init_value[2]),opt_pareto,hessian = T)

  #result
  par.ests <- c(exp(fit$par[1]),fit$par[2])
  varcov <- diag(c(exp(fit$par[1]),1))%*%solve(fit$hessian)%*%diag(c(exp(fit$par[1]),1))
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nllh.final = fit$value)
  names(out$par.ests) <- c("shape","scale")
  names(out$par.ses) <- c("shape","scale")
  class(out) <- "ltpareto"
  return(out)
}
