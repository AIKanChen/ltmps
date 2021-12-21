#' @title  Left truncated MLE estimation of the Generalized Log Moyal distribution
#' @description Apply the left-truncated MLE estimation to the Generalized Log Moyal distribution.
#' @param data Empirical data that is left-truncated.
#' @param init_value  The initial value of left-truncated Generalized Log Moyal distributions parameters mean and sigma
#' @param u The left truncated point.
#' @importFrom  stats optim
#' @importFrom pracma erfc
#' @importFrom evir rgpd
#' @return a \code{list} contains the estimation result of parameters mean and sigma, standard error, and covariance matrix of estimated parameters.
#' @export
#'
#' @examples
#' rand = rgpd(10000,-0.2,5,2)
#' rand0 = rand[rand>5]
#' fit =GlogM_LTMLE(data = rand, init_value = c(5,2),u=5)
#' fit$par.ests

GlogM_LTMLE<-function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  pGlogM<-function(q, mu, sigma){
    erfc(sqrt(0.5)*(mu/q)^(1/2/sigma))
  }
  dGlogM<-function(x, mu, sigma){
    (mu/x)^(1/(2*sigma))*exp(-0.5*(mu/x)^(1/sigma))/(sqrt(2*pi)*sigma*x)
  }
  #MLE estimation function
  opt_GlogM <- function(par){
    mu=exp(par[1]);sigma=exp(par[2])
    est = -sum(log(dGlogM(dat,mu,sigma)))+
      n*log(1-pGlogM(u,mu,sigma))
    return(est)
  }
  fit <- optim(log(init_value),opt_GlogM,hessian = T)

  #result
  par.ests <- exp(fit$par)
  tmp0 <- diag(exp(fit$par))
  varcov <- tmp0%*%solve(fit$hessian)%*%tmp0
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nllh.final = fit$value)
  names(out$par.ests) <- c("mu","sigma")
  names(out$par.ses) <- c("mu","sigma")
  class(out) <- "ltGlogM"
  return(out)
}
