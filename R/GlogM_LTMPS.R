#' @title  Left truncated MPS estimation of the Generalized Log Moyal distribution
#' @description Apply the left-truncated MPS estimation to the Generalized Log Moyal distribution.
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
#' fit =GlogM_LTMPS(data = rand, init_value = c(5,2),u=5)
#' fit$par.ests

GlogM_LTMPS<-function(data,init_value,u){
  dat = sort(data,decreasing = T)
  n=length(dat)
  u=u
  #MPS estimation function
  pGlogM<-function(q, mu, sigma){
    erfc(sqrt(0.5)*(mu/q)^(1/2/sigma))
  }
  dGlogM<-function(x, mu, sigma){
    (mu/x)^(1/(2*sigma))*exp(-0.5*(mu/x)^(1/sigma))/(sqrt(2*pi)*sigma*x)
  }
  opt_GlogM <- function(par){
    mu=exp(par[1]);sigma=exp(par[2])
    est = -sum(log(c(1, pGlogM(dat,mu,sigma))-pGlogM(c(dat,u),mu,sigma)))+
      (n+1)*log(1-pGlogM(u,mu,sigma))
    return(est)
  }
  fit <- optim(c(log(init_value)),opt_GlogM,hessian = T)

  #result
  par.ests <- exp(fit$par)
  tmp0 <- diag(exp(fit$par))
  varcov <- tmp0%*%solve(fit$hessian)%*%tmp0
  par.ses <- sqrt(diag(varcov))
  out <- list(n = n, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
              converged = fit$convergence, nmps.final = fit$value)
  names(out$par.ests) <- c("mu","sigma")
  names(out$par.ses) <- c("mu","sigma")
  class(out) <- "ltGlogM"
  return(out)
}
