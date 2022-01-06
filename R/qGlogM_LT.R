#' @title  Left-truncated log-Moyal Distribution
#' @description Cumulative probability, quantiles, density function of left-truncated log-Moyal distribution.
#' @usage pGlogM(q, mu, sigma)
#' @usage qGlogM(p, mu, sigma)
#' @usage dGlogM(x, mu, sigma)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param mean Shape parameter.
#' @param sigma Scale parameter.
#' @param u The left-truncated point.
#' @return Probability (pGlogM_LT), quantile (qGlogM_LT), density (dGlogM_LT) for the left-truncated log-Moyal distribution.
#' @export
#'
#' @examples
#'pGlogM_LT(7,mu=5,sigma = 0.1,u=5)
#'qGlogM_LT(0.2,mu=5,sigma = 0.1,u=5)
#'dGlogM_LT(7,mu=5,sigma = 0.1,u=5)

qGlogM_LT<-function (p,mu,sigma,u){
  if(!(max(p)<=1)&(min(p)>=0)){stop("Probability must lies between 0 and 1")}
  qGlogM<-function(p, mu, sigma){
    mu*(erfcinv(p)*sqrt(2))^(-2*sigma)
  }
  pGlogM<-function(q, mu, sigma){
    erfc(sqrt(0.5)*(mu/q)^(1/2/sigma))
  }
  p_u = pGlogM(u,mu,sigma)
  p0 = p*(1-p_u)+p_u
  qGlogM(p0,mu,sigma)
}


