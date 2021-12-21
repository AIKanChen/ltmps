#' @title  Left-truncated Generalized Log Moyal Distribution
#' @description Cumulative probability, quantiles, density generation from left-tructed Generalized Log Moyal distribution.
#' @usage pGlogM(q, mu, sigma)
#' @usage qGlogM(p, mu, sigma)
#' @usage dGlogM(x, mu, sigma)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param mu Shape parameter.
#' @param sigma Scale parameter.
#' @importFrom pracma erfcinv
#' @return Probability (pGlogM_LT), quantile (qGlogM_LT), density (dGlogM_LT) for the left-truncated Generalized Log Moyal distribution.
#' @export
#'
#' @examples
#'pGlogM_LT(7,mu=5,sigma = 0.1,u=5)
#'qGlogM_LT(0.2,mu=5,sigma = 0.1,u=5)
#'dGlogM_LT(7,mu=5,sigma = 0.1,u=5)

pGlogM_LT<-function(q,mu,sigma,u){
  if(min(q)<u){stop('The quantile is less than left truncated point')}
  pGlogM<-function(q, mu, sigma){
    erfc(sqrt(0.5)*(mu/q)^(1/2/sigma))
  }
  (pGlogM(q,mu,sigma)-pGlogM(u,mu,sigma))/
    (1-pGlogM(u,mu,sigma))
}



