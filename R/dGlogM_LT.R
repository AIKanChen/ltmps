#' @title  Left-truncated Generalized Log Moyal Distribution
#' @description Cumulative probability, quantiles, density generation from left-tructed Generalized Log Moyal distribution.
#' @usage pgev_LT(q, xi, sigma, u)
#' @usage qgev_LT(p, xi, sigma, u)
#' @usage dGlogM(x, mu, sigma)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param mean Shape parameter.
#' @param sigma Scale parameter.
#' @param u The left-truncated point.
#' @return Probability (pGlogM_LT), quantile (qGlogM_LT), density (dGlogM_LT) for the left-truncated Generalized Log Moyal distribution.
#' @export
#'
#' @examples
#'pGlogM_LT(7,mu=5,sigma = 0.1,u=5)
#'qGlogM_LT(0.2,mu=5,sigma = 0.1,u=5)
#'dGlogM_LT(7,mu=5,sigma = 0.1,u=5)


dGlogM_LT<-function (x,mu,sigma,u){
  if(min(x)<u){stop('The value is less than left truncated point')}
  dGlogM<-function(x, mu, sigma){
    (mu/x)^(1/(2*sigma))*exp(-0.5*(mu/x)^(1/sigma))/(sqrt(2*pi)*sigma*x)
  }
  pGlogM<-function(q, mu, sigma){
    erfc(sqrt(0.5)*(mu/q)^(1/2/sigma))
  }
  dGlogM(x,mu,sigma)/(1-pGlogM(u,mu,sigma))
}


