#' @title  Left-truncated Generalized Extreme Value(GEV) Distribution
#' @description Cumulative probability, quantiles, density function from left-truncated GEV distribution.
#' @usage pgev_LT(q, xi, sigma, u)
#' @usage qgev_LT(p, xi, sigma, u)
#' @usage dgev_LT(x,xi,sigma,u)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param xi Shape parameter.
#' @param sigma Scale parameter.
#' @param u The left-truncated point.
#' @importFrom evir dgev
#' @importFrom evir pgev
#' @importFrom evir qgev
#' @return Probability (pgev), quantile (qgev), density (dgev) for the left-truncated GEV distribution.
#' @export
#'
#' @examples
#' pgev_LT(7,xi=1,sigma=2,u=4)
#' qgev_LT(0.7,xi=1,sigma=2,u=4)
#' dgev_LT(7,xi=1,sigma=2,u=4)

pgev_LT<-function(q,xi,sigma,u){
  if(min(q)<u){stop('The quantile is less than left truncated point')}
  (pgev(q,xi=xi,sigma=sigma,mu=sigma/xi)-pgev(u,xi=xi,sigma=sigma,mu=sigma/xi))/
    (1-pgev(u,xi=xi,sigma=sigma,mu=sigma/xi))
}
