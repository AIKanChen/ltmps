#' @title  Left-truncated Log-normal distribution
#' @description Cumulative probability, quantiles, density generation from left-tructed Log-normal distribution.
#' @usage plnorm_LT(q,mu,sigma,u)
#' @usage qlnorm_LT(p,mu,sigma,u)
#' @usage dlnorm_LT(x,mu,sigma,u)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param mu Mean parameter.
#' @param sigma Scale parameter.
#' @param u The left truncated point.
#' @import stats
#' @return Probability (plnorm_LT), quantile (qlnorm_LT), density (dlnorm_LT) for the left-truncated Log-normal distribution.
#' @export
#'
#' @examples
#' plnorm_LT(q=10,mu=0,sigma=1,u=5)
#' qlnorm_LT(p=0.1,mu=0,sigma=1,u=5)
#' dlnorm_LT(x=6,mu=0,sigma=1,u=5)



plnorm_LT = function(q,mu,sigma,u){
  if(min(q)<u){stop('The quantile is less than left truncated point')}
  (plnorm(q,mu,sigma)-plnorm(u,mu,sigma))/
    (1-plnorm(u,mu,sigma))
}
