#' @title  Left-truncated Log-normal distribution
#' @description Cumulative probability, quantiles, density functions of left-truncated Log-normal distribution.
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



qlnorm_LT = function (p,mu,sigma,u){
  if(!(max(p)<=1)&(min(p)>=0)){stop("Probability must lies between 0 and 1")}
  p_u = plnorm(u,mu,sigma)
  p0 = p*(1-p_u)+p_u
  qlnorm(p0,mu,sigma)
}
