#' @title  Left-truncated Weibull distribution
#' @description Cumulative probability, quantiles, density generation from left-tructed Weibull distribution.
#' @usage pweibull_LT(q,shape,scale,u)
#' @usage qweibull_LT(p,shape,scale,u)
#' @usage dweibull_LT(x,shape,scale,u)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param u The left truncated point.
#' @import stats
#' @return Probability (pweibull_LT), quantile (qweibull_LT), density (dweibull_LT) for the left-truncated Weibull distribution.
#' @export
#'
#' @examples
#' pweibull_LT(q=10,shape=1,scale=4,u=5)
#' qweibull_LT(p=0.6,shape=1,scale=4,u=5)
#' dweibull_LT(x=6,shape=1,scale=4,u=5)


qweibull_LT = function (p,shape,scale,u){
  if(!(max(p)<=1)&(min(p)>=0)){stop("Probability must lies between 0 and 1")}
  p_u = pweibull(u,shape,scale)
  p0 = p*(1-p_u)+p_u
  qweibull(p0,shape,scale)
}
