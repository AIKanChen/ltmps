#' @title  Left-truncated Weibull distribution
#' @description Cumulative probability, quantiles, density function of left-truncated Weibull distribution.
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


pweibull_LT = function(q,shape,scale,u){
  if(min(q)<u){stop('The quantile is less than left truncated point')}
  (pweibull(q,shape,scale)-pweibull(u,shape,scale))/
    (1-pweibull(u,shape,scale))
}
