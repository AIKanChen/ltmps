#' @title  Left-truncated Pareto distribution
#' @description Cumulative probability, quantiles, density function of left-truncated Pareto distribution.
#' @usage ppareto_LT(q,shape,scale,u)
#' @usage qpareto_LT(p,shape,scale,u)
#' @usage dpareto_LT(x,shape,scale,u)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @param u The left truncated point.
#' @importFrom  actuar dpareto
#' @importFrom actuar ppareto
#' @importFrom actuar qpareto
#' @return Probability (ppareto_LT), quantile (qpareto_LT), density (dpareto_LT) for the left-truncated Pareto distribution.
#' @export
#'
#' @examples
#' ppareto_LT(q=8,shape=1,scale=4,u=5)
#' qpareto_LT(p=0.6,shape=1,scale=4,u=5)
#' dpareto_LT(x=8,shape=1,scale=4,u=5)


ppareto_LT = function(q,shape,scale,u){
  if(min(q)<u){stop('The quantile is less than left truncated point')}
  (ppareto(q,shape,scale)-ppareto(u,shape,scale))/
    (1-ppareto(u,shape,scale))
}
