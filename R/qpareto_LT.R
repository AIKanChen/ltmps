#' @title  Left-truncated Pareto distribution
#' @description Cumulative probability, quantiles, density generation from left-tructed Pareto distribution.
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
#' @return Probability (ppareto_LT), quantile (qpareto_LT), density (dpareto_LT) for the left-truncated Pareto distribution.
#' @export
#'
#' @examples
#' ppareto_LT(q=8,shape=1,scale=4,u=5)
#' qpareto_LT(p=0.6,shape=1,scale=4,u=5)
#' dpareto_LT(x=8,shape=1,scale=4,u=5)

qpareto_LT = function (p,shape,scale,u){
  if(!(max(p)<=1)&(min(p)>=0)){stop("Probability must lies between 0 and 1")}
  p_u = ppareto(u,shape,scale)
  p0 = p*(1-p_u)+p_u
  qpareto(p0,shape,scale)
}

