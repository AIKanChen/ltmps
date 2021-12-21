#' @title  Left-truncated Burr distribution
#' @description Cumulative probability, quantiles, density generation from left-tructed Burr distribution.
#' @usage pburr_LT(q,alpha,gamma,theta,u)
#' @usage qburr_LT(p,alpha,gamma,theta,u)
#' @usage dburr_LT(x,alpha,gamma,theta,u)
#' @param q vectors of quantiles.
#' @param p vectors of probabilities.
#' @param x vectors of values at which to evaluate density.
#' @param alpha First shape parameter, restricted positive.
#' @param gamma Second shape parameter, restricted positive.
#' @param theta Scale parameter.
#' @param u The left truncated point.
#' @importFrom actuar pburr
#' @importFrom actuar qburr
#' @importFrom actuar dburr
#' @return Probability (pburr_LT), quantile (qburr_LT), density (dburr_LT) for the left-truncated Burr distribution.
#' @export
#'
#' @examples
#' pburr_LT(q=6,alpha=1,gamma=4,theta=0.2,u=5)
#' qburr_LT(p=0.6,alpha=1,gamma=4,theta=0.2,u=5)
#' dburr_LT(x=6,alpha=1,gamma=4,theta=0.2,u=5)

pburr_LT<-function(q,alpha,gamma,theta,u){
  if(q<u){stop('The quantile is less than left truncated point')}
  (pburr(q,shape1=alpha,shape2=gamma,scale=theta)-
     pburr(u,shape1=alpha,shape2=gamma,scale=theta))/
    (1-pburr(u,shape1=alpha,shape2=gamma,scale=theta))
}
