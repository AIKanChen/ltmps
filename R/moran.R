#' @title Moran's Log Spacing Statistic
#' @description
#' @param fitted The fitted distribution probability function
#' @param j Number of parameters to be estimated
#'
#' @return Moran's Log Spacing Statistic
#' @export
#'
#' @examples

moran<-function(fitted,j){
  n <- length(fitted)
  m=n+1
  r=-digamma(1)
  sigma_m=m*(pi^2/6-1)-0.5-1/(6*m)
  rm=m*(log(m)+r)-0.5-1/(12*m)
  c1=rm-sqrt(n/2)*sigma_m
  c2=sigma_m/sqrt(2*n)
  stat = (-sum(log(c(1,fitted)-c(fitted,0)))+j/2-c1)/c2
  list(stat=stat, p_value = 1-pchisq(stat,df=n))
}
