
#' The youngest of N infections, density function
#'
#' @param N the number of infections
#' @param a the age of a cohort
#' @param FoI_a a cohort trace function
#' @param hhat a local scaling parameter for the FoIs
#' @param bday the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length a + 1
#' @export
#'
dYoN = function(N, a, FoI_a, hhat=NULL, bday=0, r=1/200){
  cdf = pYoN(N, a, FoI_a, hhat, bday, r)
  pdf = diff(cdf)
  pdf/sum(pdf)
}

#' The youngest of N infections, distribution function
#'
#' @param N the MoI
#' @param a the age of a cohort
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length a + 1
#' @export
#'
pYoN = function(N, a, FoI_a, bday=0, hhat=NULL, r=1/200){
  alpha = 0:a
  py = pAoI(alpha, a, FoI_a, hhat, bday, r)
  1-(1-py)^N
}

#' The youngest of N infections, random numbers
#'
#' @param R the number of observations
#' @param N the number of infections per person
#' @param a the host cohort age
#' @param FoI_a parameters that define an FoI function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of the AoI to return
#'
#' @return a [numeric] vector of length R
#' @export
#'
rYoN = function(R, N, a, FoI_a, bday=0, hhat=NULL, r=1/200, alphamin=0){
  matrix(rAoI(R*N,a,FoI_a,hhat,bday,r,alphamin), nrow=N, ncol=R)
}
