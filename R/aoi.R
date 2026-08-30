
#' Compute the density function for AoI
#'
#' @description
#' The density of the AoI is given by
#' \deqn{f_A(\alpha | a, \bday, h) = \frac{\int_0^a z(\alpha, a)}{m_\bday(a)}}
#'
#' @param alpha the age of an infection, \eqn{\alpha}
#' @param a cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
dAoI = function(alpha, a, FoI_a, bday=0, hhat=1, r=1/200){
  dAoIcompute = function(alpha, a, FoI_a, bday, hhat, r){
    zda(alpha, a, FoI_a, bday, hhat, r)/meanMoI(a, FoI_a,  bday, hhat, r)
  }

  if(length(alpha)==1){
    return(dAoIcompute(alpha, a, FoI_a, bday, hhat, r))
  }else{
    return(sapply(alpha, dAoIcompute, a=a, FoI_a=FoI_a, hhat=hhat, bday=bday, r=r))
  }
}

#' Compute the distribution function for AoI
#'
#' The distribution function for the AoI is given by
#' \deqn{F_A(\alpha | a, \bday, h) = \int_0^\alpha f_A(\alpha, a, \bday | h) d \alpha}
#'
#' @inheritParams dAoI
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
pAoI = function(alpha, a, FoI_a, bday=0, hhat=1, r=1/200){
  pAoIfunction = function(alpha, a, FoI_a, bday, hhat, r){
    stats::integrate(dAoI,0,alpha,a=a,FoI_a=FoI_a,bday=bday, hhat=hhat,r=r)$value
  }
  if(length(alpha)==1) {return(pAoIfunction(alpha, a, FoI_a, bday, hhat, r))} else{
    return(sapply(alpha, pAoIfunction, a=a, FoI_a=FoI_a, hhat=hhat, bday=bday, r=r))}
}

#' Random numbers for the AoI
#'
#' @description
#' Draw random numbers for the AoI from a cohort, \eqn{\hat A_\bday(a)}
#'
#' @param N the number of observations
#' @param a cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of the AoI to return
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
rAoI = function(N, a, FoI_a, bday=0, hhat=1, r=1/200, alphamin=0){
  stopifnot(N>0)
  alpha = alphamin:a
  scdf = pAoI(alpha, a, FoI_a, bday, hhat, r)
  pdf = diff(scdf)
  pdf / sum(pdf)
  sample(alpha[-length(alpha)], N, replace=T, prob=pdf)
}

#' Compute the moments for the AoI density function for a cohort of age a
#'
#' @inheritParams meanMoI
#' @param n the moment desired
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
momentAoI = function(a, FoI_a, bday=0, hhat=1, r=1/200, n=1){

  fAda = function(a, FoI_a, bday, hhat, r0, n){
    ff = function(alpha, a, FoI_a, bday, hhat, r, n){
      alpha^n*zda(alpha, a, FoI_a, bday, hhat, r)
    }
    m =  meanMoI(a,FoI_a,bday,hhat,r)
    stats::integrate(ff, 0, a,a=a,FoI_a=FoI_a, bday=bday, hhat=hhat, r=r, n=n)$value/m
  }

  if(length(a)==1){return(fAda(a, FoI_a, bday, hhat, r, n))} else{
    sapply(a, fAda,FoI_a=FoI_a,bday=bday,hhat=hhat,r=r, n=n)}
}
