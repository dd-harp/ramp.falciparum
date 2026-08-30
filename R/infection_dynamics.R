
#' Compute infection density in a cohort of humans, \eqn{z_\bday(\alpha, a |h)}
#'
#' @description
#' Given a function describing the FoI (\eqn{h(a,d)}), and a parameter
#' describing the clearance rate of infections (\eqn{r}),
#' the density of parasite clones of age \eqn{\alpha} distributed among a cohort of humans of
#' age \eqn{a} is \deqn{z_d(\alpha, a) = e^{-r \alpha} h(a-\alpha, d)}
#'
#' @param alpha the age of an infection, \eqn{\alpha}
#' @param a the age of the cohort
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#'
#' @export
#'
zda = function(alpha, a, FoI_a, bday=0, hhat=1, r=1/200){
  hhat*FoI_a(a-alpha, bday)*exp(-r*alpha)
}

#' The mean MoI in a host cohort of age \eqn{a}
#'
#' @description
#' The mean multiplicity of infection (MoI) is \deqn{m_\bday(a|h) = \int_0^a z_\bday(\alpha, a|h) d \alpha}
#'
#' @param a cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
meanMoI = function(a, FoI_a, bday=0, hhat=1, r=1/200){
  moi_f = function(a, FoI_a, bday, hhat,r){
    stats::integrate(zda, 0, a, a=a, FoI_a=FoI_a, bday=bday, hhat=hhat, r=r)$value
  }
  if(length(a)==1){return(moi_f(a,FoI_a,bday,hhat,r))} else{
    (return(sapply(a,moi_f,FoI_a=FoI_a,bday=bday,hhat=hhat,r=r)))}
}


#' Compute the true PR in a cohort as a function of age and exposure
#'
#' @description
#' The true PR is \deqn{p_\bday(a|h) = 1 - e^{-m_\bday(a|h)}}
#'
#' @param a the host cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
truePR = function(a, FoI_a, bday=0, hhat=1, r=1/200){
  1-exp(-meanMoI(a, FoI_a, bday, hhat, r))
}

