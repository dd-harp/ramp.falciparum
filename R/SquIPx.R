

#' Derivatives for the SIP queuing model
#'
#' @description
#'
#' This is the basic SIP queuing model, which starts from the model \eqn{M/M/\infty}
#' but adds clearance through treatment and chemoprotection.
#' The model tracks the MoI in a cohort of humans
#' as it ages. It assumes a time- and age-dependent hazard rate for infection,
#' called the force of infection (FoI, \eqn{h_\tau(a)}). Infections do not affect
#' each other, and each one clears independently at the rate \eqn{r}.
#'
#' Let \eqn{\zeta_i} the fraction of the population with
#' MoI = i, then
#' \deqn{\frac{d\zeta_0}{da}= -h_\tau(a) \zeta_0 + r \zeta_1}
#' and for \eqn{i\geq 1}
#' \deqn{\frac{d\zeta_i}{da}= h_\tau(a) \left( \zeta_{i-1} - \zeta_i \right)  - ri \zeta_i + r(i+1)\zeta_{i+1}}
#'
#' This function computes the derivatives in a form that can be used by [deSolve::ode].
#'
#' @param a the host age
#' @param M the state variables
#' @param pars the parameters
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#'
#' @return the derivatives as a [list]
#' @seealso [solveMMinfty]
#' @export
#'
d_SIP_het_queue_da = function(a, y, pars, FoIpar){with(as.list(c(y,pars)),{
  foi = h*FoI(a, FoIpar, tau)
  ix = 1:N
  II = y[ix]
  H = sum(II) + P + S
  x = sum(II)/H

  dP = - eta*P + (foi*rho+xi)*S + (foi*rho+xi+sigma)*sum(II) - mu*P
  dS = eta*P + r*II[1] - (foi+xi)*S - mu*S

  dI = - (r*ix + sigma + xi + mu)*II
  dI[1]  = dI[1]   + foi*(1-rho)*S

  # Natural Clearance
  dI[-N] = dI[-N] + ix[-1]*r*II[-1]

  # MoI loss due to FoI
  dI[-N] = dI[-N] - foi*II[-N]
  dI[N]  = dI[N]  - foi*rho*II[N]

  # MoI gain due to FoI
  dI[-1] = dI[-1] + foi*(1-rho)*II[-N]

  dm  = foi*(1-rho)*(1-P/H) - (r + foi*rho + sigma + xi)*m

  dm2 = foi*(1-rho)*(1-P/H+2*m) - (2*r+foi*rho+sigma+xi)*m2  + r*m

  list(c(dI, dS, dP, dm, dm2))
})}
