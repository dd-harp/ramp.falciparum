#' @title SquIP
#'
#' @description
#'
#' A model for the
#' distribution of the multiplicity of infection (MoI) in a cohort of humans
#' as it ages.
#' The model extends a basic queuing model that has been used in malaria epidemiology before called \eqn{M/M/\infty}: it adds parameters
#' to model treatment with anti-malarial drugs, which cures infections and is followed
#' by a short period of chemo-protection. The model is named after its states:
#'
#' @section **SquIP**  \eqn{=}:
#'
#' **S** - **S**usceptible \eqn{+}
#'
#' **quI** - a **qu**euing process for **I**nfected \eqn{+}
#'
#' **P** - chemo-**P**rotected
#'
#' The master equations are presented below. The software implementation of **SquIP** has two parts:
#'
#' + [dSquIP] computes the derivatives
#'
#' + [solve_SquIP] is a wrapper around [dSquIP] that
#'      - sets the value of a truncation parameter
#'      - sets up the initial values, times, and parameters
#'      - calls [deSolve::ode] to solve the system of ordinary differential equations
#'      - parses the outputs
#'
#' This software implementation truncates the infinite system of dIifferential equations. To
#' check the accuracy of the truncated system, we have derived hybrid variables.
#'
#' @section The Dynamical System:
#'
#' **Variables**
#'
#' The independent variable is age, \eqn{a.}
#'
#' The dependent variables are:
#'
#' + \eqn{I_i} is the expected number infected with \eqn{i} clones
#'
#' + \eqn{S} is the expected number uninfected and susceptible to infection \eqn{(S=I_0)}
#'
#' + \eqn{P} is the expected number uninfected and chemo-protected
#'
#' + \eqn{H = S+P+\sum_i I_i} is expected number in the cohort
#'
#' The model assumes everyone is born susceptible
#' so \eqn{H} is passed as a parameter, and the initial conditions are set to \eqn{S(0) = H} and \eqn{I_i(0) = P(0)=0.}
#'
#'
#' **Exposure** --- \eqn{h(a)}
#'
#' Let \eqn{h(a)} denote the force of infection (FoI). This model assumes
#' that each incident infection increases the MoI by exactly one. The FoI is passed as a trace function.
#'
#' **Natural Clearance** --- \eqn{r}
#'
#' The model assumes that each clone clears independently of the others at rate \eqn{r,} so if a person has \eqn{i} clones, the MoI goes down by one at the rate \eqn{ri.}
#'
#' **Treatment and Chemoprotection** --- \eqn{\xi, \rho, \sigma, \eta}
#'
#' These equations assume that individuals are treated and cured for several different reasons:
#'
#' + Everyone in the population who is not chemoprotected takes drugs at a background rate \eqn{\xi}
#'
#' + Some incident infections cause disease, and a fraction \eqn{\rho} gets treated
#'
#' + In addition to these other modes of treatment, infected individuals get treated at the higher rate \eqn{\sigma}
#'
#' After being treated and cured, individuals enter the chemoprotected class \eqn{P.} Chemoprotection is lost
#' at the rate \eqn{\eta,} and they enter the susceptible class \eqn{S.}
#'
#' **Demography** --- \eqn{\mu}
#'
#' The model assumes that all individuals die at the rate \eqn{\mu(a)}, such that
#'
#'  \deqn{\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H}
#'
#' @section Master Equations:
#'
#' The dynamics are described by an infinite set of differential equations:
#' \deqn{
#' \begin{array}{rl}
#'  dS/da  &= \eta P + r I_1 - (h + \xi + \mu) S \\
#'  dP/da  &=  \left(h \rho + \xi \right) \left(H-P \right) + \sigma \sum_i I_i   - \left(\eta + \mu \right) P \\
#'  dI_1/ da  &=  h(1-\rho)S + 2 r I_2 - (r + h + \xi + \sigma + \mu) I_1 \\
#'  dI_i/da  &=  h(1-\rho)I_{i-1} + (i+1) r I_{i+1} - (i r + h + \xi + \sigma + \mu) I_i\\
#' \end{array}
#' }
#'
#' In this software implementation, the infinite system of equations is truncated after \eqn{i=N} and
#' \deqn{dI_N/da  =  h(1-\rho)I_{N-1} - (N r + h + \xi + \sigma + \mu) I_N}
#'
#' @section Hybrid Variables:
#'
#' To check the accuracy of the truncated system, we have derived hybrid variables for the first
#' two moments of the distribution of the MoI.
#'
#' Let \eqn{m_n} denote the \eqn{n^{th}} moment of the MoI:
#' \deqn{m_n = \frac{\textstyle{1}}{\textstyle{H}}\sum_i i^n I_i}
#'
#' This software implementation computes the dynamics of the hybrid variables \eqn{m_1} and \eqn{m_2}:
#' \deqn{
#' \begin{array}{rl}
#' d m_1 /da &=  h(1-\rho)\left(1-P/H\right) - \left(r + h \rho +\sigma + \xi\right) m_1 \\
#' d m_2 /da &=  h(1-\rho)\left( 1- P/H + 2 m_1\right)  - \left(2r + h \rho + \sigma + \xi \right) m_2 + r m_1 \\
#' \end{array}
#' }
#'
#' After solving, the output is parsed and the
#' first two moments of the distribution of the MoI are computed from \eqn{I_i}.
#' The empirical moments should match \eqn{m_1} and \eqn{m_2} unless the value of the
#' truncation parameter \eqn{N} is set too low.
#'
#' @seealso [dSquIP] & [solve_SquIP] |  The model [SquIPz] extends **SquIP.** | [SIPm] is an approximating model.
#'
#'@name SquIP
NULL

#' Derivatives for [SquIP]
#'
#' @description
#'
#' This function computes the derivatives for the queuing model [SquIP],
#' in a form that can be used by [deSolve::ode].
#'
#' It also computes the derivatives for the hybrid variables \eqn{m_1} and \eqn{m_2}
#'
#' @seealso [SquIP] & [solve_SquIP]
#'
#'
#'
#' @param a the host age
#' @param y a vector of state variables
#' @param pars the parameters
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#'
#' @return the derivatives as a [list]
#' @keywords internal
#' @seealso [solveMMinfty]
#' @export
#'
dSquIP = function(a, y, pars, FoIpar){with(as.list(c(y,pars)),{
  foi = h*FoI(a, FoIpar, tau)
  ix = 1:N
  Ii = y[ix]
  H = sum(Ii) + P + S
  x = sum(Ii)/H

  dP = - eta*P + (foi*rho+xi)*S + (foi*rho+xi+sigma)*sum(Ii) - mu*P
  dS = eta*P + r*Ii[1] - (foi+xi)*S - mu*S

  dIi = - (r*ix + sigma + xi + mu)*Ii
  dIi[1]  = dIi[1]   + foi*(1-rho)*S

  # Natural Clearance
  dIi[-N] = dIi[-N] + ix[-1]*r*Ii[-1]

  # MoI loss due to FoI
  dIi[-N] = dIi[-N] - foi*Ii[-N]
  dIi[N]  = dIi[N]  - foi*rho*Ii[N]

  # MoI gain due to FoI
  dIi[-1] = dIi[-1] + foi*(1-rho)*Ii[-N]

  dm_1  = foi*(1-rho)*(1-P/H) - (r + foi*rho + sigma + xi)*m_1

  dm_2 = foi*(1-rho)*(1-P/H+2*m_1) - (2*r+foi*rho+sigma+xi)*m_2  + r*m_1

  list(c(dIi, dS, dP, dm_1, dm_2))
})}

#' @title Solve [SquIP]
#'
#' @description
#'
#' This function solves the model [SquIP] using [deSolve::ode]. It is a wrapper around the derivatives function [dSquIP]. It does the following:
#'
#' + `N` --- set a value for the truncation parameter, if no value was passed
#'
#' + `y` --- set up the the initial values vector as a named vector
#'
#' + `times` --- set up a mesh over the independent variable \eqn{a} (age, in days)
#'
#' + `parms` --- set up a named vector with the parameter values
#'
#' After solving the equations, it parses the outputs, and returns the values of the dependent variables by name. The wrapper
#'
#' @section Parased Outputs:
#'
#' After solving the equations, the solutions are parsed and returned as a named list:
#'
#' + `age` --- the ages at which the dependent variables were output
#' + `S` --- Susceptible, a vector
#' + `P` --- Chemoprotected, a vector
#' + `Ii` --- Infected \eqn{\times} MoI, a matrix
#' + `H` --- Cohort population size
#' + `m_1` --- The hybrid variable \eqn{m_1}, mean MoI
#' + `m_2` --- The hybrid variable \eqn{m_2}, the second moment of the MoI
#' + `m1`--- The mean MoI, computed from `Ii`
#' + `m2`--- The second moment of the MoI, computed from `Ii`
#' + `x` --- The prevalence of infection: \eqn{x = (H-P-S)/H}
#' + `F_r` --- The rate of natural clearance: \eqn{r I_1}
#' + `N` --- The truncation parameter
#' + `out` --- The matrix returned by `deSolve`
#'
#' @param h the force of infection
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param rho the fraction of incident cases that gets treted
#' @param sigma treatment rate for infected individuals
#' @param xi background drug taking
#' @param eta loss of chemoprotection
#' @param mu population death rate
#' @param H cohort size
#' @param Amax The maximum runtime (in days)
#' @param da The output interval (in days)
#' @param N truncation parameter (maximum MoI): if `NULL`, then set by rule
#'
#' @return a named [list] of parsed outputs
#'
#' @seealso [SquIP]
#' @export
solve_SquIP = function(h,  FoIpar, tau=0,
                       r=1/200,
                       rho=.2,
                       sigma = 1/365,
                       xi=1/365,
                       eta = 1/25,
                       mu = 0,
                       H=1000,
                       Amax=730,
                       da=1,
                       N=NULL){
  N = ifelse(is.null(N), round(max(10*h/r,20)), N)
  ages = seq(0, Amax, by = da)
  parms = c(h=h,r=r,N=N,tau=tau,rho=rho,sigma=sigma,xi=xi,eta=eta,mu=mu)
  inits = c(rep(0,N), S=H, P=0, m_1=0, m_2=0)
  out = deSolve::ode(inits, times=ages, dSquIP, parms, FoIpar=FoIpar)
  parsed_list <- parse_SquIP(parms, out)

  return(parsed_list)
}


#' Parse SquIP
#'
#' @param pars the parameters, a named list
#' @param out the orbits
#'
#' @returns a named list
#' @keywords internal
#' @export
parse_SquIP = function(pars, out){
  vars = out[,-1]
  parsed = list()
  parsed$out = out
  with(as.list(pars), {
    ix = 1:N
    parsed$age = out[,1]
    Ii = vars[,ix]
    I = rowSums(Ii)
    S  = vars[,N+1]
    P  = vars[,N+2]
    H = I + S + P
    parsed$Ii=Ii
    parsed$I=I
    parsed$S=S
    parsed$P=P
    parsed$H=H
    parsed$m_1 = vars[,N+3]
    parsed$m_2  = vars[,N+4]
    parsed$F_r = r*Ii[,1]/I
    parsed$m1  = (Ii/H) %*% ix
    parsed$m2 = (Ii/H) %*% ix^2
    parsed$x = (H-S-P)/H
    parsed$pars = as.list(pars)
    return(parsed)
})}
