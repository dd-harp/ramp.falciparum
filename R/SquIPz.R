#' @title SquIPz
#'
#' @description
#'
#' This is a model for the
#' dIistribution of the multiplicity of infection (MoI) in a cohort of humans
#' as it ages. It extends the queuing model [SquIP] by adding a Tweedie process to model
#' heterogeneous exposure: each infection event increases the MoI by the multiplicity of exposure (MoE), which is given by
#' a probability mass function
#' \eqn{F_\mbox{MoE}.} The model is named after its states and processes:
#'
#' @section **SquIPz**  \eqn{=}:
#'
#' **S** - **S**usceptible \eqn{+}
#'
#' **quI** - a **qu**euing process for **I**nfected \eqn{+}
#'
#' **P** - chemo-**P**rotected
#'
#' **z** - the Tweedie process distribution
#'
#'
#' This software implementation of **SquIP** has two parts:
#' + [dSquIPz] computes the derivatives
#' + [solve_SquIPz] is a wrapper around [dSquIPz] that
#'      - sets the value of a truncation parameter
#'      - sets up the initial values, times, and parameters
#'      - calls [deSolve::ode] to solve the system of ordinary differential equations
#'      - parses the outputs
#'
#' This software implementation truncates the infinite system of dIifferential equations. To
#' check the accuracy of the truncated system, we have derived hybrid variables.
#'
#'
#' @section The Model:
#'
#'
#' **Variables** ---
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
#' **Exposure** (\eqn{h}) ---
#'
#' Let \eqn{h} denote the force of infection (FoI). This model assumes
#' that each incident infection increases the MoI by \eqn{i} with probability \deqn{F_\mbox{MoE}(X=i) = z_i.}
#'
#' **Natural Clearance** ---
#' The model assumes that each clone clears independently of the others at rate \eqn{r,} so if a person has \eqn{i} clones, the MoI goes down by one at the rate \eqn{ri.}
#'
#' **Treatment and Chemoprotection** ---
#' These equations assume that individuals are treated and cured for several different reasons:
#'
#' + Everyone in the population takes drugs at a background rate \eqn{\xi}
#'
#' + Incident infections cause disease and a fraction \eqn{\rho} gets treated
#'
#' + In addition to other modes of treatment, infected individuals get treated at the higher rate \eqn{\sigma}
#'
#' After being treated and cured, individuals enter the chemoprotected class \eqn{P.} Chemoprotection is lost
#' at the rate \eqn{\eta,} and they enter the susceptible class \eqn{S.}
#'
#' **Demography** --
#' The model assumes that all individuals die at the rate \eqn{\mu(a)}, such that
#'
#'  \deqn{\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H}
#'
#' **Dynamics** ---
#' The dynamics are described by an infinite set of differential equations:
#' \deqn{
#' \begin{array}{rl}
#'  dS/da  &= \eta P + r I_1 - (h + \xi + \mu) S \\
#'  dP/da  &=  \left(h \rho + \xi \right) \left(H-P \right) + \sigma \sum_i I_i   - \left(\eta + \mu \right) P \\
#'  dI_i / da &=h (1-\rho) \left( S z_i + \sum_j  I_{i-j} z_j \right)  + (i+1) r I_{i+1} - \left(h + \sigma + r i + \xi + \mu \right) I_i\\
#' \end{array}
#' }
#'
#' In this software implementation, the infinite system of equations is truncated at \eqn{N} and
#' \deqn{dI_N/da  =  h (1-\rho) \sum_{i \geq N} \left( S  z_i + \sum_j  I_{i-j} z_j \right)  - \left(h + \sigma + r N + \xi + \mu \right) I_N}
#'
#' @section Hybrid Variables:
#'
#' Let \eqn{m_n} denote the  \eqn{n^{th}} moment of the MoI:
#' \deqn{m_n = \frac{\textstyle{1}}{\textstyle{H}}\sum_i i^n I_i}
#' and let \eqn{\tilde z_n} denote the moments of
#' \eqn{F_\mbox{MoE}:}
#' \deqn{\tilde z_n = \sum_i i^n z_i}
#'
#' This software implementation computes the dynamics of the hybrid variables \eqn{m_1} and \eqn{m_2}:
#' \deqn{
#' \begin{array}{rl}
#' d m_1 /da &=  h(1-\rho)\tilde z_1 \left(1-P/H \right) - \left(r + h \rho +\sigma + \xi\right) m_1 \\
#' d m_2 /da &=   h(1-\rho)\hat z_2 \left( 1-P/H\right)  + (r+2 h(1-\rho)\hat z) m  - \left(2r + h \rho + \sigma + \xi \right) m_2 \\
#' \end{array}
#' }
#' After solving, the output is parsed and the
#' first two moments of the distribution of the MoI are computed from \eqn{I_i}.
#' The empirical moments should match \eqn{m_1} and \eqn{m_2} unless the value of the
#' truncation parameter \eqn{N} is set too low.
#'
#' @seealso [SquIP] and [SIPm]
#'
#'@name SquIPz
NULL

#' Derivatives for [SquIP]
#'
#' @description
#' This computes the derivatives for the queuing model [SquIP].
#' It also computes the hybrid variables \eqn{I, m_1} and \eqn{m_2}
#'
#' @seealso [SquIP] & [solve_SquIP]
#'
#' This function computes the derivatives in a form that can be used by [deSolve::ode].
#'
#' @param a the host age
#' @param y a vector of state variables
#' @param pars the parameters
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#' @param Z the MoE matrix
#'
#' @return the derivatives as a [list]
#' @seealso [solveMMinfty]
#' @export
#'
dSquIPz = function(a, y, pars, FoIpar, Z){with(as.list(c(y,pars)),{
  foi = h*FoI(a, FoIpar, tau)
  ix = 1:N
  Ii = y[ix]
  H = sum(Ii) + P + S
  x = sum(Ii)/H

  dP = - eta*P + (foi*rho+xi)*S + (foi*rho+xi+sigma)*sum(Ii) - mu*P
  dS = eta*P + r*Ii[1] - (foi+xi)*S - mu*S

  dIi = - (r*ix + sigma + xi + mu)*Ii

  # Natural Clearance
  dIi[-N] = dIi[-N] + ix[-1]*r*Ii[-1]

  # MoI loss due to FoI
  dIi[-N] = dIi[-N] - foi*Ii[-N]
  dIi[N]  = dIi[N]  - foi*rho*Ii[N]

  # MoI gain due to FoI
  dIi = dIi + foi*(1-rho)*(Z %*% c(S,Ii[-N]))

  dm_1 = foi*(1-rho)*z_1*(1-P/H) - (r + foi*rho + sigma + xi)*m_1

  dm_2 = foi*(1-rho)*z_2*(1-P/H) + (r+2*foi*(1-rho)*z_1)*m_1 - (2*r+foi*rho+sigma+xi)*m_2

  list(c(dIi, dS, dP, dm_1, dm_2))
})}

#' @title Solve [SquIP]
#'
#' @description
#'
#' This function solves the model [SquIP]. It is a wrapper around the derivatives function [dSquIP] that calls [deSolve::ode].
#' The wrapper does the following:
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
#'
#' @param h the force of infection
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#' @param r the clearance rate for a simple infection
#' @param rho the fraction of incident cases that gets treted
#' @param sigma treatment rate for infected individuals
#' @param xi background drug taking
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param da The output frequency (in days)
#' @param N truncation parameter (maximum MoI): if `NULL`, then set by rule
#'
#' @return a named [list] of parsed outputs
#'
#' @seealso [SquIPz]
#' @export
solve_SquIPz = function(h, FoIpar, F_moe, tau=0,
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
  Z = make_Z(F_moe, N)
  ix = 1:N
  zi = F_moe(ix)
  z_1 = sum(ix*zi)
  z_2 = sum(ix^2*zi)
  ages = seq(0, Amax, by = da)
  parms = c(h=h,r=r,N=N,tau=tau,rho=rho,sigma=sigma,xi=xi,eta=eta,mu=mu,z_1=z_1,z_2=z_2)
  inits = c(rep(0,N), S=H, P=0, m_1=0, m_2=0)
  out = deSolve::ode(inits, times=ages, dSquIPz, parms, FoIpar=FoIpar, Z=Z)
  parsed_list <- parse_SquIP(parms, out)
  return(parsed_list)
}


#' Make the Z matrix
#'
#' @param F_moe the density distribution for a probability mass function
#' @param N the maximum MoI
#'
#' @returns a [matrix]
#' @export
#'
#' @examples
#' dd = function(x){dpois(x,.2)}
#' make_Z(dd, 20)
make_Z = function(d_pmf, N){
  Z = matrix(0, N, N)
  for(i in 1:(N-1)){
    zij = d_pmf(1:(N-i))
    zz =  c(rep(0,i-1), zij, 1-sum(zij))
    Z[,i] = zz
  }
  Z[N,N]=1
  return(Z)
}
