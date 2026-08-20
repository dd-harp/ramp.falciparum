#' @title SIPm
#'
#' @description
#' **SIPm** is an enhanced compartmental model that approximates the model **SquIPz**
#'
#' @section Variables:
#'
#' The independent variable is age, \eqn{a.}
#'
#' Three variables make up a compartmental model:
#'
#' + \eqn{S} is the expected number uninfected and susceptible to infection
#'
#' + \eqn{P} is the expected number uninfected and chemo-protected
#'
#' + \eqn{I} is the expected number infected
#'
#' Two variables track the distribution of the MoI
#'
#' + \eqn{m_1} is the first moment of the MoI
#'
#' + \eqn{m_2} is the second moment of the MoI
#'
#' @section Demography:
#'
#' We let \deqn{H = S+I+P} denote the expected population density of the cohort as it ages.
#' The model assumes that all individuals die at the rate \eqn{\mu}, such that
#'
#'  \deqn{\frac{\textstyle{dH}}{\textstyle{da}} = -\mu H}
#'
#' In the model, we assume that everyone is born susceptible
#' so \eqn{H} is passed as a parameter, and the initial conditions are set to \eqn{S(0) = H} and \eqn{I_i(0) = P(0)=0.}
#
#' @section Exposure:
#' We let \eqn{h} denote the force of infection (FoI), and we let \eqn{F_\mbox{MoE}} denote the distribution of the multiplicity of exposure (MoE).
#' + \eqn{\tilde z_1} is the first moment of \eqn{F_\mbox{MoE}}
#' + \eqn{\tilde z_2} is the second moment of \eqn{F_\mbox{MoE}}
#'
#' @section Clearance:
#' In the queuing model [SquIP], infections clear at the rate \eqn{r I_1.} If we think of **SIPm** as an approximating model, then
#' \deqn{r I_1 \approx  F_r(m_1, m_2) I}
#'
#' To get $F_r,$ we assume that the distribution of the MoI in these models is
#' approximately a zero-inflated negatively binomial distributed.
#' We assume
#' \deqn{
#' \frac{I_i} H\sim \mbox{NB}\left(\zeta=i| m, \phi \right)\left(1-\frac PH\right)
#' }
#' where `mu1=\eqn{m_1} is the mean and `size`\eqn{=\phi} (the size parameter)  is given by:
#' \deqn{\phi = \frac{m_1^2}{m_2-m_1^2-m_1}}
#' so
#' \deqn{
#' F_r(m, \phi) = \begin{cases}
#'  r & \mbox{if } m = 0 \\
#'  r \frac{
#'  \mbox{NB}\left(\zeta=1| m, \phi \right)
#'  }{1-
#'  \mbox{NB}\left(\zeta=0 | m, \phi \right)
#'  } \left(1 - \frac{P}{H} \right) & \mbox{if } m>0
#'  \end{cases}
#'  }
#'
#'  In effect, the approximation assumes \eqn{I_1/H} is close to the probability the MoI is equal to 1 in a zero-inflated negative binomial distribution. The distributions of the MoI are probably not exactly negatively binomially distributed. We use the Kullback-Liebler divergence to evaluate how well a negative binomial approximates the distributions in the queueing models.
#'
#' @section Treatment and Chemprotection:
#'
#' These equations assume that individuals are treated and cured for several different reasons:
#'
#' + Incident infections cause disease and a fraction \eqn{\rho} gets treated
#'
#' + Everyone in the population takes drugs at a background rate \eqn{\xi}
#'
#' + In addition to other modes of treatment, infected individuals get treated at the higher rate \eqn{\sigma}
#'
#' After being treated and cured, individuals enter the chemoprotected class \eqn{P.} Chemoprotection is lost
#' at the rate \eqn{\eta,} and they enter the susceptible class \eqn{S.}
#'
#' @section Differential Equations:
#'
#' **Hybrid Variables**
#' The compartmental dynamics are enhanced by computing a hybrid variable, a non-compartmental state variable describing the mean MoI (\eqn{m_1}), with dynamics:
#' \deqn{
#' \frac{\textstyle{dm_1}}{\textstyle{da}} =  h (1-\rho) \hat z \left( 1-\frac{\textstyle{P}}{\textstyle{H}} \right) - \left(r + h \rho + \sigma + \xi \right) m_1}
#'
#' We also compute the dynamics of the second moment of the distribution of the MoI:
#' \deqn{
#' \begin{array}{rl}
#' \frac{\textstyle{dm_2}}{\textstyle{da}}  = & h(1-\rho)\hat z_2 \left( 1- \frac{  \textstyle{P}}{\textstyle{H}}\right)  + (r+2 h(1-\rho)\hat z) m_1  \\[6pt]
#' &- \left(2r + h \rho + \sigma + \xi \right) m_2  \\
#' \end{array}}
#'
#' **Compartmental Variables**
#'
#' The infinite infection states in the family of queuing models specified by [SquIPz] are compressed down to a finite set of mutually exclusive and collectively exhaustive states -- \eqn{S}, \eqn{I}, and \eqn{P} -- with dynamics:
#' \deqn{
#' \begin{array}{rl}
#' \frac{\textstyle{dS}}{\textstyle{da}} &= F_r(m, \phi) I + \eta P - (h + \xi + \mu) S \\[6pt]
#' \frac{\textstyle{dI}}{\textstyle{da}} &= h \left(1-\rho\right) S  \\[6pt] & - \left( \rho h + \xi + \sigma + F_r(m, \phi) + \mu  \right) I \\[6pt]
#' \frac{\textstyle{dP}}{\textstyle{da}} &= (h \rho+ \xi)(H-P)  + \sigma I - (\eta+\mu) P \\[6pt]
#' \end{array}}
#'
#' @section Parameters:
#'
#' @seealso [SquIPz] and [SIPm]
#'
#'@name SIPm
NULL




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
dSIPm = function(a, y, pars, FoIpar){with(as.list(c(y,pars)),{
  foi = h*FoI(a, FoIpar, tau)
  #S = H - II - P
  H = S + II + P

  F_r = function(m, m2, r){
      size = m^2/(m2-m^2-m)
      r*dnbinom(1, mu=m, size=size)/(1-dnbinom(0, mu=m, size=size))*(1-P/H)
  }
  Fr = ifelse(m==0 | m2-m^2-m<0, r, F_r(m, m2, r))
  if(is.na(Fr)){
    size = m^2/(m2-m^2-m)
    browser()
    Fr=r
  }

  dP = (foi*rho+xi)*S + (foi*rho+xi+sigma)*II - (eta+mu)*P
  dS = eta*P + Fr*II - (foi+xi)*S - mu*S
  #dH = - mu*H
  dI = foi*(1-rho)*S - (foi*rho+xi+sigma+Fr)*II - mu*II

  dm  = foi*(1-rho)*(1-P/H)  - (r+foi*rho+sigma+xi)*m
  dm2 = foi*(1-rho)*(1-P/H+2*m)-(2*r+foi*rho+sigma+xi)*m2  + r*m

  #  sum(dI)+dS+dP -> dog
  #  if(a > 365) browser()
  list(c(dS, dP, dI, dm, dm2, Fr))
})}

#' Solve the queuing model \eqn{M/M/\infty}
#'
#' @description
#'
#' A wrapper to solve the queuing model \eqn{M/M/\infty} (see [dMoIda]).
#'
#' The function automatically sets the maximum MoI to be computed, and it
#' sets initial conditions. The equations are solved using [deSolve::ode] and returned at
#' regular intervals dt from age 0 up to Amax (in days).
#'
#' @param h the force of infection
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#' @param r the clearance rate for a simple infection
#' @param rho the fraction of incident cases that gets treted
#' @param sigma treatment rate for infected individuals
#' @param xi background drug taking
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return a [list] with the orbits by name
#' @seealso [dMoIda]
#' @export
solve_SIPm = function(h, FoIpar, tau=0,
                         r=1/200, rho=.2,
                         sigma = 1/365, xi=1/365,
                         eta = 1/25,
                         mu = 0, H=1000,
                         Amax=730, dt=1){
  tms = seq(0, Amax, by = dt)
  prms = c(h=h,tau=tau,r=r,rho=rho,sigma=sigma,xi=xi,eta=eta,mu=mu)
  inits = c(S=H, P=0, II=0, m=0, m2=0, dFr=0)
  out = deSolve::ode(inits, times=tms, dSIPm, prms, FoIpar=FoIpar)
  age = out[,1]; out = out[,-1]
  S = out[,1]
  P = out[,2]
  I = out[,3]
  m = out[,4]
  m2 = out[,5]
  Fr = diff(out[,6])
  S = H-P-I
  list(age=age, out=out, I=I, S=S, P=P, H=H, m=m, m2=m2, Fr=Fr)
}

#' Plot the output of [solveMMinfty]
#'
#' @param moi the mean moi
#' @param t the time
#' @param clr1 the color
#' @export
MoIDistPlot = function(moi, t, clr1 = "red"){
  N = dim(moi)[2]-2
  mm = 1:N -1
  plot(mm, moi[t,1:N +1], type="h", xlab = "MoI", ylab = expression(M[tau](a)), main = paste ("Age = ", t, "Days"))
}

