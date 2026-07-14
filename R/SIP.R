
#' Derivatives for the queuing model \eqn{M/M/\infty}
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
d_SIPmoi_da = function(a, y, pars, FoIpar){with(as.list(c(y,pars)),{
  foi = h*FoI(a, FoIpar, tau)
  ix = 1:N
  II = y[ix]
  H = sum(II) + P + S

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

  dm = foi*(1-rho)*(H-P)/H - (r + foi*rho + sigma + xi)*m

#  sum(dI)+dS+dP -> dog
#  if(a > 365) browser()
  list(c(dI, dS, dP, dm))
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
solveSIPmoi = function(h, FoIpar, tau=0,
                       r=1/200, rho=.2,
                       sigma = 1/365, xi=1/365,
                       eta = 1/25,
                       mu = 0, H=1000,
                       Amax=730, dt=1){
  tms = seq(0, Amax, by = dt)
  N = round(max(10*h/r,20))
  prms = c(h=h,r=r,N=N,tau=tau,rho=rho,sigma=sigma,xi=xi,eta=eta,mu=mu)
  inits = c(rep(0,N), S=H, P=0, m=0)
  out = deSolve::ode(inits, times=tms, d_SIPmoi_da, prms, FoIpar=FoIpar)
  time = out[,1]; out = out[,-1]
  II = out[,1:N]
  S = out[,N+1]
  P = out[,N+2]
  m = out[,N+3]
  out = out[,-N+3]
  H = rowSums(out)
  Fr = r*II[,1]/rowSums(II)
  moi = (II/H) %*% c(1:N)
  moi2 = (II/H) %*% c(1:N)^2
  var= moi2-moi^2
  cv = var/moi
  x = (H-S-P)/H
  list(age=time, out=out, x=x, II=II, Fr = Fr, S=S, P=P, H=H, m=m, moi=moi, cv=cv, var=var, N=N)
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
