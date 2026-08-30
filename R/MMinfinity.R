#' @title MMinfinity
#'
#' @description
#' The queuing model \eqn{M/M/\infty} tracks the MoI in a cohort of humans
#' as it ages. It assumes a time- and age-dependent hazard rate for infection,
#' called the force of infection (FoI, \eqn{h_\bday(a)}). Infections do not affect
#' each other, and each one clears independently at the rate \eqn{r}.
#'
#' Let \eqn{\zeta_i} the fraction of the population with
#' MoI = i, then
#' \deqn{\frac{d\zeta_0}{da}= -h_\bday(a) \zeta_0 + r \zeta_1}
#' and for \eqn{i\geq 1}
#' \deqn{\frac{d\zeta_i}{da}= h_\bday(a) \left( \zeta_{i-1} - \zeta_i \right)  - ri \zeta_i + r(i+1)\zeta_{i+1}}
#'
#'@name MMinfinity
NULL


#' Derivatives for the queuing model \eqn{M/M/\infty}
#'
#' @description
#'
#' This queuing model \eqn{M/M/\infty} tracks the MoI in a cohort of humans
#' as it ages. It assumes a time- and age-dependent hazard rate for infection,
#' called the force of infection (FoI, \eqn{h_\bday(a)}). Infections do not affect
#' each other, and each one clears independently at the rate \eqn{r}.
#'
#' Let \eqn{\zeta_i} the fraction of the population with
#' MoI = i, then
#' \deqn{\frac{d\zeta_0}{da}= -h_\bday(a) \zeta_0 + r \zeta_1}
#' and for \eqn{i\geq 1}
#' \deqn{\frac{d\zeta_i}{da}= h_\bday(a) \left( \zeta_{i-1} - \zeta_i \right)  - ri \zeta_i + r(i+1)\zeta_{i+1}}
#'
#' This function computes the derivatives in a form that can be used by [deSolve::ode].
#'
#' @param a the host age
#' @param M the state variables
#' @param pars the parameters
#' @param FoI_a a cohort trace function
#'
#' @return the derivatives as a [list]
#' @keywords internal
#' @seealso [solveMMinfty]
#' @export
#'
dMoIda = function(a, M, pars, FoI_a){with(as.list(c(M,pars)),{
  foi = h*FoI_a(a, bday)
  i = 1:N
  m = i-1
  dM = 0*M-(foi + r*m)*M
  dM[-N] = dM[-N] + r*m[-1]*M[-1]
  dM[-1] = dM[-1] + foi*M[-N]
  list(c(dM))
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
#' @param FoI_a a cohort trace function
#' @param r the clearance rate for a simple infection
#' @param bday the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return a [list] with the orbits by name
#' @seealso [dMoIda]
#' @export
#'
solveMMinfty = function(h, FoI_a, r=1/200, bday=0, Amax=730, dt=1){
  tms = seq(0, Amax, by = dt)
  N = round(max(4*h/r,20))
  prms = c(h=h,r=r,N=N,bday=bday)
  inits = rep(0,N)
  inits[1]=1
  out = deSolve::ode(inits, times=tms, dMoIda, prms, FoI_a=FoI_a)
  time = out[,1]; moi = out[,-1]
  m = moi %*% c(0:(N-1))
  list(time=time, moi=moi, m=m)
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
  plot(mm, moi[t,1:N +1], type="h", xlab = "MoI", ylab = expression(M[bday](a)), main = paste ("Age = ", t, "Days"))
}
