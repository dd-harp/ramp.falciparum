#' Compute the likelihood for a single parasite count datum
#'
#' @param count the parasite count
#' @param mu the expected value for the count, given age
#' @param sz the size parameter for a negative binomial distribution
#' @param par_O parameters [d_Omega]
#' @param q the sampling volume (in \eqn{\log_{10} microliters})
#' @param bvm the blood volume (in \eqn{\log_{10} microliters})
#' @returns the negative log likelihood for the observation
#' @export
llik_count = function(count, mu, sz, par_O, q=6, bvm=13){

  dB <- function(xi, count, mu, q, bvm, sz, par_O){
    pr_xi = d_Omega(xi, mu, bvm, par_O)
    pr_count = stats::dnbinom(count, mu=10^(xi-q), size=sz)
    return(pr_xi*pr_count)
  }

  lik <- stats::integrate(dB, 0, bvm,
                   count=count, mu=mu,
                   q=q, bvm=bvm,
                   sz=sz, par_O=par_O)$value
  -log(lik)
}


#' Compute the log likelihood for a set of parasites counts data
#'
#' @param prms A set of parameters to be fitted
#' @param sz The size parameter (not fitted)
#' @param alpha the AoI for the observations
#' @param counts the counts
#' @param par_Fmu an object to dispatch [Fmu]
#' @param par_O an object to dispath [d_Omega]
#' @param F_b a function to constrain the upper bound for fitting for [Fmu]
#' @param F_Sa a function to constrain the slope for fitting for [Fmu]
#' @param q the sampling volume (in \eqn{\log_{10} microliters})
#' @param bvm the blood volume (in \eqn{\log_{10} microliters})
#'
#' @returns the negative log likelihood for a set of observations
#' @export
llik_counts = function(prms, sz, alpha, counts, par_Fmu, par_O, F_b, F_Sa, q=6, bvm=13){

  llik_count_i = function(i, alpha, counts, sz, par_Fmu, par_O, q=6, bvm=13){
    mu = Fmu(alpha[i], 0, par_Fmu)
    llik_count(counts[i], mu, sz, par_O, q, bvm)
  }


  par_Fmu$tildeb = F_b(prms[1], bvm)
  par_Fmu$Sa = F_Sa(prms[2])

  par_O$pSig$aa=abs(prms[3])
  par_O$pSig$bb=abs(prms[4])

  lliks = sapply(1:length(alpha), llik_count_i,
                 alpha=alpha, counts=counts,
                 sz=sz, par_Fmu=par_Fmu, par_O=par_O,
                 q=q, bvm=bvm)
  sum(lliks)
}


#' Compute the MLE for parasite counts data
#'
#' @param alpha the ages of infection
#' @param counts the counts
#' @param inits the initial guesses for the parameters
#' @param sz the size parameter
#' @param F_b a function to constrain the upper bound for fitting for [Fmu]
#' @param F_Sa a function to constrain the slope for fitting for [Fmu]
#' @param q the sampling volume (in \eqn{\log_{10} microliters})
#' @param bvm the blood volume (in \eqn{\log_{10} microliters})
#'
#' @returns the fitted parameter values and the log likelihood
#' @export
counts_MLE = function(alpha, counts, inits, sz, F_b, F_Sa, q=6, bvm=13){
  par_Fmu = par_Fmu_base()
  par_O = par_Omega_beta()
  par_O$pSig$cc = 1

  out = stats::optim(inits, llik_counts, sz=sz,
              F_b = F_b, F_Sa=F_Sa,
              alpha=alpha, counts=counts,
              par_Fmu=par_Fmu, par_O=par_O,
              q=q, bvm=bvm)
  par = out$par
  mle = out$value

  return(c(F_b(par[1], 13), F_Sa(par[2]), abs(par[3]), abs(par[4]), mle))
}
