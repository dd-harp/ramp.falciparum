% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Omega-beta.R
\name{p_Omega.beta}
\alias{p_Omega.beta}
\title{Modified beta distribution, distribution function}
\usage{
\method{p_Omega}{beta}(xi, mu, bvm = 13, par_Omega = par_Omega_beta())
}
\arguments{
\item{xi}{a vector of probabilities for \eqn{log_{10}} parasite densities}

\item{mu}{the expected value for \eqn{log_{10}} parasite densities}

\item{bvm}{blood volume as \eqn{log_{10}} red blood cells}

\item{par_Omega}{parameters defining parasite densities as a function of the mu}
}
\value{
a \link{numeric} vector of length(xi)
}
\description{
Modified beta distribution, distribution function
}
