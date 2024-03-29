% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Omega-interface.R
\name{d_Omega}
\alias{d_Omega}
\title{The density function for parasite densities in a simple malaria infection}
\usage{
d_Omega(xi, mu, bvm = 13, par_Omega = par_Omega_beta())
}
\arguments{
\item{xi}{vector of quantiles for \eqn{\log_{10}} parasite densities}

\item{mu}{the expected value for \eqn{\log_{10}} parasite densities}

\item{bvm}{blood volume as \eqn{\log_{10}} red blood cells}

\item{par_Omega}{parameters to compute parasite densities as a function of mu}
}
\value{
a \link{numeric} vector of length(xi)
}
\description{
The density function for parasite densities in a simple malaria infection
}
