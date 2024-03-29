% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aoi.R
\name{dAoI}
\alias{dAoI}
\title{Compute the density function for AoI}
\usage{
dAoI(alpha, a, FoIpar, tau = 0, hhat = 1, r = 1/200)
}
\arguments{
\item{alpha}{the age of an infection, \eqn{\alpha}}

\item{a}{host age}

\item{FoIpar}{a compound \link{list} to compute \eqn{h_\tau(a)}}

\item{tau}{the host cohort's birthday}

\item{hhat}{scaling parameter for \link{FoI}}

\item{r}{the clearance rate for a simple infection}
}
\value{
a \link{numeric} vector of length(alpha)
}
\description{
The density of the AoI is given by
\deqn{f_A(\alpha | a, \tau, h) = \frac{\int_0^a z(\alpha, a)}{m_\tau(a)}}
}
