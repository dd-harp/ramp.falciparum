% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fmu-immune.R
\name{par_Fmu_W}
\alias{par_Fmu_W}
\title{Set up parameters for \link{Fmu.W}}
\usage{
par_Fmu_W(
  peak = 20,
  liver = 7,
  tildeb = 10.3,
  tildel = 2,
  Sa = 0.0033,
  Sw = 0.001
)
}
\arguments{
\item{peak}{The age of infection (in days) when parasite densities peak}

\item{liver}{The age of infection (in days) when parasites emerge from the liver}

\item{tildeb}{The maximum expected log10 parasite densities}

\item{tildel}{The minimum expected log10 parasite densities}

\item{Sa}{The decline in mu with respect to alpha}

\item{Sw}{The decline in mu with respect to immunity}
}
\value{
mean log10 parasite densities
}
\description{
Set up parameters for \link{Fmu.W}
}
