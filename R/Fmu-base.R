
#' Compute expected log10 parasite densities, `mu`, as a function of the age of infection `alpha`
#' @description Compute mean log10 parasites for a model with no immunity
#'
#' @inheritParams Fmu
#'
#' @return mean log10 parasite densities
#' @export
#'
Fmu.base = function(alpha, W, par_Fmu){with(par_Fmu,{
  B = tildel + (tildeb-tildel)*exp(-(Sa*(alpha-peak)))
  ix = which(alpha<peak)
  if(length(ix>0)) B[ix] = tildel+(tildeb-tildel)*alpha[ix]/peak
  B
})}

#' Set up parameters for [Fmu.base]
#'
#' @param peak The age of infection (in days) when parasite densities peak
#' @param liver The age of infection (in days) when parasites emerge from the liver
#' @param tildeb The maximum expected log10 parasite densities
#' @param tildel The minimum expected log10 parasite densities
#' @param Sa The decline in mu with respect to alpha
#'
#' @return a [list]
#' @export
#'
par_Fmu_base = function(peak=20, liver=7, tildeb=10, tildel=2, Sa=0.003){
  par = list()
  class(par) <- "base"
  par$peak=peak
  par$liver=liver
  par$tildeb=tildeb
  par$tildel=tildel
  par$Sa=Sa
  par
}
