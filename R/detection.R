
#' Detection of infection given parasitemia
#'
#' @param xi_mesh a mesh over \eqn{\log_{10}} parasite densities
#' @param xi_density the density distribution for parasites over xi_mesh
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param par_sample parameters that define a detection function
#'
#' @return the fraction of infected hosts that would test positive
#' @export
d_detect_mesh = function(xi_mesh, xi_density, bvm=13, par_sample = par_nb()){
  d_detect(xi_mesh, bvm, par_sample) -> pr_pos
  sum(pr_pos*xi_density)
}

#' The proportion of zero counts in
#'
#' @param alpha the age of a parasite infection
#' @param a host cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param svm the volume of a blood sample
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
FQ = function(alpha, a, FoI_a,
                 bday=0, hhat=1, r=1/200, svm = 1e6,
                 par_RBC = par_lRBC_static(),
                 par_Fmu = par_Fmu_base(),
                 par_Omega = par_Omega_beta(),
                 pWda = par_Wda_none(),
                 par_sample = par_nb()
                 ){

  frac_zero = function(mu, bvm, pO, pS){
    ff = function(x, mu, bvm, pO, pS){
      d_detect(x, bvm, pS)*d_Omega(x, mu, bvm, pO)
    }
    stats::integrate(ff, 0, bvm, mu=mu, bvm=bvm, pO=pO, pS=pS)$val
  }

  bvm = log10RBC(a, par_RBC)
  W = Wda(a, FoI_a, bday, hhat, pWda)
  mu = Fmu(alpha, W, par_Fmu)
  probs = sapply(mu, frac_zero, bvm=bvm, pO=par_Omega, pS=par_sample)
  return(probs)
}


#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [numeric] vector of length(a)
d_clone_detect = function(a, FoI_a, bday=0,
                          hhat=1, r=1/200,
                          par_RBC = par_lRBC_static(),
                          Fmu_par=par_Fmu_base(),
                          par_Omega = par_Omega_beta(),
                          pWda=par_Wda_none(),
                          par_sample = par_nb()){
  pD = function(a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, par_Omega, pWda, par_sample){
    Dx = function(x, a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, p, pWda, par_sample){
      d_clone_density(x, a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, par_Omega, pWda)*d_detect(x,a,par_sample)
    }
    hatb=log10RBC(a,par_Omega$pRBC)
    stats::integrate(Dx, 0, hatb, a=a, FoI_a=FoI_a, hhat=hhat,bday=bday,r=r,par_RBC, Fmu_par=par_RBC, Fmu_par,par_Omega=par_Omega, pWda=pWda, par_sample=par_sample)$value
  }
  if(length(a)==1) return(pD(a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, par_Omega, pWda, par_sample))
  return (sapply(a, pD, FoI_a=FoI_a, hhat=bday, hhat=bday, r=r, par_RBC, Fmu_par=par_RBC, Fmu_par, par_Omega=par_Omega, pWda=pWda, par_sample=par_sample))
}


#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param dx width of the mesh
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_parasite_detect = function(a, FoI_a, bday=0,
                               hhat=1, r=1/200, dx=0.1,
                               par_RBC = par_lRBC_static(),
                               Fmu_par = par_Fmu_base(),
                               par_Omega = par_Omega_beta(),
                               pWda=par_Wda_none(),
                               par_sample = par_nb()){

  lRBC = log10RBC(a, par_RBC)
  xi_mesh = seq(0, lRBC, by = dx)
  Bx = d_parasite_density(xi_mesh, a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, par_Omega, pWda)
  d_detect_mesh(xi_mesh, Bx, lRBC, par_sample)
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_parasite_detect_moi = function(a, FoI_a,bday=0,
                                  hhat=1,r=1/200,
                                  par_RBC = par_lRBC_static(),
                                  Fmu_par = par_Fmu_base(),
                                  par_Omega = par_Omega_beta(),
                                  pWda = par_Wda_none(),
                                  par_sample = par_nb()){
  moi = meanMoI(a, FoI_a, bday, hhat, r)
  D = d_clone_detect(a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, par_Omega, pWda, par_sample)
  1 - exp(-moi*D)
}


#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoI_a a cohort trace function
#' @param bday the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return detection probability
#' @export
d_moi_count = function(a, FoI_a,bday=0,
                       hhat=1,r=1/200,
                       par_RBC = par_lRBC_static(),
                       Fmu_par = par_Fmu_base(),
                       par_Omega = par_Omega_beta(),
                       pWda=par_Wda_none(),
                       par_sample = par_nb()){

  moi = meanMoI(a, FoI_a, bday, hhat, r)
  D = d_clone_detect(a, FoI_a, bday, hhat, r, par_RBC, Fmu_par, par_Omega, pWda, par_sample)
  moi*D/(1-exp(-moi*D))
}
