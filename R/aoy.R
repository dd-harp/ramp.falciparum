
#' The distribution function for the age of the youngest infection (AoY)
#'
#' @inheritParams dAoI
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
pAoY = function(alpha, a, FoI_a, bday=0, hhat=1, r=1/200){
  moi = meanMoI(a, FoI_a, bday, hhat, r)
  X = 1-exp(-moi)
  cdf = dAoI(alpha, a, FoI_a, bday, hhat, r)
  (1-exp(-moi*cdf))/X
}

#' Alternative method for computing the distribution function for the age of the youngest infection (AoY)
#' @description
#' This method computes the distribution function for the AoY by summing over AoI distributions
#'
#' @inheritParams dAoI
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
pAoY_long= function(alpha, a, FoI_a, bday=0, hhat=1, r=1/200){
  moi = meanMoI(a, FoI_a,  bday, hhat, r)
  py = dAoI(alpha, a, FoI_a, bday, hhat, r)
  ix = max(15, 5*moi)
  aoy=py*0
  for(m in 1:ix){
    aoy = aoy+stats::dpois(m,moi)*(1-(1-py)^m)
  }
  aoy = aoy/truePR(a,FoI_a,bday,hhat,r)
  return(aoy)
}


#' The density function for the age of the youngest infection (AoY)
#'
#' @inheritParams zda
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
dAoY = function(alpha, a, FoI_a, bday=0, hhat=1, r=1/200){
  # The function call
  dAoYcompute = function(alpha, a, FoI_a, bday, hhat, r){
    moi = meanMoI(a,FoI_a,bday,hhat,r)
    moi*dAoI(alpha,a,FoI_a,bday,hhat,r)*exp(-moi*pAoI(alpha,a,FoI_a,bday,hhat,r))/truePR(a,FoI_a,bday,hhat,r)
  }
  # Use sapply to call dAoYcompute multiple times
  if(length(alpha)==1) return(dAoYcompute(alpha, a,FoI_a,bday,hhat,r))
  return(sapply(alpha, dAoYcompute, a=a, FoI_a=FoI_a,hhat=hhat,bday=bday,r=r))
}


#' The random generation function for the age of the youngest infection (AoY)
#'
#' @inheritParams rAoI
#'
#' @return a [numeric] vector of length(N)
#' @export
#'
rAoY = function(N, a, FoI_a, bday=0, hhat=1, r=1/200, alphamin=0){
  minit = function(i, alphas, iix, jix){
    min(alphas[iix[i]:jix[i]])
  }
  moi = meanMoI(a, FoI_a, bday, hhat, r)
  hatm = nzPois(N, moi)
  Ny = sum(hatm)
  hatalpha = rAoI(Ny, a, FoI_a, bday, hhat, r, alphamin)
  jix = cumsum(hatm)
  iix = c(1,jix+1)
  sapply(1:N, minit, alphas=hatalpha, iix=iix, jix=jix)
}

#' Compute the moments for the AoY density function for a cohort of age a
#'
#' @inheritParams momentAoI
#'
#' @return a [numeric] vector of length(a)
#' @export
momentAoY = function(a, FoI_a, bday=0, hhat=1, r=1/200, n=1){

  ffAoYda = function(a, FoI_a, bday, hhat, r, n){
    ff = function(alpha, a, FoI_a, bday, hhat, r, n){
      alpha^n*dAoY(alpha, a, FoI_a, bday, hhat, r)
    }
    stats::integrate(ff, 0, a,a=a, FoI_a=FoI_a, bday=bday, hhat=hhat,r=r,n=n)$value
  }
  if(length(a)==1){return(ffAoYda(a, FoI_a, bday, hhat, r, n))} else{
    sapply(a, ffAoYda, FoI_a=FoI_a, bday=bday, hhat=hhat, r=r, n=n)}
}

