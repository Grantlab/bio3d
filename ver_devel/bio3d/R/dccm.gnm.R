"dccm.gnm" <- function(x, ...) {
  if (missing(x) || !inherits(x, 'gnm'))
    stop("dccm.gnm: must supply a 'gnm' object, i.e. from 'gnm.pdb()'")

  # variance-covariance matrix
  vcov <- cov.nma(x)
  
  # DCCM
  corr.mat <- vcov * 1/( sqrt(diag(vcov)) %*% t(sqrt(diag(vcov))) )
  class(corr.mat) <- c("matrix", "dccm")

  return(corr.mat)
}


"dccm.egnm" <- function(x, ...) {

  if (missing(x) || !inherits(x, 'egnm'))
  stop("dccm.egnm: must supply a 'egnm' object, i.e. from 'gnm.pdbs()'")
 
  # variance-covariance matrix
  vcovs <- cov.enma(x)

  # DCCM
  all.dccm <- apply(vcovs, 3, function(vcov) 
     vcov * 1/( sqrt(diag(vcov)) %*% t(sqrt(diag(vcov))) ) )
  all.dccm <- array(all.dccm, dim=dim(vcovs))

  avg.dccm <- rowMeans(all.dccm, dims=2)
  class(avg.dccm) <- c("matrix", "dccm")

  out <- list(all.dccm = all.dccm, avg.dccm = avg.dccm)
  return( out )

}
