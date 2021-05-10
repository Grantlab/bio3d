"core.cmap" <-
function(pdbs,
         write.pdb = FALSE,
         outfile = "core.pdb",
         cutoff = NULL, 
         refine = FALSE, 
         ncore = NULL,
         ...) {

  ##  Quickly find core positions that have the largest number
  ##  of contact with neighboring residues 

  ncore <- setup.ncore(ncore)

  if(is.pdbs(pdbs)) {
     xyz=pdbs$xyz

     xyz.inds <- gap.inspect(xyz)$f.inds 
     res.inds <- gap.inspect(pdbs$ali)$f.inds

     pdbseq = aa123(pdbs$ali[1,]); pdbnum = pdbs$resno[1,]

  } else if(is.xyz(pdbs)) {
    xyz <- pdbs

    xyz.inds <- gap.inspect(xyz)$f.inds
    res.inds <- xyz2atom(xyz.inds)

    pdbseq = rep("ALA",length(xyz.inds)/3)
    pdbnum = c(1:(length(xyz.inds)/3))

  } else {
    stop("input 'pdbs' should either be:
             a 'pdbs' object as obtained from 'read.fasta.pdb'
             or a numeric 'xyz' matrix of aligned coordinates")
  }

  if(is.null(cutoff)) cutoff <- 1/sqrt(length(res.inds))

  cmap.default <- list(dcut=10, scut=3, pcut=1, mask.lower=FALSE, ncore=ncore)
  cmap.args <- .arg.filter( cmap.default, FUN=cmap.xyz, ...)
  cm <- do.call(cmap.xyz, c(list(xyz=xyz[, xyz.inds]), cmap.args) )
  cm[is.na(cm)] <- 0

  dm.default <- list(scut=3, mask.lower=FALSE, ncore=ncore)
  dm.args <- .arg.filter( dm.default, FUN=dm.xyz, ... )
  distmat <- do.call(dm.xyz, c(list(xyz=xyz[, xyz.inds]), dm.args) )
  distmat[is.na(distmat)] <- 0

  if(!is.na(dim(distmat)[3L])) 
    lt <- apply(distmat, 3, function(x) x[lower.tri(x)])
  else 
    lt <- as.matrix(distmat[lower.tri(distmat)])

#  mask <- apply(lt, 1, function(x) all(x==0))
#  vars <- apply(lt, 1, var)
  m <- rowMeans(lt)
  vars <- rowSums((lt - m)^2)/(ncol(lt)-1)
#  vars[vars>0] <- 1/(vars[vars>0])
  # normalize vars to solve numerical problems when some vars are very small
  vmin <- min(vars[cm[lower.tri(cm)]==1])
  vmax <- max(vars[cm[lower.tri(cm)]==1])
  if(vmax > vmin) 
     vars <- (vmax - vars) / (vmax - vmin)
  else
     vars <- 1 # ignore weight
#  vars <- 0.5 * (cos(vars * pi ) + 1) * 100

  varmat <- array(0, dim=dim(cm))
  varmat[lower.tri(varmat)] <- vars

  varmat[!cm] <- 0

  # calculate eigenvector centrality based on the weighted contact network
  ei <- eigen(varmat, symmetric=TRUE)

  if(refine) {
    # find all possible "core" positions in the top 10 PCs
    rsmall <- 0.01
    U <- ei$vector[, 1:10]
    U[abs(U) < rsmall] <- 0
    bcalc <- apply(U, 2, function(x) all(x>=0) | all(x<=0))
  
    if(sum(bcalc)==0) 
       stop('No valid eigenvector to define core atoms.')
  
    core.all <- apply(U[, bcalc, drop=FALSE], 2, function(x) {
      core <- list()
      core$atom <- res.inds[abs(x) > cutoff]
      core$xyz <- atom2xyz(core$atom)
      core
    })

    error.ellipsoid<-function(pos.xyz) {
      S<-var(pos.xyz)
      prj  <- eigen(S, symmetric = TRUE)
      prj$values[prj$values < 0 & prj$values >= -1.0E-12]<-1.0E-12
      vol<-4/3*pi*prod( sqrt( prj$values ) )
      out<-list(vol=vol, U=prj$vectors, L=prj$values)
    }
    
    if(length(core.all) > 1) {
       sv <- sapply(core.all, function(x) {
          if(length(x$atom) >=3) {
            xyz <- fit.xyz(pdbs$xyz[1, ], pdbs, x$xyz, x$xyz)
            xyz <- xyz[, xyz.inds]
            sum( sapply(1:(ncol(xyz)/3), function(i) error.ellipsoid(xyz[, atom2xyz(i)])$vol) )
          } else {
            -1
          }
       })
       sv[sv==-1] <- max(sv)
       core <- core.all[[which.min(sv)]]
    } else {
      core <- core.all[[1]]
    }

  } else {

    core <- list()
    core$atom <- res.inds[abs(ei$vector[, 1]) > cutoff]
    core$xyz <- atom2xyz(core$atom)

  }

  if(length(core$atom) < 3) {
     warning('Too few core atoms to do a proper fitting.')
     core$atom <- res.inds[order(abs(ei$vector[, 1]), decreasing=TRUE)[1:3]]
     core$xyz <- atom2xyz(core$atom)
  }

  if(write.pdb)
     write.pdb(xyz=xyz[1, core$xyz], resno=pdbnum[core$atom], 
       resid=pdbseq[core$atom], file=outfile) 
 
  class(core) <- 'select' 
  return(core)
}
