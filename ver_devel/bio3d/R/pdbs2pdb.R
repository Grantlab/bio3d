"pdbs2pdb" <- function(pdbs, inds=NULL, rm.gaps=FALSE, all.atom=FALSE, ncore=NULL) {
  if(!inherits(pdbs, "pdbs")) {
    stop("Input 'pdbs' should be of class 'pdbs', e.g. from pdbaln() or read.fasta.pdb()")
  }
  
  if(all.atom && is.null(pdbs$all))
    stop("With 'all.atom=TRUE', input 'pdbs' must be obtained from read.all()")

  if(is.null(inds))
    inds <- seq(1, length(pdbs$id))

  ncore <- setup.ncore(ncore)

  ## Set indicies
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  
  all.pdbs <- mclapply(inds, function(j) {    
    ## Temporaray file
    fname <- tempfile(fileext = "pdb")

    ## Set indices for this structure only
    f.inds <- NULL
    if(rm.gaps) {
      f.inds$res <- gaps.res$f.inds
      f.inds$pos <- gaps.pos$f.inds
    }
    else {
      f.inds$res <- which(gaps.res$bin[j,]==0)
      f.inds$pos <- atom2xyz(f.inds$res)
    }
    
    if(length(f.inds$res) > 0) {
      ## Make a temporary PDB object
      if(all.atom){
        f.inds$res <- which( (pdbs$all.grpby %in% f.inds$res) & !is.gap(pdbs$all.elety[j, ]) )
        f.inds$pos <- atom2xyz(f.inds$res)
        all.chain <- vec2resno(pdbs$chain[j, ], pdbs$all.grpby)
        if(!is.null(pdbs$insert)) {
          all.insert <- vec2resno(pdbs$insert[j, ], pdbs$all.grpby)
          insert <- all.insert[f.inds$res]
        } else {
          insert <- rep('', length(f.inds$res))
        }
        xyz <- pdbs$all[j,f.inds$pos]
        resno <- pdbs$all.resno[j,f.inds$res]
        resid <- pdbs$all.resid[j,f.inds$res]
        chain <- all.chain[f.inds$res]
        elety <- pdbs$all.elety[j,f.inds$res]
        het <- pdbs$all.hetatm[[j]]
        if(!is.null(het)) {
          xyz <- c(xyz, het$xyz)
          resno <- c(resno, het$atom[, 'resno'])
          resid <- c(resid, het$atom[, 'resid'])
          chain <- c(chain, het$atom[, 'chain'])
          elety <- c(elety, het$atom[, 'elety'])
          insert <- c(insert, rep('', nrow(het$atom)))
        }
        write.pdb(pdb=NULL, xyz=xyz, resno=resno, resid=resid, 
                  chain=chain, insert=insert, elety=elety, file=fname)
      }
      else {
        if(!is.null(pdbs$insert)) {
          insert <- pdbs$insert[j, f.inds$res]
        } else {
          insert <- rep('', length(f.inds$res))
        }
        write.pdb(pdb=NULL,
                  xyz  =pdbs$xyz[j,f.inds$pos],   resno=pdbs$resno[j,f.inds$res],
                  resid=pdbs$resid[j,f.inds$res], chain=pdbs$chain[j,f.inds$res],
                  insert=insert, file=fname)
      }
      read.pdb(fname)
    } 
    else {
      NULL
    }
  }, mc.cores=ncore )

  names(all.pdbs) <- sub(".pdb$", "", basename(pdbs$id[inds]))
  return(all.pdbs)
}
