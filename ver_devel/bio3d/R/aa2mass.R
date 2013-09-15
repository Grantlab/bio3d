"aa2mass" <-
  function(pdb, inds=NULL, mass.custom=NULL, addter=TRUE) {
    
    if (missing(pdb))
      stop("must supply 'pdb' object or vector of amino acid residue names")
    
    if(class(pdb)=="pdb") {
      if(!is.null(inds)) {
        pdb <- trim.pdb(pdb, inds)
      }
      sequ <- pdb$atom[pdb$calpha,"resid"]
    }
    else {
      if(!is.null(inds))
        warning("'inds' has no effect when 'pdb' is vector")
      sequ <- pdb
    }
    
    ## Define residues masses
    if(FALSE) {
      ## MMTK (for reproduction purposes!)
      w <- c( 71.079018, 157.196106, 114.104059, 114.080689, 103.143407,
             128.131048, 128.107678,  57.05203,  137.141527, 113.159985,
             113.159985, 129.18266,  131.197384, 147.177144,  97.117044,
             87.078323, 101.105312, 186.213917, 163.176449,  99.132996)
    } else {
      ## Calcualted with atom2mass()
      w <- c( 71.080, 157.204, 114.108, 114.082, 103.150,
             128.134, 128.108,  57.054, 137.146, 113.158, 
             113.158, 129.184, 131.202, 147.172,  97.116,
             87.080, 101.106, 186.210, 163.172,  99.132)
    }
    
    aa <- c("ALA", "ARG", "ASN", "ASP", "CYS",
            "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO",
            "SER", "THR", "TRP", "TYR", "VAL")
    
    if (!is.null(mass.custom)) {
      if(class(mass.custom)!="list")
        stop("'mass.custom' must be of class 'list'")
      new.aas <- names(mass.custom)
      for(new.aa in new.aas) {
        aa <- c(aa, new.aa)
        w <- c(w, mass.custom[[ new.aa ]])
      }
    }
    
    "res2wt" <- function(x, w, aa) {
      ind <- which(aa==x)
      if(length(ind)==1)
        return(w[ind])
      else
        return(NA)
    }
    
    wts <- unlist(lapply(sequ, res2wt, w, aa))
    if(NA%in%wts) {
      inds <- which(wts%in%NA)
      unknown <- paste(unique(sequ[inds]), collapse=" ")
      stop(paste("Unknown residue type:", unknown))
    }

    if(addter) {
      wts[1] <- wts[1] + 1.008
      wts[length(wts)] <- wts[length(wts)] + 17.008
    }
    return(wts)
  }




