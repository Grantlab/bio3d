"aa2mass" <-
  function(pdb, inds=NULL, mass.custom=NULL, addter=TRUE, mmtk=FALSE) {

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
      if(any(nchar(sequ)==1))
        sequ <- aa123(sequ)
      if(any(nchar(sequ)!=3))
        stop("must supply 'pdb' object or vector of amino acid residue names")
    }

    ## Define residues masses
    if(mmtk) {
      ## MMTK (for reproduction purposes!)
      w <- c( 71.079018, 157.196106, 114.104059, 114.080689, 103.143407,
             128.131048, 128.107678,  57.05203,  137.141527, 113.159985,
             113.159985, 129.18266,  131.197384, 147.177144,  97.117044,
              87.078323, 101.105312, 186.213917, 163.176449,  99.132996)
      
      aa <- c("ALA", "ARG", "ASN", "ASP", "CYS",
              "GLN", "GLU", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO",
              "SER", "THR", "TRP", "TYR", "VAL")

      mat <- data.frame(aa3=aa, aa1=aa321(aa), aaMass=w, name=NA, formula=NA)
      rownames(mat) <- aa
    }
    else  {
      ## Read data matrix
      mat.file <- system.file(paste("matrices/aa_mass.mat",sep=""), package="bio3d")
      mat <- read.table(mat.file)
      ##return(mat)
    }
    
    ## Data frame with column names:
    ## aa3, aa1, aaMass, name, formula
    if (!is.null(mass.custom)) {
      if(class(mass.custom)!="list")
        stop("'mass.custom' must be of class 'list'")
      
      new.aas <- names(mass.custom)
      if(any(duplicated(new.aas)) ||
         any(new.aas %in% mat[,"aa3"]))
        stop("duplicate residue name(s) provided")
      
      for(new.aa in new.aas) {
        nr <- data.frame(list(aa3=new.aa, aa1="X",
                              aaMass=mass.custom[[ new.aa ]],
                              name=NA, formula=NA))
        rownames(nr) <- new.aa
        mat <- rbind(mat, nr)
      }
    }

    ## Fetch mass from data frame
    wts <- mat[sequ,"aaMass"]

    ## Check for missing masses
    if(NA%in%wts) {
      inds <- which(wts%in%NA)
      unknown <- paste(unique(sequ[inds]), collapse=" ")
      stop(paste("Unknown aminoacid identifier: ", unknown, sep=""))
    }
    
    if(addter) {
      wts[1] <- wts[1] + atom2mass("H")
      wts[length(wts)] <- wts[length(wts)] + atom2mass("O") + atom2mass("H")
    }
    return(wts)
  }
