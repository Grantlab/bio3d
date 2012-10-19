"fit.xyz" <-
function(fixed,
         mobile,
         fixed.inds  = NULL,
         mobile.inds = NULL,
         verbose = FALSE,
         pdb.path = "",
         pdbext = "",
         outpath = "fitlsq/",
         het = FALSE,
         full.pdbs=FALSE,
         ...) {


  ### Addation (Mon Jul 23 17:26:16 PDT 2007)
  if( is.null(fixed.inds) && is.null(mobile.inds) ) {
    if(is.list(mobile)) {
      fixed.inds <- intersect(which(!is.gap(fixed)),
                              gap.inspect(mobile$xyz)$f.inds )
    } else {
      fixed.inds <- intersect(which(!is.gap(fixed)),
                              gap.inspect(mobile)$f.inds )
    }
    mobile.inds <- fixed.inds
    warning(paste("No fitting indices provided, using the",
                  length(fixed.inds)/3,  "non NA positions\n"))
  }
  
  
  if (is.null(fixed.inds)) fixed.inds=which(!is.gap(fixed))
  if (is.null(mobile.inds))  mobile.inds=gap.inspect(mobile)$f.inds

  if (length(fixed.inds) != length(mobile.inds))
    stop("length of 'fixed.inds' != length of 'mobile.inds'")

  if(!is.vector(fixed) || !is.numeric(fixed))
    stop("input 'fixed' should be a numeric vector")

  if(is.vector(mobile)) {   # INPUT is a single vector
    if(!is.numeric(mobile))
      stop("input 'mobile' should be numeric")

    if( any(is.na(fixed[fixed.inds])) ||
       any(is.na(mobile[mobile.inds])) ) {
      stop(" NA elements selected for fitting (check indices)")
    }
    fit <- rot.lsq(xx=mobile,
                   yy=fixed,
                   xfit=mobile.inds,
                   yfit=fixed.inds,
                   verbose=verbose)
    return(fit)
  } else {
    if(is.list(mobile)) {      # INPUT is a list object
      if(!is.numeric(mobile$xyz))
        stop("non numeric input 'mobile$xyz'")
      
      if( any(is.na(fixed[fixed.inds])) ||
         any(is.na(mobile$xyz[,mobile.inds])) ) {
        stop(" NA elements selected for fitting (check indices)")
      }

      
      fit <- t( apply(mobile$xyz, 1, rot.lsq,
                      yy = fixed,
                      xfit = mobile.inds,
                      yfit = fixed.inds,
                      verbose=verbose))

      if(full.pdbs) {        # FULL PDB fitting and output
        core.inds.atom = mobile.inds[seq(3,length(mobile.inds),by=3)]/3
        dir.create(outpath, FALSE)

        if(pdb.path=="") {
          full.files  <- file.path(paste(mobile$id, pdbext, sep=""))
        } else {
          full.files  <- file.path(pdb.path, paste(mobile$id, pdbext, sep=""))
        }

        for(i in 1:length(mobile$id)) {

###          pdb  <- read.pdb( paste(pdb.path,"/",mobile$id[i],pdbext,sep=""), ... )
          pdb  <- read.pdb( full.files[i], ... )
          res.resno  <- mobile$resno[i,core.inds.atom]
          res.chains <- mobile$chain[i,core.inds.atom]
          chains <- unique(res.chains[!is.na(res.chains)])

          if(length(chains)==0) {
            string <- paste("///",
                      paste(mobile$resno[i,core.inds.atom],collapse = ","),
                            "///CA/", sep="")
            inds <- atom.select(pdb, string,
                                verbose=verbose ,rm.insert=TRUE)$xyz

          } else {
            if(length(chains)==1) {
              string <- paste("//",chains,"/",
                              paste(res.resno, collapse = ","),
                              "///CA/", sep="")
              inds <- atom.select(pdb, string,
                                  verbose=verbose ,rm.insert=TRUE)$xyz
            } else {
              # indices for each chain
              inds <- NULL
              for(j in 1:length(chains)) {
                string <- paste("//",chains[j],"/",
                                paste(res.resno[ res.chains==chains[j] ],
                                      collapse = ","),
                                "///CA/", sep="")
                inds <- c(inds, atom.select(pdb, string,
                                            verbose=verbose,
                                            rm.insert=TRUE)$xyz)
              }
            }
          }
          pdb.xyz <- pdb$xyz
          if (het) 
            pdb.xyz <- c(pdb.xyz,
                         as.numeric(t(pdb$het[,c("x","y","z")])))

          if(length(inds) > length(fixed.inds)) {
            warning("Looks like we have a multi-chain pdb with no chain id: ignoring extra indices\n\t")
            inds <- inds[1:length(fixed.inds)]
          }
          
          xyz.fit <- rot.lsq(xx=pdb.xyz,
                             yy=fixed,
                             xfit=inds, # sort!!
                             yfit=fixed.inds)

          write.pdb(xyz = xyz.fit, pdb = pdb, het = het, 
                    file = paste(outpath, basename(mobile$id[i]),
                      "_flsq.pdb",sep = "") )

        }
      }
      return(fit)
    } else {
      if(full.pdbs)
        warning("Need 'mobile' list object for 'full.pdbs=TRUE'")
      if(is.matrix(mobile)) {       # INPUT is a matrix
        if(!is.numeric(mobile))
          stop("input 'mobile' should be numeric")

        if( any(is.na(fixed[fixed.inds])) ||
           any(is.na(mobile[,mobile.inds])) ) {
          stop("error: NA elements selected for fitting")
        }
        fit <- t( apply(mobile, 1, rot.lsq,
                        yy = fixed,
                        xfit = mobile.inds,
                        yfit = fixed.inds,
                        verbose=verbose))
        return(fit)
      }
    }
  }
}
