

"align" <-
  function(fixed, mobile, fixed.inds  = NULL, mobile.inds = NULL,
           write.pdbs = TRUE, max.cycles = 10, cutoff = 0.5,
           seqaln.file="aln.fa" ) {
    
    if (missing(fixed)) 
      stop("align: must supply 'pdb' object, i.e. from 'read.pdb'")
    if (missing(mobile)) 
      stop("align: must supply 'pdb' object, i.e. from 'read.pdb'")
    
    if(class(fixed)!="pdb")
      stop("align: 'fixed' must be of type 'pdb'")
    if(class(mobile)!="pdb")
      stop("align: 'mobile' must be of type 'pdb'")
    
    ## a and b are now flipped to stay consistent with function fit.xyz
    ## if indices are provided, make new PDB entities
    if ( !is.null(fixed.inds) ) {
      if(length(fixed.inds$atom)<2)
        stop("align: insufficent atom indices for fitting")
      
      b <- NULL
      b$atom <- fixed$atom[fixed.inds$atom, ]
      b$xyz <- fixed$xyz[fixed.inds$xyz]
      b$calpha <- as.logical(b$atom[,"elety"] == "CA")
    }
    else {
      b <- fixed
    }
    
    if ( !is.null(mobile.inds) ) {
      if(length(mobile.inds$atom)<2)
        stop("align: insufficent atom indices for fitting")
      
      a <- NULL
      a$atom <- mobile$atom[mobile.inds$atom, ]
      a$xyz <- mobile$xyz[mobile.inds$xyz]
      a$calpha <- as.logical(a$atom[,"elety"] == "CA")
    }
    else {
      a <- mobile
    }

    "xyz2atom" <- function(num) {
      tmp <- seq(3, length(num), by=3)
      num[tmp]/3
    }

    "dist" <- function(c) {
      sqrt((c[1]-c[4])**2 + (c[2]-c[5])**2 + (c[3]-c[6])**2)
    }
    
    "resi.dev" <- function(xyz.a, xyz.b, cycle=1, cutoff = 0.5) {
      k <- matrix(xyz.a, ncol=3, byrow=T)
      l <- matrix(xyz.b, ncol=3, byrow=T)
      
      devs <- apply( cbind(k,l), 1, "dist")
      m <- median(devs)
      std <- sd(devs)
      
      cut <- m + (2*std)
      inds <- which( devs > cut )
      
      if ( std < cutoff ) {
        return( NULL )
      }
      else {
        print(paste("Mean: ", round(m,2),
                    " Std: ", round(std,2),
                    " Cut: ", round(cut,2)))
        print(paste( length(inds), " atoms rejected during cycle ", cycle))
        return(inds)
      }
    }
    

    ## Parse the two PDB files
    pdb.list <- NULL
    pdb.list[[1]] <- a
    pdb.list[[2]] <- b
    
    ## Sequence alignment
    s <- lapply(pdb.list, seq.pdb)
    s <- t(sapply(s, `[`, 1:max(sapply(s, length))))
    s[is.na(s)] <- "-"
    s <- seqaln(s, id = c("mobile", "fixed"), file = seqaln.file)
    gaps <- gap.inspect(s$ali)

    ## Messy code - re-write
    ## for pdb 1
    i <- 1
    pdbseq <- aa321(a$atom[a$calpha, "resid"])
    aliseq <- toupper(s$ali[i, ])
    tomatch <- gsub("X", "[A-Z]", aliseq[!is.gap(aliseq)])
    start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])), 
                           collapse = ""), text = paste(pdbseq, collapse = ""))[1]
    
    nseq <- rep(NA, length(aliseq))
    ali.res.ind <- which(!is.gap(aliseq))
    
    ali.res.ind <- ali.res.ind[1:length(pdbseq)]
    nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))

    a$atom <- cbind(a$atom, index=seq(1, nrow(a$atom)))
    ca.ali.a <- a$atom[a$calpha, ][nseq, ]
    res.nu.a <- ca.ali.a[, "resno"]
    coords.a <- as.numeric(t(ca.ali.a[, c("x", "y", "z")]))
    res.ch.a <- ca.ali.a[, "chain"]
    res.at.a <- ca.ali.a[, "eleno"]
    at.inds.a <- ca.ali.a[, "index"]

    ## for pdb 2
    i <- 2
    pdbseq <- aa321(b$atom[b$calpha, "resid"])
    aliseq <- toupper(s$ali[i, ])
    tomatch <- gsub("X", "[A-Z]", aliseq[!is.gap(aliseq)])
    start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])), 
                           collapse = ""), text = paste(pdbseq, collapse = ""))[1]
    
    nseq <- rep(NA, length(aliseq))
    ali.res.ind <- which(!is.gap(aliseq))
    
    ali.res.ind <- ali.res.ind[1:length(pdbseq)]
    nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))

    b$atom <- cbind(b$atom, index=seq(1, nrow(b$atom)))
    ca.ali.b <- b$atom[b$calpha, ][nseq, ]
    res.nu.b <- ca.ali.b[, "resno"]
    coords.b <- as.numeric(t(ca.ali.b[, c("x", "y", "z")]))
    res.ch.b <- ca.ali.b[, "chain"]
    res.at.b <- ca.ali.b[, "eleno"]
    at.inds.b <- ca.ali.b[, "index"]

    ## Fetch indices for fitting
    inds.a <- rep(FALSE, length(a$xyz))
    inds.b <- rep(FALSE, length(b$xyz))
    
    at.a <- as.numeric(at.inds.a[gaps$f.inds])
    at.b <- as.numeric(at.inds.b[gaps$f.inds])

    ## xyz indices
    inds.aa <- atom2xyz( at.a )
    inds.bb <- atom2xyz( at.b )

    ## logical vector
    inds.a[inds.aa] <- TRUE
    inds.b[inds.bb] <- TRUE

    ## Perform the initial fitting
    fit <- rot.lsq(a$xyz, b$xyz, xfit=inds.a, yfit=inds.b)
    rmsd.init <- rmsd(fit, b$xyz, a.inds=inds.aa, b.inds=inds.bb)
    print(paste("RMSD (", length(gaps$f.inds), " atoms): ", rmsd.init, sep=""))
    
    if ( write.pdbs )
      write.pdb(a, xyz=fit, file=paste("b.fit.", 0, ".pdb", sep=""))
    
    ## Refinement process 
    rmsd.all <- c(rmsd.init)
    for ( i in seq(1,max.cycles) ) {

      ## Find residues with large structural deviation
      exc <- resi.dev(fit[inds.aa], b$xyz[inds.bb], cycle = i, cutoff = cutoff)

      if ( is.null(exc) ) {
        break
      }
      else {
        ## Remove atoms for new round of fitting
        exc <- atom2xyz(exc)
        
        tmp<-seq(1,length(inds.a))
        exc.a <- tmp[which(inds.a)][exc]
        inds.a[exc.a] <- FALSE
        
        tmp <- seq(1,length(inds.b))
        exc.b <- tmp[which(inds.b)][exc]
        inds.b[exc.b] <- FALSE

        ## Fit based on new indices
        fit <- rot.lsq(a$xyz, b$xyz, xfit=inds.a, yfit=inds.b)

        if ( write.pdbs )
          write.pdb(a, xyz=fit, file=paste("b.fit.", i, ".pdb", sep=""))


        ## Calculate RMSD 
        inds.aa <- which(inds.a)
        inds.bb <- which(inds.b)
        tmp.rmsd <- rmsd(fit, b$xyz, a.inds=inds.aa, b.inds=inds.bb)

        num.resi <- length(which(inds.a))/3
        print(paste("RMSD (", num.resi, " of ", length(gaps$f.inds), " atoms): ", tmp.rmsd, sep=""))

        rmsd.all <- c(rmsd.all, tmp.rmsd)
      }
    }

    if ( write.pdbs )
      write.pdb(b, file="a.fit.pdb")

    ## Re-map indices back to the original provided PDB structures
    if ( !is.null(fixed.inds) ) {
      b.inds.new <- NULL
      b.inds.new$atom <- xyz2atom(inds.bb)

      ## This will then be the indices for the entire PDB
      ## can be used to fit the entire complex
      full.pdbs.inds <- fixed.inds$atom[b.inds.new$atom]
      
      ## map back to indices for the entire PDB given
      atoms <- fixed$atom[full.pdbs.inds, ]
      b.inds.new$atom <- full.pdbs.inds
      b.inds.new$xyz <- atom2xyz(b.inds.new$atom)
    }
    else {
      b.inds.new <- NULL
      b.inds.new$atom <- xyz2atom(inds.bb)
      b.inds.new$xyz <- inds.bb
    }
    
    if ( !is.null(mobile.inds) ) {
      a.inds.new <- NULL
      a.inds.new$atom <- xyz2atom(inds.aa)
      full.pdbs.inds <- mobile.inds$atom[a.inds.new$atom]
      atoms <- mobile$atom[full.pdbs.inds, ]
      a.inds.new$atom <- full.pdbs.inds
      a.inds.new$xyz <- atom2xyz(a.inds.new$atom)
    }
    else {
      a.inds.new <- NULL
      a.inds.new$atom <- xyz2atom(inds.aa)
      a.inds.new$xyz <- inds.aa
    }
    
    out <- list("a.inds"=b.inds.new, "b.inds"=a.inds.new, rmsd=rmsd.all)
    
    return(out)     
  }


    
 
