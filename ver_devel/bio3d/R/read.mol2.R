"read.mol2" <-
  function (file, maxlines = -1L)
  {
    if (missing(file)) {
      stop("read.mol2: please specify a MOL2 'file' for reading")
    }
    if (!is.numeric(maxlines)) {
      stop("read.mol2: 'maxlines' must be numeric")
    }
    toread <- file.exists(file)
    if (!toread) {
      stop("No input MOL2 file found: check filename")
    }
    
    atom.format <- c("eleno", "elena", "x", "y",  "z", 
                     "elety", "resno", "resid", "charge",
                     "statbit")
    
    trim <- function(s) {
        s <- sub("^ +", "", s)
        s <- sub(" +$", "", s)
        s[(s == "")] <- NA
        s
      }

    split.line <- function(x) {
      tmp <- unlist(strsplit(x, split=" "))
      inds <- which(tmp!="")
      return(tmp[inds])
    }

    ## Read and parse mol2 file
    raw.lines <- readLines(file, n = maxlines)
    
    mol.start <- grep("@<TRIPOS>MOLECULE", raw.lines)
    atom.start <- grep("@<TRIPOS>ATOM", raw.lines)
    num.mol <- length(mol.start)

    if (!num.mol>0) {
        stop("read.mol2: mol2 file contains no molecules")
    }

    ## Fetch molecule names and info
    mol.names <- raw.lines[mol.start+1]
    mol.info <- trim( raw.lines[mol.start+2] )
    mol.info <- as.numeric(unlist(lapply(mol.info, split.line)))

    ## mol.info should contain num_atoms, num_bonds, num_subs, num_feat, num_sets
    mol.info <- matrix(mol.info, nrow=num.mol, byrow=T)
 
    num.atoms <- as.numeric(mol.info[,1])
    atom.end <- atom.start + num.atoms

    ## Build a list containing ATOM record indices
    se <- matrix(c(atom.start, atom.end), nrow=length(atom.start))
    atom.indices <- lapply(1:num.mol, function(d) seq(se[d,1]+1, se[d,2]))

    ## Check if file consist of identical molecules
    same.mol <- TRUE
    mol.first <- NULL
    
    mols <- list()
    for ( i in 1:num.mol ) {
      raw.atom <- raw.lines[ atom.indices[[i]] ]
      
      ## Split by space
      atoms.list <- lapply(raw.atom, split.line)
      ncol <- length( atoms.list[[1]] )
      atom <- matrix(unlist(atoms.list), byrow=T,
                     ncol = ncol, 
                     dimnames = list(NULL, atom.format[1:ncol]))

      ## Same molecules as the previous ones?
      mol.str <- paste(atom[,"elena"], collapse="")

      if ( i==1 ) {
        mol.first <- mol.str
      }
      else if (mol.str != mol.first) {
        same.mol <- FALSE
      }

      ## Store data
      xyz <- as.numeric(t(atom[, c("x", "y", "z")]))
      out <- list("atom" = atom, "xyz" = xyz, 
                  "info" = mol.info[i,], "name" = mol.names[i])
      mols[[i]] <- out
    }

    ## If identical molecules
    if ( length(unique(num.atoms)) == 1 && same.mol == TRUE ) {
      xyz <- t(sapply(lapply(mols, function(x) x$xyz), rbind))
      
      if ( num.mol == 1 )
        xyz <- c(xyz)
      
      out <- list("atom" = atom, "xyz" = xyz,
                  "info" = mol.info[1,], "name" = mol.names[1])
    }
    else {
      out <- mols
    }
        
    return(out)
  }
  
