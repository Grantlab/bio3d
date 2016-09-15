read.pdb <- function(file, maxlines = -1, multi=FALSE, rm.insert=FALSE, rm.alt=TRUE,
                     ATOM.only=FALSE, hex=FALSE, verbose=TRUE) {
  
    cl <- match.call()
    
    if(missing(file)) {
        stop("read.pdb: please specify a PDB 'file' for reading")
    }
    
    if(!is.logical(multi)) {
        stop("read.pdb: 'multi' must be logical TRUE/FALSE")
    }
    
    ##- Check if file exists locally or on-line
    putfile <- NULL
    if(substr(file,1,4)=="http") {
        ## cpp function can not read from http
        putfile <- tempfile(fileext=".pdb")
        rt <- try(download.file(file, putfile, method='internal', quiet = !verbose), silent=TRUE)
        if(inherits(rt, "try-error")) {
            file.remove(putfile)
            stop("File not found at provided URL")
        }
        else {
            file <- putfile
        }
    }
    
    toread <- file.exists(file)
    if(toread & basename(file) != file) {
        file <- normalizePath(file)
    }
    
    ## Check for 4 letter code and possible on-line file
    if(!toread) {
        if(nchar(file)==4) {
            cat("  Note: Accessing on-line PDB file\n")
            file <- get.pdb(file, path=tempdir(), verbose=FALSE)
        }
        else {
            stop("No input PDB file found: check filename")
        }
    }
    
    ## parse PDB file with cpp function
    pdb <- .read_pdb(file, multi=multi, hex=hex, maxlines=maxlines, atoms_only=ATOM.only)

    ## remove temp file if we downloaded it above
    if(!is.null(putfile)) {
        file.remove(putfile)
    }
    
#    if(verbose)
#        cat(" ", pdb$header, "\n")
    pdb$header <- NULL
    
    if(!is.null(pdb$error))
        stop(paste("Error in reading PDB file", file))
    else
        class(pdb) <- c("pdb", "sse")
    
    ## convert xyz to matrix
    if(pdb$models > 1)
        pdb$xyz <- matrix(pdb$xyz, nrow=pdb$models, byrow=TRUE)
    
    pdb$xyz <- as.xyz(pdb$xyz)
    pdb$models <- NULL

    ## set empty strings to NA
    pdb$atom[pdb$atom==""] <- NA
    
    ## give names to seqres
    names(pdb$seqres) <- pdb$seqres_chain
    pdb$seqres_chain <- NULL
    
    ## fix stuff that should be NULL instead of empty vectors
    if(length(pdb$helix$start) > 0) {
        ## names are set in cpp func
        ## names(pdb$helix$start) <- pdb$helix$inserti
        ## names(pdb$helix$end) <- pdb$helix$inserte
    }
    else {
        pdb$helix <- NULL
    }
    
    if(length(pdb$sheet$start) > 0) {
        ## names(pdb$sheet$start) <- pdb$sheet$inserti
        ## names(pdb$sheet$end) <- pdb$sheet$inserte
        
        ##- remove repeated records for the same strand (e.g. in 1NH0)
        pa <- paste(pdb$sheet$start, pdb$sheet$inserti, pdb$sheet$chain, sep='_')
        keep.inds <- which(!duplicated(pa))
        pdb$sheet <- lapply(pdb$sheet, '[', keep.inds)

        pdb$sheet$inserti <- NULL
    }
    else {
        pdb$sheet <- NULL
    }
    
    if(!length(pdb$seqres) > 0)
        pdb$seqres <- NULL
    
    if(!length(pdb$remark350) > 0)
        pdb$remark350 <- NULL

    ## Remove 'Alt records'
    if (rm.alt) {
        if ( sum( !is.na(pdb$atom$alt) ) > 0 ) {
            first.alt <- sort( unique(na.omit(pdb$atom$alt)) )[1]
            cat(paste("   PDB has ALT records, taking",first.alt,"only, rm.alt=TRUE\n"))
            alt.inds <- which( (pdb$atom$alt != first.alt) ) # take first alt only
            if(length(alt.inds)>0) {
                pdb$atom <- pdb$atom[-alt.inds,]
                pdb$xyz <- trim.xyz(pdb$xyz, col.inds=-atom2xyz(alt.inds))
            }
        }
    }
    
    ## Remove 'Insert records'
    if (rm.insert) {
        if ( sum( !is.na(pdb$atom$insert) ) > 0 ) {
            cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
            insert.inds <- which(!is.na(pdb$atom$insert)) # rm insert positions
            pdb$atom <- pdb$atom[-insert.inds,]
            pdb$xyz <- trim.xyz(pdb$xyz, col.inds=-atom2xyz(insert.inds))
        }
    }
    
    if(any(duplicated(pdb$atom$eleno)))
        warning("duplicated element numbers ('eleno') detected")
    
    ## construct c-alpha attribute
    ca.inds <-  atom.select.pdb(pdb, string="calpha", verbose=FALSE)
    pdb$calpha <- seq(1, nrow(pdb$atom)) %in% ca.inds$atom
    
    ##- Parse REMARK records for storing symmetry matrices to 
    ##  build biological unit by calling 'biounit()'
    remark <- .parse.pdb.remark350(pdb$remark350)
    pdb$remark350 <- NULL
    pdb$remark <- remark
    
    ## set call
    pdb$call <- cl
    
    ## finished
    return(pdb)
}


##- parse REMARK records for building biological unit ('biounit()')
.parse.pdb.remark350 <- function(x) {

    raw.lines <- x

    # How many lines of REMARK 350?
    remark350 <- grep("^REMARK\\s+350", raw.lines)
    nlines <- length(remark350)

    # How many distinct biological unit?
    biolines <- grep("^REMARK\\s+350\\s+BIOMOLECULE", raw.lines)
    nbios <- length(biolines)

    if(nbios == 0) {
#       warning("REMARK 350 is incomplete.")
       return(NULL)
    }

    # End line number of each biological unit
    biolines2 <- c(biolines[-1], remark350[nlines])

    # How the biological unit was determined?
    method <- sapply(1:nbios, function(i) {
       author <- intersect(grep("^REMARK\\s+350\\s+AUTHOR DETERMINED BIOLOGICAL UNIT", raw.lines),
                            biolines[i]:biolines2[i])
       if(length(author) >= 1) return("AUTHOR")
       else return("SOFTWARE")
    } )
    # Get chain IDs to apply the transformation
    chain <- lapply(1:nbios, function(i) {
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines),
                            biolines[i]:biolines2[i])
       if(length(chlines) >= 1) {
          chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
          chs <- unlist(strsplit(chs, split=","))
       }
       else {
#          warning(paste("Can't determine chain IDs from REMARK 350 for biological unit",
#              i, sep=""))
          chs = NA
       }
       return(chs)
    } )
    if(any(is.na(chain))) return(NULL)

    mat <- lapply(1:nbios, function(i) {
       # Get transformation matrices
       mtlines <- intersect(grep("^REMARK\\s+350\\s+BIOMT", raw.lines),
                            biolines[i]:biolines2[i])
       # Get chain ID again: different trans matrices may be applied to different chains
       chlines <- intersect(grep("^REMARK\\s+350\\s+APPLY THE FOLLOWING TO CHAINS", raw.lines),
                            biolines[i]:biolines2[i])
       chs <- gsub("\\s*", "", sub("^.*:", "", raw.lines[chlines]))
       chs <- strsplit(chs, split=",")

       if(length(mtlines) == 0 || length(mtlines) %% 3 != 0) {
#          warning("Incomplete transformation matrix")
          mat <- NA
       }
       else {
          mat <- lapply(seq(1, length(mtlines), 3), function(j) {
             mt <- matrix(NA, 3, 4)
             for(k in 1:3) {
                vals <- sub("^REMARK\\s+350\\s+BIOMT[123]\\s*", "", raw.lines[mtlines[j+k-1]])
                vals <- strsplit(vals, split="\\s+")[[1]]
                mt[k, ] <- as.numeric(vals[-1])
             }
             mt
          } )
          chs.pos <- findInterval(mtlines[seq(1, length(mtlines), 3)], chlines)
          names(mat) <- sapply(chs[chs.pos], paste, collapse=" ") ## apply each mat to specific chains
       }
       return(mat)
    } )
    if(any(is.na(mat))) return(NULL)

    out <- list(biomat = list(num=nbios, chain=chain, mat=mat, method=method))
    return(out)
}

