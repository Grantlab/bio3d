read.cif <- function(file, maxlines = -1, multi=FALSE,
                     rm.insert=FALSE, rm.alt=TRUE, verbose=TRUE) {

    cl <- match.call()
    warning("helix/sheet records will not be parsed in this version of the code")
    
    if(missing(file)) {
        stop("read.cif: please specify a CIF 'file' for reading")
    }
    
    if(!is.logical(multi)) {
        stop("read.cif: 'multi' must be logical TRUE/FALSE")
    }
    
    ##- Check if file exists locally or on-line
    putfile <- NULL
    if(substr(file,1,4)=="http") {
        ## cpp function can not read from http
        putfile <- tempfile(fileext=".cif")
        rt <- try(download.file(file, putfile, quiet = !verbose), silent=TRUE)
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
            cat("  Note: Accessing on-line CIF file\n")
            file <- get.pdb(file, path=tempdir(), quiet=TRUE, format="cif")
        }
        else {
            stop("No input CIF file found: check filename")
        }
    }
    
    ## parse CIF file with cpp function
    pdb <- .read_cif(file, maxlines=maxlines, multi=multi)
    
    ## remove temp file if we downloaded it above
    if(!is.null(putfile)) {
        file.remove(putfile)
    }

    #if(verbose)
    #    cat(" ", pdb$header, "\n")
    #pdb$header <- NULL
    
    if(!is.null(pdb$error))
        stop(paste("Error in reading CIF file", file))
    else
        class(pdb) <- c("pdb") ## this should be cif perhaps?
    
    if(pdb$models > 1)
        pdb$xyz <- matrix(pdb$xyz, nrow=pdb$models, byrow=TRUE)
    
    pdb$xyz <- as.xyz(pdb$xyz)
    
    ## set empty strings to NA
    pdb$atom[pdb$atom==""] <- NA
    pdb$atom[pdb$atom=="?"] <- NA
    pdb$atom[pdb$atom=="."] <- NA
    pdb$models <- NULL


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
    
    ## set call
    pdb$call <- cl
    
    return(pdb)
}

