
vmd.pdbs <- function(pdbs, col=NULL, as="tube", writepdbs=TRUE, file="R.vmd", ...) {
  
  allowed <- c("tube", "trace", "cartoon", "lines")
  if(!as %in% allowed) {
    stop(paste("input argument 'as' must be either of:",
               paste(allowed, collapse=", ")))
  }
  
  if(!is.null(col) & !inherits(col, "core")) {
    if(length(col) == 1) {
      allowed <- c("index", "index2", "rmsf") #, "gaps")
      if(!col %in% allowed) {
        stop(paste("input argument 'col' must be either of:",
                   paste(allowed, collapse=", ")))
      }
    }
    else {
      if(!is.numeric(col)) {
        stop("col must be a numeric vector with length equal to the number of structures in the input pdbs object")
      }
      
      if(length(col) != length(pdbs$id)) {
        stop("col must be a vector with length equal to the number of structures in input pdbs")
      }
    }
  }

  tdir <- "."
  pdbdir <- paste(tdir, "pdbs", sep="/")
  if(!file.exists(pdbdir) & writepdbs)
    dir.create(pdbdir)
  
  vmdfile <- tempfile(tmpdir=tdir, fileext=".vmd")
  ids <- basename.pdb(pdbs$id)

  ## include stuff in the b-factor column
  bf <- NULL; ch <- NULL;
    
    ## default coloring: by index of ncol(pdbs$ali)
    if(col[1] == "index2") {
        bf <- 1:ncol(pdbs$ali)/ncol(pdbs$ali)
    }
    
    if(!is.null(col)) {
        ## RMSF coloring
        if(col[1] == "rmsf") {
            bf <- rmsf(pdbs$xyz)
        }
        
        ## color by core
        if(inherits(col, "core")) {
            core <- col
            ch <- rep("A", ncol(pdbs$ali))
            ch[core$atom] = "B"
        }
    }

  if(writepdbs) {
    ## use all-atom PDBs if they exist
    if(all(file.exists(pdbs$id))) {
      allatom <- TRUE
      files <- pdbs$id
      
      ## align all-atom PDBs to pdbs$xyz
      for(i in 1:length(pdbs$id)) {
        pdb <- read.pdb(files[i])
        sele <- atom.select(pdb, "calpha")
        gaps <- is.gap(pdbs$xyz[i,])
        pdb$xyz <- fit.xyz(pdbs$xyz[i, !gaps], pdb$xyz,
                           fixed.inds = 1:length(pdbs$xyz[i, !gaps]),
                           mobile.inds = sele$xyz)
        fn <- paste0(pdbdir, "/", ids[i], ".pdb")
        
        ## store new b-factor column to PDB
        bf2 <- NULL
        if(!is.null(bf)) {
          gaps <- is.gap(pdbs$ali[i,])
          bf2 <- pdb$atom$b*0
          bf2[sele$atom] <- bf[!gaps]
        }
        
        ## store new chain column to PDB
        ch2 <- NULL
        if(!is.null(ch)) {
            gaps <- is.gap(pdbs$ali[i,])
            resid1 <- paste(pdbs$resno[i, !gaps], pdbs$chain[i, !gaps], sep="-")
            resid2 <- paste(pdb$atom$resno, pdb$atom$chain, sep="-")
            inds <- resid2 %in% resid1
            
            ch2 <- rep("A", nrow(pdb$atom))
            ch3 <- vec2resno(ch, resid2[ inds ])
            ch2[ inds ] <- ch3
        }
        
        write.pdb(pdb, b=bf2, chain=ch2, file=fn)
        files[i] <- fn
      }
    }
    else {
      ## use pdbs$xyz to build CA-atom PDBs
      allatom <- FALSE
      files <- rep(NA, length(pdbs$id))
      for(i in 1:length(pdbs$id)) {
        pdb <- pdbs2pdb(pdbs, inds=i)[[1]]
        fn <- paste0(pdbdir, "/", ids[i], ".pdb")
        
        ## store new b-factor column to PDB
        bf2 <- NULL
        if(!is.null(bf)) {
          gaps <- is.gap(pdbs$ali[i,])
          bf2 <- bf[!gaps]
        }

        ## chains
        ch2 <- NULL
        if(!is.null(ch)) {
          gaps <- is.gap(pdbs$ali[i,])
          ch2 <- ch[!gaps]
        }
        
        write.pdb(pdb=pdb, b=bf2, chain=ch2, file=fn)
        files[i] <- fn
      }
    }
    
    if(!allatom) {
      if(!as %in% c("tube", "trace")) {
        warning("'as' set to 'tube' for c-alpha only structures")
        as <- "tube"
      }
    }
  }

    colby <- "Beta"
    if(inherits(col, "core"))
        colby <- "Chain"
    else {
        if(col=="index")
            colby <- "Index"
    }
    
    
  ## write script
    lines <- c()
    
    for(i in 1:length(files)) {
        
        lines <- c(lines,
                   paste("mol new", files[i],
                         "type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all"),
                   "mol delrep 0 top", 
                   paste("mol representation", as),
                   paste("mol color", colby), 
                   "mol material EdgyShiny",
                   "mol addrep top",
                   "mol drawframes top 0 {now}")
    }
    
    cat(lines, file=vmdfile, sep="\n")
    file.copy(vmdfile, file, overwrite=TRUE)
    unlink(vmdfile)
    message(paste("VMD state file written", file))
    invisible(file)
    
}


