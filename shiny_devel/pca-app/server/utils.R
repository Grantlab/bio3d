
## first 4 chars should be upper
## any chainId should remain untouched - see e.g. PDB ID 3R1C
format_pdbids <- function(acc, casefmt=toupper) {
  mysplit <- function(x) {
    str <- unlist(strsplit(x, "_"))
    if(length(str)>1)
      paste(casefmt(str[1]), "_", str[2], sep="")
    else
      casefmt(str[1])
  }

  out <- unlist(lapply(acc, mysplit))
  return(out)
}

## pdb2lum_A.ent.gz --> 2lum_A
## note that chain ID could be longer than 1 char
pdbfilename2label <- function(fn) {
  fn <- basename(fn)
  parsefn <- function(x) {
    a <- unlist(strsplit(x, "_"))
    b <- unlist(strsplit(a, "\\."))
    paste(substr(a[1], 4, 7), b[2], sep="_")
  }
  return(unlist(lapply(fn, parsefn)))
}

##- Remove leading and trailing spaces from character strings
trim.character <- function(s) {
  s <- sub("^ +", "", s)
  s <- sub(" +$", "", s)
  s[(s=="")]<-""
  return(s)
}

##- generates a random string, e.g. 11ca4f6ea112
randstr <- function() {
  return(basename(tempfile(pattern="")))
}

## should be moved to bio3d
## determines SSE data for a pdbs object
"pdbs2sse" <- function(pdbs, ind=1, rm.gaps=FALSE, exefile="dssp") {
  ind <- ind[1]
  if(file.exists(pdbs$id[ind]))
    id <- pdbs$id[ind]
  #else if(file.exists(rownames(pdbs$ali)[ind]))
  #  id <- rownames(pdbs$ali)[ind]

  sse.aln <- NULL
  pdb.ref <- try(read.pdb(id), silent=TRUE)

  if(inherits(pdb.ref, "try-error"))
    pdb.ref <- try(read.pdb(substr(basename(id), 1, 4)), silent=TRUE)

  gaps.res <- gap.inspect(pdbs$ali)

  sse.ref <- NULL
  if(!inherits(pdb.ref, "try-error"))
    sse.ref <- try(dssp(pdb.ref, exefile=exefile), silent=TRUE)

  if(!inherits(sse.ref, "try-error") & !inherits(pdb.ref, "try-error")) {
    if(rm.gaps) {
      resid <- paste0(pdbs$resno[ind, gaps.res$f.inds], pdbs$chain[ind, gaps.res$f.inds])
    }
    else {
      resid <- paste0(pdbs$resno[ind, ], pdbs$chain[ind, ])
    }

    ## Helices
    if(length(sse.ref$helix$start) > 0) {
      resid.helix <- unbound(sse.ref$helix$start, sse.ref$helix$end)
      resid.helix <- paste0(resid.helix, rep(sse.ref$helix$chain, sse.ref$helix$length))
      inds        <- which(resid %in% resid.helix)

      ## inds points now to the position in the alignment where the helices are
      new.sse <- bounds( seq(1, length(resid))[inds] )
      if(length(new.sse) > 0) {
        sse.aln$helix$start  <- new.sse[,"start"]
        sse.aln$helix$end    <- new.sse[,"end"]
        sse.aln$helix$length <- new.sse[,"length"]
      }
    }

    ## Sheets
    if(length(sse.ref$sheet$start) > 0) {
      resid.sheet <- unbound(sse.ref$sheet$start, sse.ref$sheet$end)
      resid.sheet <- paste0(resid.sheet, rep(sse.ref$sheet$chain, sse.ref$sheet$length))
      inds        <- which(resid %in% resid.sheet)
      
      new.sse <- bounds( seq(1, length(resid))[inds] )
      if(length(new.sse) > 0) {
        sse.aln$sheet$start  <- new.sse[,"start"]
        sse.aln$sheet$end    <- new.sse[,"end"]
        sse.aln$sheet$length <- new.sse[,"length"]
      }
    }

    ## SSE vector
    sse <- rep(" ", length(resid))
    if(length(sse.aln$helix$start) > 0) {
      for(i in 1:length(sse.aln$helix$start))
        sse[sse.aln$helix$start[i]:sse.aln$helix$end[i]] <- "H"
    }

    if(length(sse.aln$sheet$start) > 0) {
      for(i in 1:length(sse.aln$sheet$start))
        sse[sse.aln$sheet$start[i]:sse.aln$sheet$end[i]] <- "E"
    }

    sse.aln$sse <- sse
  }
  else {
    msg <- NULL
    if(inherits(pdb.ref, "try-error"))
      msg = c(msg, paste("File not found:", pdbs$id[1]))
    if(inherits(sse.ref, "try-error"))
      msg = c(msg, "Launching external program 'DSSP' failed")

    warning(paste("SSE cannot be drawn", msg, sep="\n  "))
  }

  return(sse.aln)
}

'col2hex' <- function(cname) {
    colMat <- col2rgb(cname)
    rgb(red = colMat[1, ]/255,
        green = colMat[2, ]/255,
        blue = colMat[3, ]/255)
}

## hacked version summary.pdb()
pdbsum <- function(object, printseq=FALSE, pdbid = NULL, chainid = NULL, ...) {

  ## Print a summary of basic PDB object features

  if( !is.pdb(object) ) {
    stop("Input should be a pdb object, as obtained from 'read.pdb()'")
  }

  ## Multi-model check and total atom count
  nmodel <- nrow(object$xyz)
  if( is.null(nmodel) ) {
    ntotal <- length(object$xyz)/3
    nmodel = 1
  } else {
    ntotal <- length(object$xyz[1,])/3
  }

  nxyz <- length(object$xyz)
  nres <- sum(object$calpha)
  chains <- unique(object$atom[,"chain"])

  all.inds <- atom.select(object, "all", verbose=FALSE)$atom
  prot.inds <- atom.select(object, "protein", verbose=FALSE)$atom
  nuc.inds <- atom.select(object, "nucleic", verbose=FALSE)$atom
  other.inds <- all.inds[! (all.inds %in% c(prot.inds, nuc.inds)) ]
  
  nprot <-length(prot.inds)
  nnuc <-length(nuc.inds)
  nresnuc <- length(unique(
    paste(object$atom$chain[nuc.inds], object$atom$insert[nuc.inds], object$atom$resno[nuc.inds], sep="-")))
                                                
  het <- object$atom[other.inds,]
  nhet.atom <- nrow(het)
  
  if(is.null(nhet.atom) | nhet.atom==0) {
    nhet.atom <- 0
    nhet.res <- 0
    hetres <- "none"
  } else { 
  	hetres.resno <- apply(het[,c("chain","resno","resid")], 1, paste, collapse=".")
  	nhet.res <- length(unique(hetres.resno))
	hetres.nres <- table(het[,c("resid")][!duplicated(hetres.resno)])
  	hetres <- paste( paste0( names(hetres.nres), " (",hetres.nres, ")"), collapse=", ")
  }

  if((nprot+nnuc+nhet.atom) != ntotal)
    warning("nPROTEIN + nNUCLEIC + nNON-PROTEIN + nNON-NUCLEIC != nTotal")

 
  #cat("\n Call:  ", paste(deparse(object$call), sep = "\n", collapse = "\n"),
                                        #      "\n", sep = "")

  #cat("\n Call: ", 
  #    "\n   pdb <- ", "read.pdb('", pdbid, "')",
  #    "\n   pdb <- ", "trim(pdb, chain = '", chainid, "')",
  #    "\n", sep = "")

  cat("\n RCSB PDB ID:", pdbid)
  cat("\n Selected Chain ID:", chainid,
      "\n")
  
  s <- paste0("\n Total Models#: ", nmodel, 
              "\n Total Atoms#: ", ntotal, ##",  XYZs#: ", nxyz,
              
              "\n Chains#: ", length(chains),
              "  (values: ", paste(chains, collapse=" "),")",

              "\n\n Protein Atoms#: ", nprot,
              "\n  (residues/Calpha atoms#: ", nres,")",
              
              "\n Nucleic acid Atoms#: ", nnuc,
              "\n  (residues/phosphate atoms#: ", nresnuc,")",

             "\n\n Non-protein/nucleic Atoms#: ", nhet.atom,
             " \n  (residues: ", nhet.res, ")",
             "\n Non-protein/nucleic resid values: \n  [ ", hetres," ]",
             "\n\n")
              
  cat(s)

  if(printseq) {
    ##protein
    if(nres>0) {
      prot.pdb <- trim.pdb(object, as.select(prot.inds))
      aa <- pdbseq(prot.pdb)
      if(!is.null(aa)) {
        if(nres > 225) {
          ## Trim long sequences before output
          aa <- c(aa[1:225], "...<cut>...", aa[(nres-3):nres])
        }
        aa <- paste("     ",  gsub(" ","", 
                                   strwrap( paste(aa,collapse=" "), 
                                           width=120, exdent=0) ), collapse="\n")
        cat("   Protein sequence:\n", aa, "\n\n", sep="")
      }
    }
    
    ## nucleic
    if(nresnuc>0) {
      na.pdb <- trim.pdb(object, as.select(nuc.inds))
      aa <- paste(object$atom$chain[nuc.inds], object$atom$insert[nuc.inds],
                  object$atom$resno[nuc.inds], object$atom$resid[nuc.inds],
                  sep="-")
      aa <- aa[!duplicated(aa)]
      aa <- unlist(lapply(strsplit(aa, "-"), function(x) x[4]))
      aa <- .aa321.na(aa)
      if(nresnuc > 225) {
        ## Trim long sequences before output
        aa <- c(aa[1:225], "...<cut>...", aa[(nresnuc-3):nresnuc])
      }
      aa <- paste("     ",  gsub(" ","", 
                                 strwrap( paste(aa,collapse=" "), 
                                         width=120, exdent=0) ), collapse="\n")
      cat("   Nucleic acid sequence:\n", aa, "\n\n", sep="")
    }
  }
    
  #i <- paste( attributes(object)$names, collapse=", ")
  #cat(strwrap(paste(" + attr:",i,"\n"),width=45, exdent=8), sep="\n")

  invisible( c(nmodel=nmodel, natom=ntotal, nxyz=nxyz, nchains=length(chains),  
               nprot=nprot, nprot.res=nres, nother=nhet.atom, nother.res=nhet.res) )


}


mk.pfam.tbl <- function(aa1) {
  ##-- Annotate a sequence with PFAM (online)
  ##    Used for single sequence input annotation.
  ##
  ## pdb <- read.pdb('5p21')
  ## aa <- pdbseq(pdb)
  ## mk.pfam.tbl(aa)

  ##- Check for invalid residue types
  aa.protein <- c("-", "X", "A", "C", "D", "E", "F", "G", "H", "I", "K", 
                  "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  aa1.match <- aa1 %in% aa.protein
  aa1.valid <- ifelse(all(aa1.match), "YES", "NO")

  ##- Report on top PFAM match (with link)
  fam <- hmmer(seq=aa1, type="hmmscan", db="pfam")
  fam.txt <- paste0(fam[1,"name"], " (acc: ", fam[1,"acc"],")")
  fam.url <- paste0("http://pfam.sanger.ac.uk/family/",fam[1,"acc"])

  ##- Simple length check and exclusionn of long sequences
  aa1.len <- length(aa1)
  #if(aa1.len > 500) { stop("Sequence length beyond our limts!") }

  return( c("length"=aa1.len, "valid"=aa1.valid, "pfam"=fam.txt, 
    "pfam.url"=fam.url, "pfam.e"=fam[1,"evalue"]) )
}


vec_is_sorted <- function(v) {
  return(sum(sort(v) == v) == length(v))
}

split_height <- function(hc) {
  i <- which.max(diff(hc$height))
  split_height <- (hc$height[i] * 0.7 + hc$height[i+1] * 0.3)
  return(split_height)
}

cutreeBio3d <- function(hc, minDistance=0.1, k=NA) {
  sh <- split_height(hc)

  message(paste("k is", k))
  
  if(vec_is_sorted(hc$height) && max(diff(hc$height)) >= minDistance) {
    c <- stats::cutree(hc, h=sh)
  }
  else {
    c <- rep(1, length(hc$order))
  }

  autok <- length(unique(c))
  
  if(!is.na(k)) {
    c <- stats::cutree(hc, k=k)
  }
  
  k <- length(unique(c))
  
  return(list(grps=c, h=sh, k=k, autok=autok))
}


