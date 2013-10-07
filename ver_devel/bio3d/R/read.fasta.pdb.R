"read.fasta.pdb" <-
function(aln, pdb.path=".", pdbext="", ncore=1, nseg.scale=1, ...) {

  # Parallelized by multicore package (Fri Apr 26 17:58:26 EDT 2013)
  if(ncore > 1) {
     oops <- require(multicore)
     if(!oops)
        stop("Please install the multicore package from CRAN")

     options(cores = ncore)

     # Issue of serialization problem
     # Maximal number of cells of a double-precision matrix
     # that each core can serialize: (2^31-1-61)/8
     R_NCELL_LIMIT_CORE = 2.68435448e8
     R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE

     if(nseg.scale < 1) {
        warning("nseg.scale should be 1 or a larger integer\n")
        nseg.scale=1
     }
  }

  files  <- file.path(pdb.path, paste(aln$id, pdbext,sep=""))

  ##cat(files,sep="\n")
  toread <- file.exists(files)

  ## check for online files
  toread[ substr(files,1,4)=="http" ] <- TRUE

  
  if(all(!toread))
    stop("No corresponding PDB files found")

#  coords <- NULL; res.nu <- NULL
#  res.bf <- NULL; res.ch <- NULL
  blank <- rep(NA, ncol(aln$ali))
 
  mylapply <- lapply
  if(ncore > 1) mylapply <- mclapply
 
#  for (i in 1:length(aln$id)) {
  retval <- mylapply(1:length(aln$id), function(i) {
    coords <- NULL; res.nu <- NULL
    res.bf <- NULL; res.ch <- NULL
    cat(paste("pdb/seq:",i,"  name:", aln$id[i]),"\n")

    if(!toread[i]) {
      warning(paste("No PDB file found for seq", aln$id[i],
              ": (with filename) ",files[i]), call.=FALSE)
      coords <- rbind(coords, rep(blank,3))
      res.nu <- rbind(res.nu, blank)
      res.bf <- rbind(res.bf, blank)
      res.ch <- rbind(res.ch, blank)
      
    } else {
      pdb <- read.pdb( files[i], verbose=FALSE, ... )
      pdbseq  <- aa321(pdb$atom[pdb$calpha,"resid"])
      aliseq  <- toupper(aln$ali[i,])
      tomatch <- gsub("X","[A-Z]",aliseq[!is.gap(aliseq)])
      
      ##-- Search for ali residues (1:15) in pdb
      start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])),collapse=""),
                           text = paste(pdbseq,collapse=""))[1]
      if (start.num == -1) {
        stop("Starting residues of sequence not found in PDB")
      }

      ##-- Numeric vec, 'nseq', for mapping aln to pdb
      nseq <- rep(NA,length(aliseq))
      ali.res.ind <- which(!is.gap(aliseq))
      if( length(ali.res.ind) > length(pdbseq) ) {
        warning(paste(aln$id[i],
         ": sequence has more residues than PDB has Calpha's"))
        ali.res.ind <- ali.res.ind[1:length(pdbseq)] ## exclude extra
        tomatch <-  tomatch[1:length(pdbseq)]        ## terminal residues
      }
      nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))

      ##-- Check for miss-matchs
      match <- aliseq != pdbseq[nseq] 
      if ( sum(match, na.rm=TRUE) >= 1 ) {
        mismatch.ind <- which(match)
        mismatch <- cbind(aliseq, pdbseq[nseq])[mismatch.ind,]
        n.miss <- length(mismatch.ind)

        if(sum(mismatch=="X") != n.miss) { ## ignore masked X res        
          details <- rbind(aliseq, !match,
                           pdbseq[nseq],
                           pdb$atom[pdb$calpha,"resno"][nseq])
          rownames(details) = c("aliseq","match","pdbseq","pdbnum")
          msg <- paste("ERROR:", aln$id[i],
                       "alignment and pdb sequences do not match")
          cat(msg,"\n"); print(details); cat(msg,"\n")
          print( cbind(details[,mismatch.ind]) )
          stop(msg)
        }
      }
      
      ##-- Store nseq justified PDB data
      ca.ali <- pdb$atom[pdb$calpha,][nseq,]
      coords <- rbind(coords, as.numeric( t(ca.ali[,c("x","y","z")]) ))
      res.nu <- rbind(res.nu, ca.ali[, "resno"])
      res.bf <- rbind(res.bf, as.numeric( ca.ali[,"b"] ))
      res.ch <- rbind(res.ch, ca.ali[, "chain"])
#    } # end for
#  } # end else
    } # end else
    return (list(coords=coords, res.nu=res.nu, res.bf=res.bf, res.ch=res.ch))
  } ) # end mylapply
  if(ncore > 1) readChildren()
  retval <- do.call(rbind, retval)
  coords <- matrix(unlist(retval[, "coords"]), nrow=length(aln$id), byrow=TRUE)
  res.nu <- matrix(unlist(retval[, "res.nu"]), nrow=length(aln$id), byrow=TRUE)
  res.bf <- matrix(unlist(retval[, "res.bf"]), nrow=length(aln$id), byrow=TRUE)
  res.ch <- matrix(unlist(retval[, "res.ch"]), nrow=length(aln$id), byrow=TRUE)

  rownames(aln$ali) <- aln$id
  rownames(coords) <- aln$id
  rownames(res.nu) <- aln$id
  rownames(res.bf) <- aln$id
  rownames(res.ch) <- aln$id
  
  out<-list(xyz=coords, resno=res.nu, b=res.bf,
            chain = res.ch, id=aln$id, ali=aln$ali)
  class(out)="3dalign"
  return(out)
}

