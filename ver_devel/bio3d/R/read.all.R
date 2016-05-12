"read.all" <-
function(aln, prefix ="", pdbext="", sel=NULL, ncore=NULL, ...) {

  ## Usage:
  ## sel <- c("N", "CA", "C", "O", "CB", "*G", "*D",  "*E", "*Z")
  ## pdbs.all <- read.all(aln, sel=sel)
 
  ncore <- setup.ncore(ncore)
 
  files  <- paste(prefix, aln$id, pdbext,sep="")

  ##cat(files,sep="\n")
  toread <- file.exists(files)

  ## check for online files
  toread[ substr(files,1,4)=="http" ] <- TRUE


  if(all(!toread))
    stop("No corresponding PDB files found")


  blank <- rep(NA, ncol(aln$ali))
  
#  for (i in 1:length(aln$id)) {
  rtn <- mclapply(1:length(aln$id), function(i) {

    cat(paste("pdb/seq:",i,"  name:", aln$id[i]),"\n")

    if(!toread[i]) {
      warning(paste("No PDB file found for seq", aln$id[i],
              ": (with filename) ",files[i]), call.=FALSE)
      coords <- rep(blank,3)
      res.nu <- blank
      res.bf <- blank
      res.ch <- blank
      res.id <- blank
      ## all atom data
      coords.all <- NULL
      elety.all <- NULL; resid.all <- NULL; resno.all <- NULL
      ##
      ##coords.all
      ##
    } else {
      pdb <- read.pdb( files[i], verbose=FALSE, ... )

      ## Currently only works for protein.
      ## Consider developing for ligand etc. in future.
      pdb <- trim(pdb, 'protein') 

      pdbseq  <- aa321(pdb$atom[pdb$calpha,"resid"])
      aliseq  <- toupper(aln$ali[i,])
      tomatch <- gsub("X","[A-Z]",aliseq[aliseq!="-"])
      
      ##-- Search for ali residues (1:15) in pdb
      start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])),collapse=""),
                           text = paste(pdbseq,collapse=""))[1]
      if (start.num == -1) {
        stop("Starting residues of sequence not found in PDB")
      }

      ##-- Numeric vec, 'nseq', for mapping aln to pdb
      nseq <- rep(NA,length(aliseq))
      ali.res.ind <- which(aliseq != "-")
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
                           pdb$atom[pdb$calpha,"resno"][nseq] )
                           #### calpha[,"resno"][nseq] ) ###- typo??
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
      coords <- as.numeric( t(ca.ali[,c("x","y","z")]) )
      res.nu <- ca.ali[, "resno"]
      res.bf <- as.numeric( ca.ali[,"b"] )
      res.ch <- ca.ali[, "chain"]
      res.id <- ca.ali[, "resid"]

      raw <- store.atom(pdb)
      if(is.null(sel)) {
        coords.all <- as.numeric( raw[c("x","y","z"),,nseq] ) 
        elety.all <- c(raw[c("elety"),,nseq])
        resid.all <- c(raw[c("resid"),,nseq])
        resno.all <- c(raw[c("resno"),,nseq])


      } else {
        coords.all <- as.numeric( raw[c("x","y","z"), sel, nseq] )
        elety.all <- c(raw[c("elety"),sel,nseq])
        resid.all <- c(raw[c("resid"),sel,nseq])
        resno.all <- c(raw[c("resno"),sel,nseq])
      } 
##      raw <- store.main(pdb)
##      b <- cbind(b, raw[,,nseq])

    } # end else 
    list(coords=coords, coords.all=coords.all, res.nu=res.nu, res.bf=res.bf,
         res.ch=res.ch, res.id=res.id, elety.all=elety.all, resid.all=resid.all,
         resno.all=resno.all)
  }, mc.cores=ncore) # end for

  coords <- do.call( rbind, unname(sapply(rtn, '[', 'coords')) )
  res.nu <- do.call( rbind, unname(sapply(rtn, '[', 'res.nu')) )
  res.bf <- do.call( rbind, unname(sapply(rtn, '[', 'res.bf')) )
  res.ch <- do.call( rbind, unname(sapply(rtn, '[', 'res.ch')) )
  res.id <- do.call( rbind, unname(sapply(rtn, '[', 'res.id')) )
  coords.all <- do.call( rbind, unname(sapply(rtn, '[', 'coords.all')) )
  elety.all <- do.call( rbind,  unname(sapply(rtn, '[', 'elety.all')) )
  resid.all <- do.call( rbind,  unname(sapply(rtn, '[', 'resid.all')) )
  resno.all <- do.call( rbind,  unname(sapply(rtn, '[', 'resno.all')) )

  rownames(aln$ali) <- aln$id
##  out<-list(xyz=coords, resno=res.nu, b=res.bf,
##            chain = res.ch, id=aln$id, ali=aln$ali)
  out<-list(xyz=coords, all=coords.all, resno=res.nu, b=res.bf,
            chain = res.ch, id=aln$id, ali=aln$ali, resid=res.id,
            all.elety=elety.all, all.resid=resid.all, all.resno=resno.all)

  atm <- rep( rep(sel,each=3), ncol(aln$ali))
  colnames(out$all) = atm
  atm <- rep( sel, ncol(aln$ali))
  colnames(out$all.elety) = atm
  colnames(out$all.resid) = atm
  colnames(out$all.resno) = atm
  
  class(out)=c("pdbs", "fasta")
  return(out)
  
}

