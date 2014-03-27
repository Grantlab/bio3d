"tmalign" <- function(files, ref=1) {
  ref=1

  alns <- list()
  j=1
  for(i in 1:length(files)) {
    if(!file.exists(files[i]))
      stop(paste(files[i], "does not exist"))
    
    if(i!=ref) {
      alns[[j]] <- tmalign.pair(files[i], files[ref])
      j=j+1
    }
  }

  cat("\nCombining pairwise alignments ...\n")
  
  aln <- alns[[1]]
  for ( i in 2:length(alns)) {
    aln <- .tmalignCombine(aln, alns[[i]])
  }

  out <- NULL
  out$ali <- aln
  out$id <- files

  rownames(out$ali) <- files
  return(out)
}

"tmalign.pair" <- function(a,b, exefile="TMalign", verbose=TRUE) {
  
  outfile <- tempfile()
  if(!is.pdb(a) | !is.pdb(b) ) {
    if(!file.exists(a))
      stop(paste("file", a, "not found"))
    if(!file.exists(b))
      stop(paste("file", b, "not found"))
    filea <- a
    fileb <- b
  }
  else {
    filea <- tempfile()
    fileb <- tempfile()
    
    write.pdb(a, file=filea)
    write.pdb(b, file=fileb)
  }
  
  cmd <- paste(exefile, "-A", filea, "-B", fileb, ">", outfile)

  if(verbose)
    print(cmd)
  
  system(cmd)

  rawlines <- readLines(outfile)
  seqa <- unlist(strsplit(rawlines[34], ""))
  seqb <- unlist(strsplit(rawlines[36], ""))

  seqs <- seqbind(seqa, seqb)
  return(seqs)
}



.tmalignCombine <- function(a, b) {
  ## a[1,] and b[1,] is the same sequence
  ## find b[2,] to match a[1,]
  da <- dim(a); db <- dim(b);

  new <- seqbind(a,b[1,])
  dims <- dim(new)
  new[dims[1],]=NA

  gaps.a <- gap.inspect(a)
  resnos.a <- matrix(NA, ncol=ncol(gaps.a$bin), nrow=nrow(gaps.a$bin))

  for(i in 1:nrow(gaps.a$bin)) {
      resnos.a[i, which(gaps.a$bin[i,]==0)] = seq(1, length(which(gaps.a$bin[i,]==0)))
  }

  gaps.b <- gap.inspect(b)
  resnos.b <- matrix(NA, ncol=ncol(gaps.b$bin), nrow=nrow(gaps.b$bin))
  for(i in 1:nrow(gaps.b$bin)) {
    resnos.b[i, which(gaps.b$bin[i,]==0)] = seq(1, length(which(gaps.b$bin[i,]==0)))
  }

  
  new <- seqbind(resnos.a, resnos.b[1,])
  dims <- dim(new)
  new[dims[1],]=NA


  j <- 0
  inds <-  which(!is.na(resnos.b[2, ]))
  for( i in 1:length(inds) ) {
      k=inds[i]

      res.a <- as.character(resnos.b[1, k])
      res.b <- as.character(resnos.b[2, k])
      

      if(!is.na(res.a)) {
          ind <- which(new[1, ]==res.a)
          new[dims[1], ind] = res.b

          j <- as.numeric(which(new[1, ]==res.a))
      }

      if(is.na(res.a)) {
          tmp <- matrix(NA, nrow=dims[1], ncol=1)
          tmp[dims[1], 1]=res.b

          if(j>0)
              tmp=cbind(new[, 1:j], tmp)

          new=cbind(tmp, new[, (j+1):ncol(new)])
          j=j+1
      }
  }


  for(i in 1:(nrow(new)-1)) {
    new[i, !is.gap(new[i,])] = a[i, !is.gap(a[i, ])]
  }
  new[dims[1], !is.gap(new[dims[1],])] = b[2, !is.gap(b[2, ])]

  new[is.na(new)]="-"
  gaps=gap.inspect(new)
  new = new[, which(gaps$col<dims[1])]
  
  if(dims[1]>3) 
    return(.tmalignGapRefine(new))
  else 
    return(new)
}


.tmalignGapRefine <- function(aln) {

  tmpfasta <- tempfile()  
  gaps <- gap.inspect(aln)
  dims <- dim(aln)

  if(length(gaps$t.inds)==0)
    return(aln)
  
  bs <- bounds(gaps$t.inds)
  bs=bs[which(bs[,"length"]>1),]

  hmm <- list()
  for ( i in 1:nrow(bs)) {
    tmpaln <- NULL
    tmpaln$ali=aln[, bs[i,"start"]:bs[i,"end"]]
    rownames(tmpaln$ali) <- paste("seq", 1:dims[1], sep="")
    tmpaln$id=rownames(tmpaln$ali)

    seq.dims <- apply(tmpaln$ali, 1, function(x) length(x[!is.gap(x)]))
    if(max(seq.dims)>1)
      tmp <- seqaln(tmpaln, outfile=paste("hmm", i, sep=""))$ali
    else {
      tmp = apply(tmpaln$ali, 1, function(x) if(all(is.gap(x))) return(NA) else return(x[!is.gap(x)]))
      tmp = matrix(tmp, ncol=1)
    }

    tmp[is.gap(tmp)]="-"
    rownames(tmp) <- paste("seq", 1:dims[1], sep="")
    
    do.ent <- FALSE
    if(do.ent) {
      entro <- entropy(tmp)$H.10.norm
      su <- sum(entro)/length(entro)
      if(su>0.2)
        hmm[[i]]=tmp
      else
        hmm[[i]]=tmpaln$ali
    }
    else {
      hmm[[i]]=tmp
    }
  }
  
  newaln <- aln
  for ( i in 1:nrow(bs)) {
    newaln[,bs[i, "start"]:bs[i, "end"]]=NA
    newaln[,bs[i, "start"]:(bs[i, "start"]+ncol(hmm[[i]])-1) ]=hmm[[i]]
  }

  newaln[is.na(newaln)]="-"
  gaps = gap.inspect(newaln)
  return(newaln[, which(gaps$col<dims[1])])

}
