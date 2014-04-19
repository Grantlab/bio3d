print.fasta <- function(x, width=NULL, col.inds=NULL, numbers=TRUE, ...) {

  ##-- Print sequence alignment in a nice formated way
  ##    source("print.aln.R")
  ##    x<-read.fasta("poo.fa")
  ##    print.aln(x)
  ##
  ##    file <- system.file("examples/kif1a.fa",package="bio3d")
  ##    aln  <- read.fasta(file)
  ##    print.aln(aln, width=40)
  ##    print.aln(aln, width=60)
  ##    print.aln(aln, col=c(10,12,14,22,90:100), numbers=F)
  ## ToDo:
  ##   Does not work if alignment contains only one position (one seq?)
  ##     y=x; y$ali=x$ali[,1]


  if(class(x) %in% c("fasta", "3dalign")) {
    id <- x$id
    ali <- as.matrix(x$ali)
  } else {
    ali <- as.matrix(x)
    id <- rownames(ali)
  }

  ##- Trim to 'col.inds' if provided 
  if(!is.null(col.inds)) {
    ali <- ali[,col.inds, drop=FALSE]
  }

  ##- Format sequence identifiers
  ids.nchar <- max(nchar(id))+3 ## with a gap of 3 spaces btwn id and sequence
  ids.format <- paste0("%-",ids.nchar,"s")
  ids <- sprintf(ids.format, id)

  ## Format for annotation printing (see below)
  pad.format <- paste0("%+",(ids.nchar+1),"s")

  ##- Scale 'width' of output if not specified in input call
  tput.col <- 85  ## typical terminal width from system("tput cols")
  if(is.null(width)) {
    width <- tput.col - ids.nchar - 4
   }

  ## Make sure we end on a 10 block
  width <- floor(width/10)*10

  ##- Work out sequence block widths
  nseq <- length(ids)
  nres <- ncol(ali)

  block.start <- seq(1, nres, by=width)
  if(nres < width) {
    block.end <- nres
  } else {
   block.end <- unique(c( seq(width, nres, by=width), nres))
  }
  nblocks <- length(block.start)

  block.annot  <- rep(" ", width)
  block.annot[ c(1,seq(10, width, by=10)) ] = "."

  blocks <- matrix(NA, ncol=nblocks, nrow=nseq) 
  for(i in 1:nblocks) {
    ##- Sequence block
    positions <- block.start[i]:block.end[i]
    blocks[,i] <- paste0(ids, apply(ali[, positions, drop=FALSE], 1, paste, collapse=""))

    ##-- Formated Printing of annotations (numbers & ticks) and sequence blocks
    if(numbers) {
      ##- Annotations for each sequence block
      annot = block.annot[1:length(positions)]
      annot[length(annot)] = block.end[i]
      annot[1] = sprintf(pad.format, block.start[i])
      cat(paste(annot, collapse=""),"\n")
    }

    ##- Sequence block
    cat(blocks[,i], sep="\n")

    ##- Ticks + numbers again
    if(numbers) {
      cat(paste(annot, collapse=""),"\n\n")
    } else{ cat("\n") }
  }

  ## Attribute summary
  j <- paste(attributes(x)$names, collapse = ", ")
  cat(strwrap(paste(" + attr:", j, "( dim(x$ali) =", 
    paste(dim(ali),collapse="x"),")\n"), 
    width = 70, exdent = 8), sep = "\n")

  ##invisible(blocks) ## Can be useful for plot.fasta() later!!
}

