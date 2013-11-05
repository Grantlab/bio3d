summary.cna <- function(object, ...) {

  ## summary.cna(net)
  ## y <- summary.cna(net, file="tmp.tbl", col.names=FALSE, append=T)
  ## or just write 'y' to file!
  
  if( !"cna" %in% class(object) ) {
    stop("Input should be a cna network object")
  }

  size <- table(object$raw.communities$membership)
  id   <- names(size)
  memb <- sapply(id, function(i) { which(object$raw.communities$membership==i) } )
  ## NOTE: Perhaps the memb should be names() of which inds
  ##       rather than the inds themselves as it is curently?
  ##memb <- sapply(id, function(i) { names(which(object$membership==i)) } )

  ##- Format as condensed vector for printing
  if( is.numeric(unlist(memb)) ) {
    members <- rep(NA, length(id))
    for(i in 1:length(id)) {
      b <- bounds(memb[[i]])[,c("start","end"),drop=FALSE]
      members[i] <- paste0("c(",
         paste( apply(b, 1, paste, collapse=":"), collapse=", "), ")")
    }
  } else {
    ##- non numeric vectors can not be condensed
    members <- unlist(lapply(memb, paste, collapse=", "))
  }
  
  ## Print to terminal
  ###write.table(cbind(id, size, members), sep="\t",
  ###            row.names=FALSE, quote=FALSE)
  ## NOTE: Maybe strwrap(), sprintf() or format() could help 
  ##       indent the second line of the above print when members
  ##       is over a certain length?
  ## ALSO: Might be useful to sumarize other attributes in object

  ## Print with groupings (clustered membership) and ordered by member size
  z <- cbind(Group=object$clustered.communities$members, id, size, members)
  a <- z[order(z[,"Group"], -as.numeric(z[,"size"]), as.numeric(z[,"id"])),]

  b <- bounds(a[,"Group"], dup.inds=T)[,c("start","end"), drop=FALSE]
  a[ duplicated(a[,"Group"]),"Group"] = " "
  
  for(i in 1:nrow(b)) {
    inds <- unbound(b[i,1], b[i,2] )
    write.table(a[inds,,drop=FALSE], na=" ", sep="\t",
                row.names=FALSE, quote=FALSE, ...)

    ## Summary of groups (i.e. clustered membership)
    cat(paste("   -- Sum: ",nrow(a[inds,,drop=FALSE]),
              "(communities) |", sum(as.numeric(a[inds,"size",drop=FALSE])),
              "(members) --",sep=" "),"\n\n")
  }

  ## Output silently as a list
  tbl <- data.frame( apply(z[,1:3],2,as.numeric), members=z[,4],
                    stringsAsFactors=FALSE )
  y <- list("id"=id, "size"=size, "members"=memb, "tbl"=tbl)
}

