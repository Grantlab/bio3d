membership.comparison <- function(x, important.members=important.members, plot=FALSE, comm.color=1, sse=NULL, dash.lines=FALSE){

  ## Check if x has 'cna' class
  if(!is.numeric(x)){
    stop("Input 'x' object must be numeric")
  }

  if(is.vector(x)){
    y <- matrix(x, nrow=1, ncol=length(x))
    colnames(y) <- names(x)
    x <- y
  }
  
  if(is.null(rownames(x))){
    warning("No row names on the membership matrix! Please add them\n")
    rownames(x) <- c(1:dim(x)[1])
  }

  ## Initializing vectors for the plot
  memb.matrix <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[2])
  
  for(i in 1:dim(x)[1]){
    ## Determine the biggest community in the important.members group
    size <- table(x[i,important.members])

    biggest <- as.numeric(names(which(size==max(size))))

    if(length(which(size==max(size))) > 1){
      warning(paste0("Communities ", paste0(names(which(size == max(size))), collapse=" and "), " have the same number of elements in the specified residue subset."))
      print(paste0("Only community ", names(which(size == max(size))[1])," will be plotted in method ", rownames(x)[i]))
      biggest <- as.numeric(names(which(size==max(size))))[1]
    }
      
    big.comm.inds <- as.numeric(names(which(x[i,]==biggest)))
    other.comm.inds <- as.numeric(names(which(x[i,]!=biggest)))
 
    memb.matrix[i,big.comm.inds] <- 1  
   
    memb.matrix[i,other.comm.inds] <- 0
  }

  rownames(memb.matrix) <- rownames(x)
    
  ## Let's make the plot...
  if(plot){
    
    if(is.numeric(comm.color)){
      color <- vmd.colors(comm.color)[comm.color]
    }
    if(is.character(comm.color)){
      color <- comm.color
    }
  
    par(mar=c(5,9,1,0.5))
    image(1:ncol(memb.matrix), 1:nrow(memb.matrix), t(memb.matrix), col=c("white",color), ylab="", axes=FALSE, xlab="Residue No", cex.lab=1.2)
    axis(BELOW<-1, at=seq(1,ncol(memb.matrix),25), labels=seq(1,ncol(memb.matrix),25), cex.axis=0.8)
    axis(LEFT <-2, at=seq(1,nrow(memb.matrix),1),labels=rownames(memb.matrix),
         las= HORIZONTAL<-1,
         cex.axis=0.8,
         )
 
    if(!is.null(sse)){
      sse.draw(sse, ylim=c(0,dim(x)[1]))
    }
  
    if(dash.lines){
      for(b in 1:(nrow(memb.matrix)-1)){
        abline(h=(b+0.5), lty=2)
      }
    }
  }
  
  return(memb.matrix)
}
 
