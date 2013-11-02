dssp.block <- function(cij, sse.membership, collapse.method="max"){

  if(dim(cij)[1] != dim(cij)[2]){
    stop("The cij matrix is not squared")
  }
  
  ## Definition of function to collapse the cij matrix
  contract.matrix <- function(cij.abs, membership,## membership=comms$membership,
                              collapse.method="max"){
    
    ## Function to collapse a NxN matrix to an mxm matrix
    ##  where m is the communities of N. The collapse method
    ##  can be one of the 'collapse.options' below
    collapse.options=c("max", "median", "mean", "trimmed")
    collapse.method <- match.arg(tolower(collapse.method), collapse.options)

    ## Fill a 'collapse.cij' community by community matrix
    node.num <- max(membership)
    collapse.cij <- matrix(0, nrow=node.num, ncol=node.num)
    inds <- pairwise(node.num)

    for(i in 1:nrow(inds)) {
      comms.1.inds <- which(membership==inds[i,1])
      comms.2.inds <- which(membership==inds[i,2])
      submatrix <- cij.abs[comms.1.inds, comms.2.inds]

      ## Use specified "collapse.method" to define community couplings
      collapse.cij[ inds[i,1], inds[i,2] ] = switch(collapse.method,
                  max = max(submatrix),
                  median = median(submatrix),
                  mean = mean(submatrix),
                  trimmed = mean(submatrix, trim = 0.1))
    }

    ## Copy values to lower triangle of matrix and set class and colnames
    collapse.cij[ inds[,c(2,1)] ] = collapse.cij[ inds ]
    colnames(collapse.cij) <- 1:ncol(collapse.cij)
    class(collapse.cij) <- c("dccm", "matrix")
    return(collapse.cij)
  }

  if(class(sse.membership) == "dssp"){
    ## Helix initialization
    start <- sort(sse$helix$start)
    end <- sort(sse$helix$end)

    ## Cutoff on helix length
    helix.inds <- (end - start) > 4
  
    sse$helix$start <- start[helix.inds]
    sse$helix$end <- end[helix.inds]

    names(sse$helix$start) <- paste0("H",c(1:length(sse$helix$start)))
    names(sse$helix$end) <- paste0("H",c(1:length(sse$helix$end)))
  
    ## Sheet inizialization
    start <- sort(sse$sheet$start)
    end <- sort(sse$sheet$end)

    ## Cutoff on helix length
    sheet.inds <- (end - start) > 4
  
    sse$sheet$start <- start[sheet.inds]
    sse$sheet$end <- end[sheet.inds]

    names(sse$sheet$start) <- paste0("S",c(1:length(sse$sheet$start)))
    names(sse$sheet$end) <- paste0("S",c(1:length(sse$sheet$end)))
  
    start <- sort(c(sse$helix$start, sse$sheet$start))
    end <- sort(c(sse$helix$end, sse$sheet$end))

    loop.num <- 1

    ## Add first loop, in case the ss does not start with an helix of a strand
    if(start[1] != 1){
      end <- c(start[1]-1,end)
      start <- c(1, start)
      
      names(start)[1] <- paste0("L", loop.num)
      names(end)[1] <- paste0("L", loop.num)
      
      loop.num <- loop.num + 1
      
    }

    sse.elements <- length(start)

    ## Add the loop's extremes to the start and end vectors
    for(i in 1:sse.elements){
   
      if((i == sse.elements) & (end[length(end)] != length(sse$sse))){
        start <- sort(start)
        end <- sort(end)
      
        start <- c(start,(end[length(end)]+1))
        end <- c(end,length(sse$sse))
      
        names(start)[length(start)] <- paste0("L",loop.num)
        names(end)[length(end)] <- paste0("L",loop.num)

      }
      else if((i < sse.elements) & end[i] != start[i+1]-1){
        start <- c(start, (end[i]+1))
        end <- c(end, (start[i+1]-1))
      
        names(start)[length(start)] <- paste0("L",loop.num)
        names(end)[length(end)] <- paste0("L",loop.num)

        loop.num <- loop.num + 1
      }
    }
 
    ## Order the start and end vectors
    start <- sort(start)
    end <- sort(end)

    ## Create the sse membership vector
    sse.block <- rep(0,length(sse$sse))
    for(i in 1:length(start)){
      sse.block[start[i]:end[i]] <- i
      names(sse.block)[start[i]:end[i]]  <- names(start)[i]
    }
  }
  else if(is.numeric(sse.membership)){
    
    if(length(sse.membership) != dim(cij)[1]){
      stop("The length of the 'sse.membership' is different from the cij matrix dimensions")
    }
    
    sse.block <- sse.membership
  }

  ## Collapsing cij...
  cij.collapsed <- contract.matrix(cij, sse.block, collapse.method=collapse.method)

  rownames(cij.collapsed) <- unique(names(sse.block))
  colnames(cij.collapsed) <- unique(names(sse.block))
  
  output <- list("sse.membership"=sse.block, "cij.collapsed"=cij.collapsed)
  
  return(output)
}
