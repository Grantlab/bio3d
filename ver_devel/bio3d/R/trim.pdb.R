"trim.pdb" <-
  function(pdb, inds=NULL, sse=TRUE) {
    if(class(pdb)!="pdb")
      stop("input 'pdb' must be a list object as returned from 'read.pdb'")

    if(is.null(inds))
      stop("no selection indices provided")

    if(!is.list(inds))
      stop("selection indices must be provided i.e. from 'atom.select'")

    if(is.null(inds$atom) || is.null(inds$xyz))
      stop("selection indices must be provided i.e. from 'atom.select'")

    ## Trim main components
    atom <- pdb$atom[inds$atom,]
    xyz <-  pdb$xyz[inds$xyz]
    calpha <- as.logical(atom[,"elety"]=="CA")

    ## Multi-model records
    if(!is.null(pdb$xyz.models))
      xyz.models <- pdb$xyz.models[,inds$xyz]
    else
      xyz.models <- NULL
    
    sse.unbound <- function(sse) {
      ex <- NULL; ch <- NULL; ty <- NULL;
      for(i in 1:length(sse$start)) {
        tmp <- sse$start[i]:sse$end[i]
        ex <- c(ex, tmp)
        ch <- c(ch, rep(sse$chain[i], length(tmp)))
        ty <- c(ty, rep(sse[[4]][i], length(tmp)))
      }
      ub <- cbind(ex, ch, ty) ## Matrix of unbound nums chains types
      ub <- cbind(ub, apply(ub, 1, function(x) paste(x[1:2], collapse="_")))
      colnames(ub) <- c("nums", "chains", "type", "strid")
      return(ub)
    }

    bounds.inds <- function(nums) {
      if (length(nums) == 1) {
        return(1)
      }
      
      bounds.inds <- 1
      nums.start <- nums[1]
      diff.i <- 1
      for (i in 2:length(nums)) {
        if ((nums[i] - diff.i) != nums.start) {
          bounds.inds <- c(bounds.inds, i)
          nums.start <- nums[i]
          diff.i <- 1
        }
        else {
          diff.i <- diff.i + 1
        }
      }
      return(bounds.inds)
    }

    helix <- NULL; sheet <- NULL;

    if(sse) {
      ## What's left after trimming
      trimmed <- atom[calpha,c("resno", "chain")]

      ## Build a string identifier for easy comparison
      strid <- apply(trimmed, 1, function(x) paste(x, collapse="_"))
      trimmed <- cbind(trimmed, strid)
      colnames(trimmed) <- c("nums", "chains", "strid")
      
      ## Helices
      if(length(pdb$helix$start)>0) {
        tmp <- sse.unbound(pdb$helix)
        ht <- tmp[tmp[,"strid"] %in% trimmed[,"strid"], 1:3]
        new.helix <- matrix(ht, ncol=3, byrow=FALSE)
        colnames(new.helix) <- c("nums", "chains", "type")
        
        if(nrow(new.helix)>0) {
          nums <- as.numeric(new.helix[,"nums"])
          bnds <- bounds(nums, dup.inds=FALSE, pre.sort=FALSE)
          
          tmp.inds <-  bounds.inds(nums)
          chain <- new.helix[tmp.inds,"chains"]
          type <- new.helix[tmp.inds,"type"]
          bnds <- cbind(bnds, chain, type)
          
          helix$start = as.numeric(bnds[,"start"])
          helix$end = as.numeric(bnds[,"end"])
          helix$chain = as.vector(bnds[,"chain"])
          helix$type = as.vector(bnds[,"type"])
        }
      }
      
      ## Sheets
      if(length(pdb$sheet$start)>0) {
        tmp <- sse.unbound(pdb$sheet)
        ht <- tmp[tmp[,"strid"] %in% trimmed[,"strid"], 1:3]
        new.sheet <- matrix(ht, ncol=3, byrow=FALSE)
        colnames(new.sheet) <- c("nums", "chains", "type")
        
        if(nrow(new.sheet)>0) {
          nums <- as.numeric(new.sheet[,"nums"])
          bnds <- bounds(nums, dup.inds=FALSE, pre.sort=FALSE)
          
          tmp.inds <-  bounds.inds(nums)
          chain <- new.sheet[tmp.inds,"chains"]
          type <- new.sheet[tmp.inds,"type"]
          bnds <- cbind(bnds, chain, type)
          
          sheet$start = as.numeric(bnds[,"start"])
          sheet$end = as.numeric(bnds[,"end"])
          sheet$chain = as.vector(bnds[,"chain"])
          sheet$sense = as.vector(bnds[,"type"])
        }
      }
    }

    output<-list(atom=atom,
                 het=pdb$het, ## return unmodified
                 helix=helix, 
                 sheet=sheet, 
                 seqres=pdb$seqres, ## return unmodified
                 xyz=xyz,
                 xyz.models=xyz.models,
                 calpha = calpha)
    
    class(output) <- "pdb"
    return(output)
  }
