atom2ele <- function(...)
  UseMethod("atom2ele")

atom2ele.default <- function(x, elety.custom=NULL, rescue=TRUE, ...){
  if(!is.null(elety.custom)) {
    if(!all(c("name","symb") %in% names(elety.custom)))
      stop("'elety.custom' must contains 'name' and 'symb' components")
    inds <- unlist(lapply(elety.custom, is.factor))
    elety.custom[inds] <- lapply(elety.custom[inds], as.character)
  }
  atom.index <- rbind(elety.custom[,c("name","symb")], bio3d::atom.index[,c("name","symb")])
  # Why atom names starting by "H" are directly converted to "H" as follow?
  # x[substr(x,1,1) == "H"] <- "H"
  symb <- atom.index[match(x, atom.index[,"name"]), "symb"]
  is.unknown <- is.na(symb)
  if(any(is.unknown)) {
    if(rescue) {
        ## vector of unknown elements
        unknowns <- unique(x[is.unknown])
        symb2 <- rep(NA, length(unknowns))
        names(symb2) <- unknowns

        ## format element names before matching
        spl <- strsplit(unknowns, "")
        totest <- lapply(spl, function(b) {
            ## remove numbering from atom name (e.g. FE2, C4A)
            inds <- grep("[^0-9]", b)
            if(length(inds) == 0) return(b)
            j <- bounds(inds)[1, c("start", "end")]
            b <- b[j[1]:j[2]]

            ## First char to upper, remaining to lower case
            new <- NULL
            for(i in 1:length(b)){
                new <- c(new, ifelse(i==1, toupper(b[i]), tolower(b[i])))
            }
            return(paste(new, collapse=""))
        })

        ## match with bio3d::elements$symb
        totest <- unlist(totest)
        if(any(totest %in% bio3d::elements$symb)) {
            symb2[unknowns[totest %in% bio3d::elements$symb]] <-
                totest[totest %in% bio3d::elements$symb]
        }

        ## try with first character to see if it matches elements$symb
        if(any(is.na(symb2))) {
            na.inds <- which(is.na(symb2))
            totest <- toupper(substr(names(symb2[na.inds]), 1, 1))

            if(any(totest %in% bio3d::elements$symb)) {
                rplc <- names(symb2[na.inds])[totest %in% bio3d::elements$symb]
                symb2[rplc] <- totest[totest %in% bio3d::elements$symb]
            }
        }
        
        ## stop with error in case of un-mapped elements
        if(any(is.na(symb2))) {
            stop("\telements could not be determined for: ",
                 paste(names(symb2)[is.na(symb2)], collapse=", "))
        }

        ## inform user on mapped elements
        warning(paste("\n\tmapped element ", names(symb2), " to ", symb2, sep=""))
        
        ## include matched elements to original symbol vector
        symb[is.unknown] <- symb2[ x[is.unknown] ]
    }
    else {
      stop("\telements could not be determined for: ",
           paste(unique(x[is.unknown]), collapse=", "))
    }
  }
  symb <- unlist(symb)
  return(symb)
}

atom2ele.pdb <- function(pdb, inds=NULL, ...){
  if(!is.null(inds))
    pdb <- trim(pdb, inds)
  atom.names <- pdb$atom[,"elety"]
  return(atom2ele.default(atom.names, ...))
}
