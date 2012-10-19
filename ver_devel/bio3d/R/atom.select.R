"atom.select" <-
function(pdb, string=NULL,
         chain=NULL, resno=NULL, resid=NULL,
         eleno=NULL, elety=NULL,
         verbose=TRUE, rm.insert=FALSE) {

  
  ## Version 0.3 ... Thu May 12 17:38:25 PDT 2011
  ##
  ##  Removed old "pdb.summary()" section and replaced
  ##   with call to new "print.pdb()" function.
  ##
  ## Version 0.2 ... Tue Jun 23 18:15:28 PDT 2009
  ##
  ##  Builds selection string from components 
  ##   chain, resno, resid, eleno, elety
  ##
  # Version 0.1 ... Tue Mar 21 10:58:43 PST 2006
  #
  #   Prints a summary of 'pdb' makeup if called 
  #    without a selection 'string'.
  #   Also added 'string' shortcuts "calpha", 
  #    "back", "backbone" and "cbeta"
  #
  # Version 0.0 ... Fri Mar 17 14:45:37 PST 2006
  #
  #  Description:
  #   Atom selection function
  #    Losely based on PyMol selection Macro:
  #   see: http://www.pymolwiki.org/index.php/Selection_Macros
  #
  #   String Selection Syntax:
  #    "//A/130:142///N,CA,C,O/"
  #    "/segid/chain/resno/resid/eleno/elety/"
  #
  #  E.g.
  #   # read a PDB file
  #   pdb<-read.pdb("1bg2.pdb")
  #   # print a structure summary
  #   atom.select(pdb)
  #   # select all C-alpha atoms from resno 65 to 143
  #   ca.inds   <- atom.select(pdb, "///65:143///CA/")
  #   # or all C-alphas
  #   ca.inds   <- atom.select(pdb, "calpha")
  #   # more examples
  #   inds<-atom.select(pdb, "//A/130:142///N,CA,C,O/")


  if (missing(pdb)) {
    stop("atom.select: must supply 'pdb' object, e.g. from 'read.pdb'")
  }

  
  sel.txt2nums <- function (num.sel.txt) {
    
    # Splitting functions for numbers

    num1<-unlist(strsplit(num.sel.txt, split=","))
    num2<-suppressWarnings( as.numeric(num1) )
    # comma split still may have "10:100" = NA in num2
    tosplit <- num1[ which(is.na(num2)) ]
    num3 <- unlist(strsplit(tosplit, split=":"))
    # pair-up num3 to make num4
    num4<-NULL; i<-1
    while (i < length(num3) ) {
      num4 <- c(num4, as.numeric(num3[i]):
                as.numeric(num3[i+1]) )
      i<-i+2
    }
    # join and order num2 with num4 
    return( sort(unique(c(na.omit(num2),num4))) )
  }

  sel.txt2type <- function (type.sel.txt) {
    
    # Splitting functions for characters
    
    type1 <- unlist(strsplit(type.sel.txt, split=","))
    # split on coma and remove white space
    return( gsub(" ","",type1) )
  }

  ##- No selection string, just print PDB summary
  if(is.null(string)) {
    if(is.null(c(chain,resno,resid,eleno,elety))) {
      print.pdb(pdb)  ## -- Edit May 12 2011
      return()
    }
    ## Build selection string from components
    ## /segid/chain/resno/resid/eleno/elety/
    string <- paste("//",paste(chain, collapse=","),"/",
                    paste(resno, collapse=","),"/",
                    paste(resid, collapse=","),"/",
                    paste(eleno, collapse=","),"/",
                    paste(elety, collapse=","),"/",sep="")
    
    
    rm(chain,resno,resid,eleno,elety)
  } ##else {

  # string shortcuts
    if (string=="calpha" || string=="CA") {
      string= "//////CA/"
    }
    if (string=="cbeta" || string=="CB") {
      string= "//////N,CA,C,O,CB/"
    }
    if (string=="backbone" || string=="back") {
      string= "//////N,CA,C,O/"
    }
    if (string=="all") {
      string= "///////"
    }
    ## - Edit Jan 17 2008
    if (string=="h") {
      h.atom <- which( substr(pdb$atom[,"elety"], 1, 1) %in% "H" )
      match <- list(atom=h.atom, xyz=atom2xyz(h.atom))
      class(match) <- "select"
      if(verbose) 
        cat(paste(" *  Selected", length(h.atom), "hydrogen atoms *\n"))
      return(match)
    }
    if (string=="noh") {
      noh.atom <- which( !substr(pdb$atom[,"elety"], 1, 1) %in% "H" )
      match <- list(atom=noh.atom, xyz=atom2xyz(noh.atom))
      class(match) <- "select"
      if(verbose) 
        cat(paste(" *  Selected", length(noh.atom), "non-hydrogen atoms *\n"))
      return(match)
    }
    
    
    ## - end edit
    
    # main function  
    sel <- unlist(strsplit(string, split = "/"))

    if (sel[1] == "") { # full selection string (starts with "/")
      sel <- sel[-1]
      if(length(sel) != 6) {
        print("missing elements, should be:")
        print("/segid/chain/resno/resid/eleno/elety/")
      }

      names(sel) <- c("segid","chain", "resno","resid","eleno","elety")
      #print(sel)

      blank <- rep(TRUE, nrow(pdb$atom) )
      sel.inds <- NULL

    # SEGID
      if(sel["segid"] != "") { 
        sel.inds <- cbind(sel.inds,
                          segid=is.element( pdb$atom[,"segid"],
                            sel.txt2type( sel["segid"] )) )
      } else {  sel.inds <- cbind(sel.inds, segid=blank)  }

    # CHAIN
      if(sel["chain"] != "") {
        sel.inds <- cbind(sel.inds,
                          chain=is.element( pdb$atom[,"chain"],
                            sel.txt2type( sel["chain"] )) )        
      } else { sel.inds <- cbind(sel.inds, chain=blank)  }
  
    # RESNO
      if(sel["resno"] != "") {
        
        rn <- sel.txt2nums( sel["resno"] )
        if(is.numeric(rn) & length(rn)==0) {
          # check for R object 
          rn <- get(gsub(" ","",sel["resno"]))
          
        }
        
        sel.inds <- cbind(sel.inds,
                          resno=is.element( as.numeric(pdb$atom[,"resno"]),
                            rn))
                            #sel.txt2nums( sel["resno"] )) )
      } else {  sel.inds <- cbind(sel.inds, resno=blank)  }

    # RESID
      if(sel["resid"] != "") {
        sel.inds <- cbind(sel.inds,
                          resid=is.element(pdb$atom[,"resid"],
                            sel.txt2type( sel["resid"] )) )
      } else {  sel.inds <- cbind(sel.inds, resid=blank)  }

    # ELENO
      if(sel["eleno"] != "") {
        sel.inds <- cbind(sel.inds,
                          eleno=is.element(as.numeric(pdb$atom[,"eleno"]),
                            sel.txt2nums( sel["eleno"] )) )
      } else {  sel.inds <- cbind(sel.inds, eleno=blank)  }

    # ELETY
      if(sel["elety"] != "") {
      ##  cat( sel["elety"] ,"\n" ) ### glob2rx
      #if(any(i <- grep("*", sel["elety"]))) {
      #  print("WARN: no wild card '*' matching, yet")
      #}
        
        sel.inds <- cbind(sel.inds,
                          elety=is.element(pdb$atom[,"elety"],
                            sel.txt2type( sel["elety"] )) )
      } else {  sel.inds <- cbind(sel.inds, elety=blank)  }

      match.inds <- ( (apply(sel.inds, 1, sum, na.rm=TRUE)==6) )

      if (rm.insert) { # ignore INSERT records
        insert <- which(!is.na(pdb$atom[,"insert"]))
        match.inds[insert] <- FALSE
      }
      # return XYZ indices
      xyz.inds <- matrix(1:length( pdb$atom[,c("x","y","z")] ),nrow=3,byrow=FALSE)
      xyz.inds <- as.vector(xyz.inds[,match.inds])

      if (verbose) {
        sel <- rbind( sel, apply(sel.inds, 2, sum, na.rm=TRUE) )
        rownames(sel)=c("Stest","Natom"); print(sel)
        cat(paste(" *  Selected a total of:",sum(match.inds),
                  "intersecting atoms  *"),sep="\n")
      }
      
      match <- list(atom=which(match.inds), xyz=xyz.inds)
      class(match) <- "select"
      return(match)
    }

}

