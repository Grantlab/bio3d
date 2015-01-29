".is.protein" <- function(pdb, query.names = c("C", "CA", "N", "O")) {
  return(.residue.type(pdb, query.names))
}

".is.nucleic" <- function(pdb, query.names = c("O3'", "C3'", "C4'", "C5'", "O5'")) {
  return(.residue.type(pdb, query.names))
}

".residue.type" <- function(pdb, query.names=NULL) {
  
  ## nucleic can be both XX' and XX*
  pdb$atom$elety <- gsub("[*]", "'", pdb$atom$elety)
  pdb$atom$insert[is.na(pdb$atom$insert)] <- ""
  pdb$atom$chain[is.na(pdb$atom$chain)]   <- ""
  
  ## make unique residue identifiers
  resid <- paste(pdb$atom$chain, pdb$atom$insert, pdb$atom$resno, sep="-")
  resid.unq <- unique(resid)
  
  ## selection vector
  sele <- rep(FALSE, length(resid))
  
  for( res in resid.unq ) {
    inds <- which(resid==res)
    sele[inds] <- all(query.names %in% pdb$atom$elety[inds])
  }
  
  return(sele)
}


"atom.select.pdb" <-
function(pdb, string=NULL,
         chain=NULL, resno=NULL, resid=NULL,
         eleno=NULL, elety=NULL, type=NULL,
         verbose=TRUE, rm.insert=FALSE, inverse=FALSE, ...) {

  ## Version 0.0 ... Fri Mar 17 14:45:37 PST 2006
  ##
  ##  Description:
  ##   Atom selection function
  ##    String Selection Syntax:
  ##    "//A/130:142///N,CA,C,O/"
  ##    "/segid/chain/resno/resid/eleno/elety/"
  ##
  ##  E.g.
  ##   # read a PDB file
  ##   pdb<-read.pdb("1bg2")
  ##   # print a structure summary
  ##   print.pdb(pdb)
  ##   # select all C-alpha atoms from resno 65 to 70
  ##   ca.inds   <- atom.select(pdb, "///65:70///CA/")
  ##
  ##   # or use C-alpha string shortcut
  ##   ca.inds <- atom.select(pdb, "calpha")
  ##   ca.inds <- atom.select(pdb, "calpha", resno=65:70)
  ##
  ##   # more examples
  ##   inds<-atom.select(pdb, "//A/130:142///N,CA,C,O/")
  ##   inds<-atom.select(pdb, chain="A", resno=130:142, elety="N,CA,C,O")
  ##
  ##   # or using string shortcut
  ##   inds <- atom.select(pdb, "back", resno=65:70)


  cl <- match.call()

  sel.txt2nums <- function(num.sel.txt) {
    ##- Splitting function for numbers  - split on coma & colon
    num1 <- unlist( strsplit(num.sel.txt, split=",") )
    num2 <- suppressWarnings( as.numeric(num1) )

    if( any(is.na(num2)) ) {
      ## Split range elements (e.g. may have "10:100" = NA in num2)
      num3 <- unlist(strsplit( num1[ which(is.na(num2)) ], split=":"))
      ## Pair-up range elements in num3 to make num4
      num3 <- matrix(as.numeric(num3),nrow=2)
      num4 <- unbound(num3[1,], num3[2,])
      num2 <- c(num2, num4)
    }
    return( sort(unique(c(na.omit(num2)))) )
  }

  sel.txt2type <- function(type.sel.txt, na=FALSE) {
    ##- Splitting function for characters
    if(is.na(type.sel.txt) || (type.sel.txt==" ")) {
      ## Blank or NA selections will return NA
      sel <- NA
    } else {
      ## Split larger strings on coma & remove white space
      sel <- gsub(" ","", unlist(strsplit(type.sel.txt, split=",")) )
      if(na) { sel[sel=="NA"]=NA }
    }
    return(sel)
  }

  ##-- Parse string and return the selection
  parse.string <- function(pdb, string, verbose, rm.insert, type="") {

     ##-- We have input selection string
     aa <- unique(pdb$atom[,"resid"])
     hoh <- c("H2O", "OH2", "HOH", "HHO", "OHH", "SOL", "WAT",
              "TIP", "TIP2", "TIP3", "TIP4")
     
     protein.atoms <- NULL; not.protein.atoms <- NULL;
     if(string %in% c("protein", "calpha", "cbeta", "backbone", "back", "notprotein", "ligand")) {
       protein.atoms     <- .is.protein(pdb)
       not.protein.atoms <- !protein.atoms
     }
     
     nucleic.atoms <- NULL; not.nucleic.atoms <- NULL;
     if(string %in% c("nucleic", "notnucleic", "ligand")) {
       nucleic.atoms     <- .is.nucleic(pdb)
       not.nucleic.atoms <- !nucleic.atoms
     }

     water.atoms <- NULL
     if(string %in% c("water", "notwater", "ligand"))
       water.atoms       <- pdb$atom$resid %in% hoh
     
     ligand.atoms <- NULL; water.atoms <- NULL;
     if(string %in% c("ligand"))
       ligand.atoms      <- !protein.atoms & !nucleic.atoms & !water.atoms
     
     
     if (string=="h") {
       h.atom <- which( substr(pdb$atom[,"elety"], 1, 1) %in% "H" )
       match <- list(atom=h.atom, xyz=atom2xyz(h.atom))
       class(match) <- "select"
       if(verbose)
         cat(paste(" *  Selected", length(h.atom), "hydrogen atoms *\n"))
       return(match)
     }
     if (string=="noh") {
      ## gsub fix added to catch non-standard 1HH1, 2HH1, 3HH1 etc... (Fri Jan 24 EST 2014)
       noh.atom <- which( !substr( gsub("^[123]", "",pdb$atom[,"elety"]) , 1, 1) %in% "H" )
       match <- list(atom=noh.atom, xyz=atom2xyz(noh.atom))
       class(match) <- "select"
       if(verbose)
         cat(paste(" *  Selected", length(noh.atom), "non-hydrogen atoms *\n"))
       return(match)
     }

     ##-- Check for string 'shortcuts'
     i <- switch(string,
                 all = "///////",
                 calpha   = paste("/////", paste(pdb$atom$eleno[protein.atoms], collapse=","), "/CA/",sep=""),
                 cbeta    = paste("/////", paste(pdb$atom$eleno[protein.atoms], collapse=","), "/CB/",sep=""),
                 backbone = paste("/////", paste(pdb$atom$eleno[protein.atoms], collapse=","), "/N,CA,C,O/",sep=""),
                 back     = paste("/////", paste(pdb$atom$eleno[protein.atoms], collapse=","), "/N,CA,C,O/",sep=""),

                 protein    = paste("/////", paste(pdb$atom$eleno[protein.atoms], collapse=","), "//",sep=""),
                 notprotein = paste("/////", paste(pdb$atom$eleno[!protein.atoms], collapse=","), "//",sep=""),

                 nucleic    = paste("/////", paste(pdb$atom$eleno[nucleic.atoms], collapse=","), "//",sep=""),
                 notnucleic = paste("/////", paste(pdb$atom$eleno[!nucleic.atoms], collapse=","), "//",sep=""),

                 ligand   = paste("/////", paste(pdb$atom$eleno[ligand.atoms], collapse=","), "//",sep=""),
                 water    = paste("/////", paste(pdb$atom$eleno[water.atoms], collapse=","), "//",sep=""),
                 notwater = paste("/////", paste(pdb$atom$eleno[!water.atoms], collapse=","), "//",sep=""),

                 NA)

     if(is.na(i)) {
       ##- No string shortcut match

       if(!substr(string,1,1)=="/") {
         ## Check if we have a valid selection sting
         stop("Not a valid selection string shortcut.\n\t Please use one of:
            'calpha' 'cbeta' 'backbone'
            'protein' 'notprotein' 'ligand'
            'nucleic' 'notnucleic'
            'water' 'notwater'
            'h' 'noh'\n
       Or valid selection string:
            /segid/chain/resno/resid/eleno/elety/ \n")
       }

     } else {
       ##- Use string shortcut from switch function
       ###if(verbose)
       ###  cat(paste("\n Using selection 'string' keyword shortcut:",string, "=", i, "\n\n"))
       string = i
     }

     sel <- unlist(strsplit(string, split = "/"))

     if (sel[1] == "") { # full selection string (starts with "/")
       sel <- sel[-1]
       if(length(sel) != 6) {
         print("missing elements, should be:\n/segid/chain/resno/resid/eleno/elety/")
       }
       ## Add type ATOM/HETATM
       sel <- c(sel, paste(type,collapse=","))
       names(sel) <- c("segid","chain","resno","resid","eleno","elety","type")
       ##print(sel)

       blank <- rep(TRUE, nrow(pdb$atom) )
       sel.inds <- NULL

       ## SEGID
       if(sel["segid"] != "") {
         sel.inds <- cbind(sel.inds,
                           segid=is.element( pdb$atom[,"segid"],
                             sel.txt2type( sel["segid"] )) )
       } else {  sel.inds <- cbind(sel.inds, segid=blank)  }

       ## CHAIN
       if(sel["chain"] != "") {
         sel.inds <- cbind(sel.inds,
                           chain=is.element( pdb$atom[,"chain"],
                             sel.txt2type( sel["chain"], na=TRUE )) )
       } else { sel.inds <- cbind(sel.inds, chain=blank)  }

       ## RESNO
       if(sel["resno"] != "") {
         rn <- sel.txt2nums( sel["resno"] )
         if(is.numeric(rn) & length(rn)==0) {
           ## check for R object
           rn <- get(gsub(" ","",sel["resno"]))
         }
         sel.inds <- cbind(sel.inds,
                           resno=is.element( as.numeric(pdb$atom[,"resno"]),
                             rn))
       } else {  sel.inds <- cbind(sel.inds, resno=blank)  }

       ## RESID
       if(sel["resid"] != "") {
         sel.inds <- cbind(sel.inds,
                           resid=is.element(pdb$atom[,"resid"],
                             sel.txt2type( sel["resid"] )) )
       } else {  sel.inds <- cbind(sel.inds, resid=blank)  }

       ## ELENO
       if(sel["eleno"] != "") {
         sel.inds <- cbind(sel.inds,
                           eleno=is.element(as.numeric(pdb$atom[,"eleno"]),
                             sel.txt2nums( sel["eleno"] )) )
       } else {  sel.inds <- cbind(sel.inds, eleno=blank)  }

       ## ELETY
       if(sel["elety"] != "") {
         sel.inds <- cbind(sel.inds,
                           elety=is.element(pdb$atom[,"elety"],
                             sel.txt2type( sel["elety"] )) )
       } else {  sel.inds <- cbind(sel.inds, elety=blank)  }

       ##- TYPE ATOM/HETATM
       if(type != "") {
         sel.inds <- cbind(sel.inds,
                           elety=is.element(pdb$atom[,"type"],
                             type ) )
       } else {  sel.inds <- cbind(sel.inds, elety=blank)  }

       ##- Take intersection of all seven components
       match.inds <- ( (apply(sel.inds, 1, sum, na.rm=TRUE)==7) )

       if (rm.insert) { # ignore INSERT records
         insert <- which(!is.na(pdb$atom[,"insert"]))
         match.inds[insert] <- FALSE
       }
       ## return XYZ indices
       xyz.inds <- atom2xyz(which(match.inds))

       if (verbose) {
         sel <- rbind( sel, apply(sel.inds, 2, sum, na.rm=TRUE) )
         rownames(sel)=c("Stest","Natom"); ###print(sel)
         cat(paste(" *  Selected a total of:",sum(match.inds),
                   "intersecting atoms  *"),sep="\n")
       }

       match <- list(atom=which(match.inds), xyz=xyz.inds, call = cl)
       class(match) <- "select"
       return(match)
     }
  } #END parse.string


  ##-- Check on input
  if (missing(pdb)) {
    stop("atom.select: must supply 'pdb' input object, e.g. from 'read.pdb'")
  }

  got.string <- TRUE
  got.component <- TRUE
  if(is.null(string)) { got.string <- FALSE }
  if(is.null(c(chain,resno,resid,eleno,elety,type)))  { got.component <- FALSE }

  ## No selection string or component, then just print PDB summary and exit
  if( !any(got.string,  got.component) ) {
    return( print.pdb(pdb) )
  }

  ##
  ##-- Main function
  ##

  ## Replace missing chain " " with NA values to aid selection
  pdb$atom$chain[pdb$atom$chain==" "]= NA

  sel1 <- NULL
  sel2 <- NULL
  if(got.string) {
    if(verbose)
       cat("\n Build selection from input string\n")
    sel1 <- parse.string(pdb, string, verbose, rm.insert)
  }

  if(got.component) {

    ##- Build selection string from input components
    ##   /segid/chain/resno/resid/eleno/elety/
    string <- paste("//",paste(chain, collapse=","),"/",
                    paste(resno, collapse=","),"/",
                    paste(resid, collapse=","),"/",
                    paste(eleno, collapse=","),"/",
                    paste(elety, collapse=","),"/",sep="")
    rm(chain,resno,resid,eleno,elety)
    if(verbose)
       cat("\n Build selection from input components\n")

    ##- Check on 'type="ATOM"'' component
    if(!is.null(type)) {
      if( !all(type %in% c("ATOM", "HETATM")) ) {
        warning("Ignoring input 'type' component (should be one of'ATOM' or 'HETATM' only)")
        type = ""
      }
      } else { type = ""}

    sel2 <- parse.string(pdb, string, verbose, rm.insert, type=type)
  }

  ##- Combine selections from input string and components
  match <- combine.sel(sel1, sel2, op="AND", verbose=verbose)

  ##- Inverse selection
  if(inverse) {
    sele.all <- atom.select(pdb, "all", verbose=verbose)$atom

    sele.new        <- NULL
    sele.new$atom   <- which(!(sele.all %in% match$atom))
    sele.new$xyz    <- atom2xyz(sele.new$atom)
    class(sele.new) <- "select"

    match <- sele.new
    rm(sele.new)

    if(verbose) {
      cat("\n Selection inversed\n")
      cat(paste(" *  Selected a total of:",length(match$atom),
                "atoms  *"), sep="\n")
    }
  }

  return(match)
}
