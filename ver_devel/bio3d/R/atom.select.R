"atom.select" <-
function(pdb, string=NULL,
         chain=NULL, resno=NULL, resid=NULL,
         eleno=NULL, elety=NULL,
         verbose=TRUE, rm.insert=FALSE) {
  ## Version 0.6 ... Wed Jul  3 12:25:18 EDT 2013
  ##
  ##  Include the function to take the intersection
  ##  of string and component (Xinqiu)
  ##
  ## Version 0.5 ... Fri Nov 30 13:10:24 EST 2012
  ##
  ##  Re-worte large parts of the keyword matching
  ##   code to make it slightly more transparent
  ##
  ## Version 0.4 ... Thu Nov 29 15:27:03 EST 2012
  ##
  ##  Added string="protein" and a few other keywords
  ##
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
  ## Version 0.1 ... Tue Mar 21 10:58:43 PST 2006
  ##
  ##   Prints a summary of 'pdb' makeup if called 
  ##    without a selection 'string'.
  ##   Also added 'string' shortcuts "calpha", 
  ##    "back", "backbone" and "cbeta"
  ##
  ## Version 0.0 ... Fri Mar 17 14:45:37 PST 2006
  ##
  ##  Description:
  ##   Atom selection function
  ##    Losely based on PyMol selection Macro:
  ##   see: http://www.pymolwiki.org/index.php/Selection_Macros
  ##
  ##   String Selection Syntax:
  ##    "//A/130:142///N,CA,C,O/"
  ##    "/segid/chain/resno/resid/eleno/elety/"
  ##
  ##  E.g.
  ##   # read a PDB file
  ##   pdb<-read.pdb("1bg2.pdb")
  ##   # print a structure summary
  ##   atom.select(pdb)
  ##   # select all C-alpha atoms from resno 65 to 143
  ##   ca.inds   <- atom.select(pdb, "///65:143///CA/")
  ##   # or all C-alphas
  ##   ca.inds   <- atom.select(pdb, "calpha")
  ##   # more examples
  ##   inds<-atom.select(pdb, "//A/130:142///N,CA,C,O/")


  
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
  
  sel.txt2type <- function(type.sel.txt) {
    ##- Splitting function for characters - split on coma & remove white space
    return( gsub(" ","", unlist(strsplit(type.sel.txt, split=",")) ) )
  }
  
  ##-- Parse string and return the selection
  parse.string <- function(pdb, string, verbose, rm.insert) {

     ##-- We have input selection string
     aa <- unique(pdb$atom[,"resid"])
     prot.aa <- c("ALA", "CYS", "ASP", "GLU", 
                  "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", 
                  "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR", 
                  "SEP", "TPO", "MLY", "MSE", "IAS", "ABA", "CSO", "CSD", 
                  "CYM", "CME", "CSX", "CMT", "CYX", "HIE", "HIP", "HID", 
                  "HSD", "HSE", "HSP", "DDE", "MHO", "ASX", "CIR", "PFF")
     
     hoh <- c("H2O", "OH2", "HOH", "HHO", "OHH", "SOL", "WAT",
              "TIP", "TIP2", "TIP3", "TIP4")

     not.prot.aa <- (aa[!aa %in% prot.aa])
     ligand.aa <- not.prot.aa[!not.prot.aa %in% hoh]
     if(length(not.prot.aa) == 0) { not.prot.aa = "xxxxx" }
     if(length(ligand.aa) == 0) { ligand.aa = "xxxxx" }
     
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
                 calpha = "//////CA/",
                 cbeta = "//////N,CA,C,O,CB/",
                 backbone = "//////N,CA,C,O/",
                 back = "//////N,CA,C,O/",
                 all = "///////",
                 protein = paste("////",paste(prot.aa, collapse=","), "///",sep=""),
                 notprotein = paste("////", paste(not.prot.aa, collapse=","), "///",sep=""),
                 ligand = paste("////", paste(ligand.aa, collapse=","), "///",sep=""),
                 water = paste("////",paste(hoh, collapse=","), "///",sep=""),
                 notwater = paste("////", paste(aa[!aa %in% hoh], collapse=","), "///",sep=""),
                 #noh = "!H*",
                 #h = "H*",
                 NA)
     
     if(is.na(i)) {
       ##- No string shortcut match
       
       if(!substr(string,1,1)=="/") {
         ## Check if we have a valid selection sting
         stop("Not a valid selection string shortcut.\n\t Please use one of:
        'calpha' 'cbeta' 'backbone'
        'protein' 'notprotein' 'ligand'
        'water' 'notwater'
        'h' 'noh'\n
       Or valid selection string:
         /segid/chain/resno/resid/eleno/elety/ \n")
       }

     } else {
       ##- Use string shortcut from switch function
       if(verbose)
         cat(paste("\n Using selection 'string' keyword shortcut:",string, "=", i, "\n\n"))
       string = i
     }
     
     sel <- unlist(strsplit(string, split = "/"))
   
     if (sel[1] == "") { # full selection string (starts with "/")
       sel <- sel[-1]
       if(length(sel) != 6) {
         print("missing elements, should be:\n/segid/chain/resno/resid/eleno/elety/")
       }
   
       names(sel) <- c("segid","chain","resno","resid","eleno","elety")
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
                             sel.txt2type( sel["chain"] )) )        
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
   
       match.inds <- ( (apply(sel.inds, 1, sum, na.rm=TRUE)==6) )
       ## In future, could take the inverse here for NOT selection
       
       if (rm.insert) { # ignore INSERT records
         insert <- which(!is.na(pdb$atom[,"insert"]))
         match.inds[insert] <- FALSE
       }
       ## return XYZ indices
       #xyz.inds <- matrix(1:length( pdb$atom[,c("x","y","z")] ),nrow=3,byrow=FALSE)
       #xyz.inds <- as.vector(xyz.inds[,match.inds])
       xyz.inds <- atom2xyz(which(match.inds))
   
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
  } #END parse.string


  ##-- Check on input
  if (missing(pdb)) {
    stop("atom.select: must supply 'pdb' input object, e.g. from 'read.pdb'")
  }

  got.string <- TRUE
  got.component <- TRUE
  if(is.null(string)) { got.string <- FALSE }
  if(is.null(c(chain,resno,resid,eleno,elety)))  { got.component <- FALSE }

  ## No selection string or component, then just print PDB summary
  if( !any(got.string,  got.component) ) {
    return( print.pdb(pdb) )
  }

## Modified for combining string and component, (Jul 3, 2013)
#  ## Selection string and component, then component is ignored
#  if(all(got.string, got.component)) {
#    warning("Selection string AND selection component given: using first input string only!")
#    ## Could Change this to combine string and component here!!
#  } else {

  ##-- Main function  

  sel1 <- NULL
  sel2 <- NULL
  if(got.string) {
    if(verbose) 
       cat("\nBuild selection from input string\n\n")
    sel1 <- parse.string(pdb, string, verbose, rm.insert)
  }

  if(got.component) {

    ## Build selection string from input components
    ## /segid/chain/resno/resid/eleno/elety/
    string <- paste("//",paste(chain, collapse=","),"/",
                    paste(resno, collapse=","),"/",
                    paste(resid, collapse=","),"/",
                    paste(eleno, collapse=","),"/",
                    paste(elety, collapse=","),"/",sep="")
    rm(chain,resno,resid,eleno,elety)
    if(verbose) 
       cat("\nBuild selection from input components\n\n")
    sel2 <- parse.string(pdb, string, verbose, rm.insert)
  }

  if(!is.null(sel1) && !is.null(sel2))
     if(verbose)
        cat("\nCombine selections from input string and components\n\n")

  match <- combine.sel(sel1, sel2, op="AND", verbose=verbose)
  
  return(match)
}

