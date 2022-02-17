"convert.pdb" <-
function(pdb, type = c("original", "pdb", "charmm", "amber", "gromacs"),
         renumber=FALSE, first.resno=1, first.eleno=1,
         consecutive=TRUE,
         rm.h=TRUE, rm.wat=FALSE, verbose=TRUE) {

  ##-- Check the requested output format is one of 'type.options' 
  type <- match.arg(type)

  ##-- Water and hydrogen removal
  inds <- NULL
  if(rm.wat) {
    inds <- atom.select(pdb, "notwater", verbose=FALSE)
    if(verbose){
      cat(paste("\t Retaining", length(inds$atom),"non-water atoms\n"))
    }
  }
  if(rm.h) {
    inds <- combine.select(inds, atom.select(pdb, "noh", verbose=FALSE), verbose=FALSE) 
    if(verbose){
      cat(paste("\t Retaining", length(inds$atom),"non-hydrogen atoms\n"))
    } 
  }
  if(!is.null(inds)){
    nrm <- nrow(pdb$atom) - length(inds$atom)
    if( nrm > 0) {  
      if(verbose){
        cat(paste("\t Removing a total of", nrm," atoms\n"))
      } 
      pdb <- trim.pdb(pdb, inds)
    }
  }

  
  ##-- Renumbering of residues and atoms
  if(renumber) {
    if(verbose){
      cat(paste("\t Renumbering residues ( from",first.resno,") and atoms ( from",first.eleno,")\n"))
    } 

    ## Assign consecutive atom numbers 
    pdb$atom[,"eleno"] <- seq(first.eleno, length=nrow(pdb$atom))

    ## Determine chain start and end indices
    s.ind <- which(!duplicated(pdb$atom[,"chain"]))
    e.ind   <- c(s.ind[-1]-1, nrow(pdb$atom))

    ##- Assign new (consecutive) residue numbers for each chain
    prev.chain.res = 0  ## Number of residues in previous chain
    for (i in 1:length(s.ind)) {
      ## Combination of resno and insert code define a residue (wwpdb.org)
      insert = pdb$atom[s.ind[i]:e.ind[i], "insert"]
      insert[is.na(insert)] = ""
      resno0 <- paste0(pdb$atom[s.ind[i]:e.ind[i], "resno"], insert)

      ## Ordered table of residue occurrences
      tbl <- table(resno0)[unique(resno0)]
      n.chain.res <- length(tbl)


      new.nums <- (first.resno+prev.chain.res):(first.resno+n.chain.res-1+prev.chain.res)
      pdb$atom[s.ind[i]:e.ind[i],"resno"] <- rep(new.nums, tbl)

      ## SSE
      if(length(pdb$helix)>0) {
#         chs = unique(pdb$helix$chain)
         ch = pdb$atom$chain[s.ind][i]
         if(ch %in% pdb$helix$chain) {
            h.start <- pdb$helix$start[pdb$helix$chain %in% ch]
            h.insert <- names(h.start)
            t.inds = match(paste0(h.start, h.insert), unique(resno0))
            pdb$helix$start[pdb$helix$chain %in% ch] = new.nums[t.inds]

            h.end <- pdb$helix$end[pdb$helix$chain %in% ch]
            h.insert <- names(h.end)
            t.inds = match(paste0(h.end, h.insert), unique(resno0))
            pdb$helix$end[pdb$helix$chain %in% ch] = new.nums[t.inds]
         }
      }

      if(length(pdb$sheet)>0) {
#         chs = unique(pdb$sheet$chain)
         ch = pdb$atom$chain[s.ind][i]
         if(ch %in% pdb$sheet$chain) {
            s.start <- pdb$sheet$start[pdb$sheet$chain %in% ch]
            s.insert <- names(s.start)
            t.inds = match(paste0(s.start, s.insert), unique(resno0))
            pdb$sheet$start[pdb$sheet$chain %in% ch] = new.nums[t.inds]
            
            s.end <- pdb$sheet$end[pdb$sheet$chain %in% ch]
            s.insert <- names(s.end)
            t.inds = match(paste0(s.end, s.insert), unique(resno0))
            pdb$sheet$end[pdb$sheet$chain %in% ch] = new.nums[t.inds]
         }
      }
    
      if(consecutive) {
        ## Update prev.chain.res for next iteration 
        prev.chain.res = prev.chain.res + n.chain.res
      }
    }
    
    ## clean up 'insert'
    pdb$atom$insert <- as.character(NA)
    if(length(pdb$helix)>0) {
      names(pdb$helix$start) <- NULL
      names(pdb$helix$end) <- NULL
    }
    if(length(pdb$sheet)>0) {
      names(pdb$sheet$start) <- NULL
      names(pdb$sheet$end) <- NULL
    }
  }


  ##-- Format conversion
  if(type != "original") {
    if(verbose){
      cat(paste0("\t Converting to '", type, "' format\n"))
    }
    ## residue and atom types from PDB
    restype <- unique(pdb$atom[,"resid"])

    #eletype <- unique(pdb$atom[,"elety"])
    ## In future could determine 'input type' based on resid/elety

    ##- Check for non-standard residue names
    if(verbose){
      not.prot.inds <- atom.select(pdb, "notprotein", verbose=FALSE)$atom
      if(length(not.prot.inds) > 0) { 
        not.prot.res <- paste(unique(pdb$atom[not.prot.inds, "resid"]), collapse = " ")
        cat(paste("\t Non-standard residue names present (",not.prot.res,")\n") )
      }
    }

    ##- Convert HIS resid
    his <- matrix( c("HIS", "HSD", "HID","HISA",
                     "HIS", "HSE", "HIE","HISB",
                     "HIS", "HSP", "HIP","HISH"),
                  nrow=3, byrow=TRUE,
                  dimnames = list(c("d","e","b"),
                    c("pdb","charmm","amber","gromacs")) )
  
    type.inds <-  (colnames(his) %in% type)
    conv.inds <- !(colnames(his) %in% c(type,"pdb"))

    his.d.ind <- (pdb$atom[,"resid"] %in% his["d", !type.inds ])
    his.e.ind <- (pdb$atom[,"resid"] %in% his["e", conv.inds ])
    his.b.ind <- (pdb$atom[,"resid"] %in% his["b", conv.inds ])
  
    pdb$atom[his.d.ind,"resid"] <- his["d", type.inds ]
    pdb$atom[his.e.ind,"resid"] <- his["e", type.inds ]
    pdb$atom[his.b.ind,"resid"] <- his["b", type.inds ]

    ##- Convert ILE CD1 to CD elety and remove chainID
    if (type=="charmm") {
      ile.ind <- atom.select(pdb, resid="ILE", elety="CD1", verbose=FALSE)$atom
      pdb$atom[ile.ind,"elety"] <- "CD"

      pdb$atom[,"chain"]=NA ## strip chain ID

      ## Could also add a SEGID via call to chain.pdb() function
      pdb$atom[,"segid"] <- chain.pdb(pdb)

    } else {
      ile.ind <- atom.select(pdb, resid="ILE", elety="CD", verbose=FALSE)$atom
      pdb$atom[ile.ind,"elety"] <- "CD1"
    }
  } ## END type != "original" (conversion)
  

  ##-- Convert hydrogen atom types (unfinished!)
  if(!rm.h) { 
    if(type=="pdb") {  
      pdb$atom[ pdb$atom[,"elety"]=="HN", "elety"] = "H"
      ###!!! ADD Many MORE ATOM TYPE CONVERSIONS HERE !!!###
    }
    if(verbose){
      warning(paste("\t Additional hydrogen elety names may need converting.",
                    "\t   N.B. It is often best to remove hydrogen (rm.h=TRUE)", 
                    "\t        before building systems for simulation",sep="\n"))
    }  
    ## Add other atom name conversions here as the need arises...

  } 

  return(pdb)
}

