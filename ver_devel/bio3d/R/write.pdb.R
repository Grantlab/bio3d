"write.pdb" <-
function (pdb = NULL,
          file = "R.pdb",
          xyz = pdb$xyz,
          resno = NULL,
          resid = NULL,
          eleno = NULL,
          elety = NULL,
          chain = NULL,
          o = NULL,
          b = NULL,
          segid = NULL,
          elesy = NULL,
          charge = NULL,
          het = FALSE,
          append = FALSE,
          verbose =FALSE,
          chainter = FALSE,
          end = TRUE, 
          print.segid = FALSE) {

  if(is.null(xyz) || !is.numeric(xyz))
    stop("write.pdb: please provide a 'pdb' object or numeric 'xyz' coordinates")

  if(any(is.na(xyz)))
    stop("write.pdb: 'xyz' coordinates must have no NA's.")
  
  if (is.vector(xyz)) {
    natom <- length(xyz)/3
    nfile <- 1
  } else if (is.matrix(xyz)) {
    natom <- ncol(xyz)/3 
    nfile <- nrow(xyz)
    if (verbose) {
      cat("Multiple rows in 'xyz' will be interperted as multimodels/frames\n")
    }
  } else {
    stop("write.pdb: 'xyz' or 'pdb$xyz' must be either a vector or matrix")
  }

  card <- rep("ATOM", natom)
  
  if (!is.null(pdb)) {
    if(natom == 1)
      ## make sure we are a matrix
      pdb$atom <- t(as.matrix(pdb$atom))
    
    if (het) 
      card <- c( rep("ATOM", nrow(pdb$atom)), rep("HETATM", nrow(pdb$het)) )
    if (is.null(resno)) {
      resno = pdb$atom[, "resno"]
      if (het) { resno = c(resno, pdb$het[, "resno"]) }}
    
    if (is.null(resid)) {
      resid = pdb$atom[, "resid"]
      if (het) { resid = c(resid, pdb$het[, "resid"]) }}
    
    if (is.null(eleno)) {
      eleno = pdb$atom[, "eleno"]
      if (het) { eleno = c(eleno, pdb$het[, "eleno"]) }}
    
    if (is.null(elety)) {
      elety = pdb$atom[, "elety"]
      if (het) { elety = c(elety, pdb$het[, "elety"]) }}
    
    if (is.null(chain)) {
      chain = pdb$atom[, "chain"]
      if (het) { chain = c(chain, pdb$het[, "chain"]) }}
    
    if (is.null(o)) {
      o = pdb$atom[, "o"]
      if (het) { o = c(o, pdb$het[, "o"]) }}
    
    if (is.null(b)) {
      b = pdb$atom[, "b"]
      if (het) { b = c(b, pdb$het[, "b"]) }}
     
    if (any(is.na(o))) {      o = rep("1.00", natom) }
    if (any(is.na(b))) {      b = rep("0.00", natom) }
    #if (any(is.na(chain))) { chain = rep(" ", natom) }
    chain[is.na(chain)]= ""
    
    if (is.null(segid)) {
       segid = pdb$atom[, "segid"]
       if (het) { segid = c(segid, pdb$het[, "segid"]) }}
    segid[is.na(segid)] = ""
 
    if (is.null(elesy)) {
       elesy = pdb$atom[, "elesy"]
       if (het) { elesy = c(elesy, pdb$het[, "elesy"]) }}
    elesy[is.na(elesy)] = ""
   
    if (is.null(charge)) {
       charge = pdb$atom[, "charge"]
       if (het) { charge = c(charge, pdb$het[, "charge"]) }}
   
  } else {
    if (is.null(resno)) resno = c(1:natom)
    if (is.null(resid)) resid = rep("ALA", natom)
    if (is.null(eleno)) eleno = c(1:natom)
    if (is.null(elety)) elety = rep("CA", natom)
    if (is.null(chain)) chain = rep("", natom)
    if (is.null(o))         o = rep("1.00",natom)
    if (is.null(b))         b = rep("0.00", natom)
    if (is.null(segid)) segid = rep("", natom)
    if (is.null(elesy)) elesy = rep("", natom)
    if (is.null(charge)) charge = rep("", natom)
  }

  
  if (!is.logical(append)) 
    stop("write.pdb: 'append' must be logical TRUE/FALSE")
  
  if (length(as.vector(xyz))%%3 != 0) {
    stop("write.pdb: 'length(xyz)' must be divisable by 3.")
  }
  check.lengths <- sum(length(resno), length(resid), length(eleno),
                       length(elety), length(o), length(b), length(segid), 
                       length(elesy), length(charge))
  if (check.lengths%%natom != 0) {
    stop("write.pdb: the lengths of all input vectors != 'length(xyz)/3'.")
  }

  
  o <- as.numeric(o)
  b <- as.numeric(b)
  eleno <- as.character(eleno)
  resno <- as.character(resno)
  charge <- as.character(charge)
  charge[is.na(charge)] = ""
  if(!print.segid) segid = rep("", natom)
  ## Inserted Jul 8th 2008 for adding TER between chains
  ter.lines <- (which(!duplicated(chain))[-1] - 1)


  
  atom.print <- function(card = "ATOM", eleno, elety, alt = "",
        resid, chain = "", resno, insert = "", x, y, z, o = "1.00",
        b = "0.00", segid = "", elesy = "", charge = "") {
    
    format <- "%-6s%5s  %-3s%1s%-3s%1s%1s%4s%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s"
    if (nchar(elety) > 3) {
#    if (nchar(elety) >= 3) {
#      if ((substr(elety, 2, 2) == "H") | (substr(elety, 1, 1) == "H")) {
        format <- "%-6s%5s %-4s%1s%-3s%1s%1s%4s%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s"
#      }
    }
    sprintf(format, card, eleno, elety, alt, resid, "", chain,
            resno, insert, "", x, y, z, o, b, "", segid, elesy, charge)
  } 


  
  if(nfile==1) {
    coords <- matrix(round(as.numeric(xyz), 3), ncol = 3, byrow = TRUE)
    if (verbose) {
      cat(paste("Writing 1 frame with",natom,"atoms "))
    }
    lines <- NULL
    ii = 0
    teleno <- as.numeric(eleno)
    for (i in 1:natom) {
       
      lines <- rbind(lines, atom.print( card = card[i],
                                       eleno = as.character(teleno[i] + ii),
                                       elety = elety[i],
                                       resid = resid[i],
                                       chain = chain[i],
                                       resno = resno[i],
                                       x = coords[i, 1],
                                       y = coords[i, 2],
                                       z = coords[i, 3],
                                       o = o[i], b = b[i],
                                       segid = segid[i],
                                       elesy = elesy[i],
                                       charge = charge[i]))
      
      ## Inserted Jul 8th 2008 for adding TER between chains 
      ## Modified to be consistent to PDB format v3.3
      if(chainter) {
        if(i %in% ter.lines) {
#          lines <- rbind(lines, "TER   ")
          ii = ii + 1
          lines <- rbind(lines, sprintf("%-6s%5s%6s%3s%1s%1s%4s%1s", 
               "TER", as.character(teleno[i] + ii), "", resid[i], "", chain[i], resno[i], ""))
        }
      }
          
    }
    ## Changed cat() for write.table() as sugested by Joao Martins <joao.martins@env.ethz.ch>
    ##cat(lines, file = file, sep = "\n", append = append)
    write.table(lines, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = append)
    if(chainter) {
          ii = ii + 1
          cat(sprintf("%-6s%5s%6s%3s%1s%1s%4s%1s\n", "TER", as.character(teleno[i] + ii), "", 
              resid[i], "", chain[i], resno[i], ""), file = file, append = TRUE)
    }
    if(end) {
      cat("END   ", file = file, append = TRUE)
    }
    
  } else {
    if (verbose) {
      cat(paste("Writing",nfile,"frames with",natom,"atoms"),"\n")
      cat("Frame Progress (x50) ")
    }
    if(!append) unlink(file) 
    for (j in 1:nfile) {
      coords <- matrix(round(as.numeric(xyz[j,]), 3), ncol = 3, byrow = TRUE)
      lines <- NULL
      ii = 0
      teleno <- as.numeric(eleno)
      for (i in 1:natom) {
        lines <- rbind(lines, atom.print( eleno = as.character(teleno[i] + ii),
                                         elety = elety[i],
                                         resid = resid[i],
                                         chain = chain[i],
                                         resno = resno[i],
                                         x = coords[i, 1],
                                         y = coords[i, 2],
                                         z = coords[i, 3],
                                         o = o[i], b = b[i],
                                         segid = segid[i],
                                         elesy = elesy[i],
                                         charge = charge[i]))

        ## Inserted Jul 8th 2008 for adding TER between chains (untested) 
        ## Modified to be consistent to PDB format v3.3
        if(chainter) {
          if(i  %in% ter.lines) {
#            lines <- rbind(lines, "TER   ")
            ii = ii + 1
            lines <- rbind(lines, sprintf("%-6s%5s%6s%3s%1s%1s%4s%1s", 
               "TER", as.character(teleno[i] + ii), "", resid[i], "", chain[i], resno[i], ""))
          }
        }
        
      }
      if (verbose) {
        if (j%%50 == 0) cat(".")
      }
      ##cat(lines, file = file, sep = "\n", append = TRUE)
      cat(sprintf("%-6s%4s%4d\n", "MODEL", " ", j), file = file, append = TRUE)
      write.table(lines, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      if(chainter) {
         ii = ii + 1
         cat(sprintf("%-6s%5s%6s%3s%1s%1s%4s%1s\n", "TER", as.character(teleno[i] + ii), "", 
              resid[i], "", chain[i], resno[i], ""), file=file, append=TRUE)
      }
      cat(sprintf("%-6s\n", "ENDMDL"), file = file, append = TRUE)
    }
    if(end) {
      cat("END   ", file = file,  append = TRUE)
    }

  }
  if (verbose) cat(" DONE","\n")
}

