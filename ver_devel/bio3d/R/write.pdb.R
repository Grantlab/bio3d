"write.pdb" <-
function (pdb = NULL,
          file = "R.pdb",
          xyz = pdb$xyz,
          type = NULL,
          resno = NULL,
          resid = NULL,
          eleno = NULL,
          elety = NULL,
          chain = NULL,
          insert= NULL,
          alt = NULL,
          o = NULL,
          b = NULL,
          segid = NULL,
          elesy = NULL,
          charge = NULL,
          append = FALSE,
          verbose =FALSE,
          chainter = FALSE,
          end = TRUE, 
          print.segid = FALSE) {

  if(is.null(xyz) || !is.numeric(xyz))
    stop("write.pdb: please provide a 'pdb' object or numeric 'xyz' coordinates")

  if(any(is.na(xyz)))
    stop("write.pdb: 'xyz' coordinates must have no NA's.")
  
  if ( is.null(nrow(xyz)) ) {
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

#  card <- rep("ATOM", natom)
  card <- type
  
  if (!is.null(pdb)) {
    if (is.null(card)) card <- pdb$atom$type
    else if(length(card) == 1) card = rep(card, natom)
    if (is.null(resno)) resno = pdb$atom[, "resno"]
    if (is.null(resid)) resid = pdb$atom[, "resid"]
    if (is.null(eleno)) eleno = pdb$atom[, "eleno"]
    if (is.null(elety)) elety = pdb$atom[, "elety"]
    if (is.null(chain)) chain = pdb$atom[, "chain"]
    else if(length(chain) == 1) chain = rep(chain, natom)
    if (is.null(insert)) insert = pdb$atom[, "insert"]
    if (is.null(alt)) alt = pdb$atom[, "alt"]
    if (is.null(o)) o = pdb$atom[, "o"]
    if (is.null(b)) b = pdb$atom[, "b"]

    if (any(is.na(o))) {      o = rep("1.00", natom) }
    if (any(is.na(b))) {      b = rep("0.00", natom) }
    #if (any(is.na(chain))) { chain = rep(" ", natom) }
    chain[is.na(chain)]= ""
    insert[is.na(insert)] = ""
    alt[is.na(alt)] = ""
    
    if (is.null(segid)) segid = pdb$atom[, "segid"]
    segid[is.na(segid)] = ""
 
    if (is.null(elesy)) elesy = pdb$atom[, "elesy"]
    elesy[is.na(elesy)] = ""
   
    if (is.null(charge)) charge = pdb$atom[, "charge"]
  } else {
    if (is.null(card)) card = rep("ATOM", natom)
    else if(length(card) == 1) card = rep(card, natom)
    if (is.null(resno)) resno = c(1:natom)
    if (is.null(resid)) resid = rep("ALA", natom)
    if (is.null(eleno)) eleno = c(1:natom)
    if (is.null(elety)) elety = rep("CA", natom)
    if (is.null(chain)) chain = rep("", natom)
    else if (length(chain) == 1) chain = rep(chain, natom)
    ##if(any(is.na(chain))) chain[is.na(chain)]= ""
    if (is.null(insert)) insert=rep("", natom)
    if (is.null(alt)) alt=rep("", natom)
    if (is.null(o))         o = rep("1.00",natom)
    if (is.null(b))         b = rep("0.00", natom)
    if (is.null(segid)) segid = rep("", natom)
    if (is.null(elesy)) elesy = rep("", natom)
    if (is.null(charge)) charge = rep("", natom)

    chain[is.na(chain)]= ""
    insert[is.na(insert)] = ""
    alt[is.na(alt)] = ""
  }

  if (!is.logical(append)) 
    stop("write.pdb: 'append' must be logical TRUE/FALSE")
  
  if (length(as.vector(xyz))%%3 != 0) {
    stop("write.pdb: 'length(xyz)' must be divisable by 3.")
  }
  check.lengths <- sum(length(card), length(resno), length(resid), 
                       length(eleno), length(elety), length(chain), 
                       length(insert), length(alt), length(o), length(b), 
                       length(segid), length(elesy), length(charge))
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
    
    format <- "%-6s%5s  %-3s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s"
    if (nchar(elety) > 3) {
#    if (nchar(elety) >= 3) {
#      if ((substr(elety, 2, 2) == "H") | (substr(elety, 1, 1) == "H")) {
        format <- "%-6s%5s %-4s%1s%-4s%1s%4s%1s%3s%8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s%2s%2s"
#      }
    }
    sprintf(format, card, eleno, elety, alt, resid, chain,
            resno, insert, "", x, y, z, o, b, "", segid, elesy, charge)
  } 

  ##### To write SSE annotations #####
  format.sse <- function(pdb) {
    # format 'sse' component in a pdb object to include resid and length
    lapply(c(helix='helix', sheet='sheet'), function(x) {
       sse <- pdb[[x]]
       if(length(sse$start) > 0) {
         ref <- pdb$atom[pdb$calpha, c("resno", "chain", "insert", "resid")]
         refkey <- paste(ref$resno, ref$chain, ref$insert, sep = "_")
         insert <- sub(' +', '', names(sse$start))
         insert[insert == ''] <- NA
         skey <- paste(sse$start, sse$chain, insert, sep="_")
         insert <- sub(' +', '', names(sse$end))
         insert[insert == ''] <- NA
         ekey <- paste(sse$end, sse$chain, insert, sep="_")
         resno = data.frame(start = sse$start, end = sse$end)
         resid = data.frame(start = ref$resid[match(skey, refkey)], 
                              end = ref$resid[match(ekey, refkey)])
         insert = data.frame(start= names(sse$start), end = names(sse$end))
         type = ifelse(x=='helix', 'type', 'sense')
         length = match(ekey, refkey) - match(skey, refkey) + 1
         return(list(resno=resno, resid=resid, chain=sse$chain, insert=insert,
                     type=sse[[type]], length=length))
       } else {
         return (NULL)
       }
     } )
  }

  if(!is.null(pdb)) {
    lines <- NULL
    sse <- format.sse(pdb)
    if(!is.null(sse$helix)) {
       format <- "%-6s %3d %3s %3s %1s %4d%1s %3s %1s %4d%1s%2s%30s %5d"
       for(i in 1:nrow(sse$helix$resno))
         lines <- rbind(lines, sprintf(format, 'HELIX', i, i, sse$helix$resid$start[i],
               sse$helix$chain[i], sse$helix$resno$start[i], sse$helix$insert$start[i],
               sse$helix$resid$end[i], sse$helix$chain[i], sse$helix$resno$end[i],
               sse$helix$insert$end[i], sse$helix$type[i], '', sse$helix$length[i])) 
    }
    if(!is.null(sse$sheet)) {
       format <- "%-6s %3d %3s%2d %3s %1s%4d%1s %3s %1s%4d%1s%2s"
       for(i in 1:nrow(sse$sheet$resno))
         lines <- rbind(lines, sprintf(format, 'SHEET', i, 'S1', nrow(sse$sheet$resno), 
               sse$sheet$resid$start[i], sse$sheet$chain[i], sse$sheet$resno$start[i], 
               sse$sheet$insert$start[i], sse$sheet$resid$end[i], sse$sheet$chain[i], 
               sse$sheet$resno$end[i], sse$sheet$insert$end[i], sse$sheet$type[i]) )
    }
    if(!is.null(lines)) {
       write.table(lines, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = append)
       if(!append) append = TRUE
    }
  }
  ##################################### 

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
                                       alt = alt[i],
                                       resid = resid[i],
                                       chain = chain[i],
                                       resno = resno[i],
                                       insert = insert[i],
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
               "TER", as.character(teleno[i] + ii), "", resid[i], "", chain[i], resno[i], insert[i]))
        }
      }
          
    }
    ## Changed cat() for write.table() as sugested by Joao Martins <joao.martins@env.ethz.ch>
    ##cat(lines, file = file, sep = "\n", append = append)
    write.table(lines, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = append)
    if(chainter) {
          ii = ii + 1
          cat(sprintf("%-6s%5s%6s%3s%1s%1s%4s%1s\n", "TER", as.character(teleno[i] + ii), "", 
              resid[i], "", chain[i], resno[i], insert[i]), file = file, append = TRUE)
    }
    if(end) {
      cat("END   \n", file = file, append = TRUE)
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
                                         alt = alt[i],
                                         resid = resid[i],
                                         chain = chain[i],
                                         resno = resno[i],
                                         insert = insert[i],
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
               "TER", as.character(teleno[i] + ii), "", resid[i], "", chain[i], resno[i], insert[i]))
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
              resid[i], "", chain[i], resno[i], insert[i]), file=file, append=TRUE)
      }
      cat(sprintf("%-6s\n", "ENDMDL"), file = file, append = TRUE)
    }
    if(end) {
      cat("END   \n", file = file,  append = TRUE)
    }
  }
  ## Write connectivity
  if(!is.null(pdb$con)) {
    pdb$con[is.na(pdb$con)] <- ""
    eleno.1 <- pdb$con$eleno.1
    eleno.1 <- split(eleno.1, eleno.1)
    eleno.1 <- unlist(lapply(eleno.1,
                             function(object)
                               rep(unique(object), ceiling(length(object)/4))))
    eleno.2 <- split(pdb$con$eleno.2, pdb$con$eleno.1)
    eleno.2 <- lapply(eleno.2,
                      function(object) {
                        v <- vector("character", 4*ceiling(length(object)/4))
                        v[1:length(object)] <- object
                        M <- matrix(v, ncol = 4, nrow = ceiling(length(object)/4))
                        return(M)
                      })
    eleno.2 <- as.data.frame(do.call(rbind, eleno.2))
    fmt <- "CONECT%5s%5s%5s%5s%5s"
    args  <- c(fmt, as.list(cbind(eleno.1, eleno.2)))
    lines <- do.call(sprintf, args)
    cat(lines, file = file, sep = "\n", append = TRUE)
  }
  ## End write connectivity
  if (verbose) cat(" DONE","\n")
}

