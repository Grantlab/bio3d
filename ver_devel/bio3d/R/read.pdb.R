"read.pdb" <-
function (file, maxlines=50000, multi=FALSE,
          rm.insert=FALSE, rm.alt=TRUE, het2atom=FALSE, verbose=TRUE) {

  if(missing(file)) {
    stop("read.pdb: please specify a PDB 'file' for reading")
  }
  if(!is.numeric(maxlines)) {
    stop("read.pdb: 'maxlines' must be numeric")
  }
  if(!is.logical(multi)) {
    stop("read.pdb: 'multi' must be logical TRUE/FALSE")
  }

  ## Check if file exists localy or online
  toread <- file.exists(file)
  if(substr(file,1,4)=="http") { toread <- TRUE }

  ## Check for 4 letter code and possible online file
  if(!toread) {
    if(nchar(file)==4) {
      file <- get.pdb(file, URLonly=TRUE)
      cat("  Note: Accessing online PDB file\n")
    } else {
      stop("No input PDB file found: check filename")
    }
  }
  
  # PDB FORMAT v2.0:    colpos,  datatype,    name,      description
  atom.format <- matrix(c(-6,     NA,          NA,       # (ATOM)
                          5,     'numeric',   "eleno",   # atom_no
                         -1,     NA,          NA,        # (blank)
                          4,     'character', "elety",   # atom_ty
                          1,     'character', "alt",     # alt_loc
                          4,     'character', "resid",   # res_na 
                          1,     'character', "chain",   # chain_id 
                          4,     'numeric',   "resno",   # res_no
                          1,     'character', "insert",  # ins_code
                         -3,     NA,           NA,       # (blank)
                          8,     'numeric',   "x",       # x
                          8,     'numeric',   "y",       # y
                          8,     'numeric',   "z",       # z
                          6,     'numeric',   "o",       # o
                          6,     'numeric',   "b",       # b
                         -6,     NA,           NA,       # (blank)
                          4,     'character', "segid"    # seg_id
                         ), ncol=3, byrow=TRUE,
                       dimnames = list(c(1:17), c("widths","what","name")) )

  split.string <- function(x) {
    # split a string 'x'
    x <- substring(x, first, last)
    x[nchar(x) == 0] <- as.character(NA)
    x
  }
  is.character0 <- function(x){length(x)==0 & is.character(x)}
  
  trim <- function (s) {
    # Remove leading and traling
    # spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-NA
    s
  }


  # finds first and last (substr positions)
  widths <-  as.numeric(atom.format[,"widths"]) # fixed-width spec  
  drop.ind <- (widths < 0) # cols to ignore (i.e. -ve)
  widths <- abs(widths)    # absolute vales for later  
  st <- c(1, 1 + cumsum( widths ))
  first <- st[-length(st)][!drop.ind] # substr start
  last <- cumsum( widths )[!drop.ind] # substr end

  # read n lines of PDB file
  raw.lines  <- readLines(file, n = maxlines)
  type <- substring(raw.lines,1,6)

  # check number of END/ENDMDL records
  raw.end <- sort(c(which(type == "END"),
                    which(type == "ENDMDL")))
  
  if (length(raw.end) > 1) {
    print("PDB has multiple END/ENDMDL records")
    if (!multi) {
      print("multi=FALSE: taking first record only")
      raw.lines <- raw.lines[ (1:raw.end[1]) ]
      type <- type[ (1:raw.end[1]) ]
    } else {
      print("multi=TRUE: 'read.dcd' will be quicker!")
    }
  }
  if ( length(raw.end) !=1 ) {
    if (length(raw.lines) == maxlines) {
      # have not yet read all the file
      print("You may need to increase 'maxlines'")
      print("check you have all data in $atom")
    }
  }

  # split by record type
  raw.header <- raw.lines[type == "HEADER"]
  raw.seqres <- raw.lines[type == "SEQRES"]
  raw.helix  <- raw.lines[type == "HELIX "]
  raw.sheet  <- raw.lines[type == "SHEET "]
  raw.atom   <- raw.lines[type == "ATOM  "]
  het.atom   <- raw.lines[type == "HETATM"]
  all.atom   <- raw.lines[type %in% c("ATOM  ","HETATM")]
  
  # also look for "TER" records
  rm(raw.lines)
  
  if (verbose) {
    if (!is.character0(raw.header)) { cat(" ", raw.header, "\n") }
  }
  ## Edit from Baoqiang Cao <bqcao@ices.utexas.edu> Nov 29, 2009
  ## Old version:
  ##  seqres <- unlist(strsplit( trim(substring(raw.seqres,19,80))," +"))
  ## New version

  seqres <- unlist(strsplit( trim(substring(raw.seqres,19,80))," +"))
  if(!is.null(seqres)) {
    seqres.ch <- substring(raw.seqres, 12, 12)
    seqres.ln <- substring(raw.seqres, 13, 17)
    seqres.in <- ( !duplicated(seqres.ch) )
    names(seqres) <- rep(seqres.ch[seqres.in], times=seqres.ln[seqres.in])
  }
  ## End Edit from Baoqiang:

  
  helix  <- list(start = as.numeric(substring(raw.helix,22,25)),
                 end   = as.numeric(substring(raw.helix,34,37)),
                 chain = trim(substring(raw.helix,20,20)),
                 type  = trim(substring(raw.helix,39,40)))

  sheet  <- list(start = as.numeric(substring(raw.sheet,23,26)),
                 end   = as.numeric(substring(raw.sheet,34,37)),
                 chain = trim(substring(raw.sheet,22,22)),
                 sense = trim(substring(raw.sheet,39,40)))
  
  # format ATOM records as a character matrix
  if (het2atom) {
    atom <- matrix(trim(sapply(all.atom, split.string)), byrow=TRUE,
                   ncol=nrow(atom.format[ !drop.ind,]), 
                   dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )
  } else {
    atom <- matrix(trim(sapply(raw.atom, split.string)), byrow=TRUE,
                   ncol=nrow(atom.format[ !drop.ind,]), 
                   dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )
  }
  
  # Alt records with m[,"alt"] != NA
  if (rm.alt) {
    if ( sum( !is.na(atom[,"alt"]) ) > 0 ) {
      ## Edited: Mon May  4 17:41:11 PDT 2009 to cope with
      ## both numeric and character ALT records
      first.alt <- sort( unique(na.omit(atom[,"alt"])) )[1]
      cat(paste("   PDB has ALT records, taking",first.alt,"only, rm.alt=TRUE\n"))
      alt.inds <- which( (atom[,"alt"] != first.alt) ) # take first alt only
      if(length(alt.inds)>0)
        atom <- atom[-alt.inds,]
    }
  }
  # Insert records with m[,"insert"] != NA
  if (rm.insert) {
    if ( sum( !is.na(atom[,"insert"]) ) > 0 ) {
      cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
      insert.inds <- which(!is.na(atom[,"insert"])) # rm insert positions
      atom <- atom[-insert.inds,]
    }
  }
  het <- matrix(trim(sapply(het.atom, split.string)), byrow=TRUE,
                ncol=nrow(atom.format[ !drop.ind,]), 
                dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )

  output<-list(atom=atom,
               het=het,
               helix=helix,
               sheet=sheet,
               seqres=seqres,
               xyz=as.numeric(t(atom[,c("x","y","z")])),
               calpha = as.logical(atom[,"elety"]=="CA"))

  class(output) <- "pdb"
  return(output)
  
}

