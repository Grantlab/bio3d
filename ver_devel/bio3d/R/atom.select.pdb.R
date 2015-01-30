".is.protein" <- function(pdb) {
  pdb$atom$insert[is.na(pdb$atom$insert)] <- ""
  pdb$atom$chain[is.na(pdb$atom$chain)]   <- ""
  resid <- paste(pdb$atom$chain, pdb$atom$insert, pdb$atom$resno, sep="-")
  
  at.ca <- resid[ pdb$atom$elety == "CA" ]
  at.o  <- resid[ pdb$atom$elety == "O" ]
  at.c  <- resid[ pdb$atom$elety == "C" ]
  at.n  <- resid[ pdb$atom$elety == "N" ]
  
  common <- intersect(intersect(intersect(at.ca, at.o), at.n), at.c)
  return(resid %in% common)
}

".is.nucleic" <- function(pdb) {
  nuc.aa <- c("A",   "U",  "G",  "C",   "T",  "I",
              "DA", "DU", "DG", "DC",  "DT", "DI")
  return(pdb$atom$resid %in% nuc.aa)
}

".is.water" <- function(pdb) {
  hoh <- c("H2O", "OH2", "HOH", "HHO", "OHH", "SOL",
           "WAT", "TIP", "TIP2", "TIP3", "TIP4")
  return(pdb$atom$resid %in% hoh)
}

".is.hydrogen" <- function(pdb) {
  return(substr( gsub("^[123]", "",pdb$atom$elety) , 1, 1) %in% "H")
}

.match.type <- function(pdb, t) {
  if(!is.character(t))
    stop("'type' must be a character vector")
  pdb$atom$type %in% t
}

.match.eleno <- function(pdb, eleno) {
  if(!is.numeric(eleno))
    stop("'eleno' must be a numeric vector")
  pdb$atom$eleno %in% eleno
}

.match.elety <- function(pdb, elety) {
  if(!is.character(elety))
    stop("'elety' must be a character vector")
  pdb$atom$elety %in% elety
}

.match.resid <- function(pdb, resid) {
  if(!is.character(resid))
    stop("'resid' must be a character vector")
  pdb$atom$resid %in% resid
}

.match.chain <- function(pdb, chain) {
  if(!is.character(chain))
    stop("'chain' must be a character vector")
  pdb$atom$chain %in% chain
}

.match.resno <- function(pdb, resno) {
  if(!is.numeric(resno))
    stop("'resno' must be a numeric vector")
  pdb$atom$resno %in% resno
}

.match.segid <- function(pdb, segid) {
  if(!is.character(segid))
    stop("'segid' must be a character vector")
  pdb$atom$segid %in% segid
}

.match.elesy <- function(pdb, elesy) {
  if(!is.character(elesy))
    stop("'elesy' must be a character vector")
  pdb$atom$elesy %in% elesy
}

atom.select.pdb <- function(pdb, string=NULL,
                            type  = NULL, eleno = NULL, elety = NULL,
                            resid = NULL, chain = NULL, resno = NULL,
                            segid = NULL, elesy = NULL, 
                            inverse = FALSE, verbose=FALSE, ...) {

  if(!is.pdb(pdb))
    stop("'pdb' must be an object of class 'pdb'")

  cl <- match.call()
  M <- rep(TRUE, nrow(pdb$atom))

  if(!is.null(string)) {
    M <- switch(string,
                all       =   M,
                protein   =  .is.protein(pdb),
                nucleic   =  .is.nucleic(pdb),
                water     =  .is.water(pdb),
                calpha    =  .is.protein(pdb)  & .match.elety(pdb, "CA"),
                cbeta     =  .is.protein(pdb)  & .match.elety(pdb, "CB"),
                backbone  =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O")),
                back      =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O")),
                ligand    = !.is.protein(pdb)  & !.is.nucleic(pdb) & !.is.water(pdb),
                h         =  .is.hydrogen(pdb),
                noh       = !.is.hydrogen(pdb),
                hetatom   =  .match.type(pdb, "HETATM"),
                atom      =  .match.type(pdb, "ATOM"),
                NA
    )
  }

  if(verbose)
    print(M)

  if(!is.null(type)) {
    M <- M & .match.elety(pdb, type)
  }
  if(!is.null(eleno)) {
    M <- M & .match.elety(pdb, eleno)
  }
  if(!is.null(elety)) {
    M <- M & .match.elety(pdb, elety)
  }
  if(!is.null(resid)) {
    M <- M & .match.resid(pdb, resid)
  }
  if(!is.null(chain)) {
    M <- M & .match.chain(pdb, chain)
  }
  if(!is.null(resno)) {
    M <- M & .match.resno(pdb, resno)
  }
  if(!is.null(segid)) {
    M <- M & .match.segid(pdb, segid)
  }
  if(!is.null(elesy)) {
    M <- M & .match.elesy(pdb, elesy)
  }
  
  if(inverse)
    sele <- as.select(which(!M))
  else
    sele <- as.select(which(M))

  sele$call <- cl
  return(sele)
}
