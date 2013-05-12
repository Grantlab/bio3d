"stride" <-
function(pdb, exepath = "") {
  
  infile  <- tempfile()
  outfile <- tempfile()
  write.pdb(pdb, file=infile)
  os1 <- .Platform$OS.type
  if(os1 == "windows") {
     shell( paste(exepath,"stride -f",outfile," ",infile,sep="") )
  } else {
     system( paste(exepath,"stride -f",outfile," ",infile,sep="") )
  }
  raw.lines <- readLines(outfile)
  type <- substring(raw.lines, 1, 3)
  unlink(c(infile, outfile))

  raw.loc <- raw.lines[type == "LOC"]
  raw.tor <- raw.lines[type == "ASG"]

  phi <- as.numeric(substring(raw.tor, 43,49))
  psi <- as.numeric(substring(raw.tor, 53,59))

  start <- as.numeric(substring(raw.loc, 23,27))
  end   <- as.numeric(substring(raw.loc, 42,45))
  chain <- substring(raw.loc, 29,29)

  sse <- substring(raw.loc, 6,9)
  
  h.ind <- sse == "Alph"
  g.ind <- sse == "310H"
  e.ind <- sse == "Stra"
  t.ind <- sse == "Turn"

  sse.type <- sse
  sse.type[h.ind] <- "H"
  sse.type[g.ind] <- "G"
  sse.type[e.ind] <- "E"
  sse.type[t.ind] <- "T"

  sheet = list(start=NULL, end=NULL, length=NULL, chain=NULL)
  helix = list(start=NULL, end=NULL, length=NULL, chain=NULL, type=NULL)
  turn = sheet
  
  if(any(h.ind | g.ind)) {
     helix=list(start = c(start[h.ind], start[g.ind]),
                end    = c(end[h.ind], end[g.ind]),
                length = ( c(end[h.ind], end[g.ind]) -
                          c(start[h.ind], start[g.ind]) + 1),
                chain  = c(chain[h.ind], chain[g.ind]),
                type   = c(sse.type[h.ind], sse.type[g.ind]))
  }
  if(any(e.ind)) {
     sheet = list(start = start[e.ind],
                end    = end[e.ind],
                length = (end[e.ind] - start[e.ind] + 1),
                chain = chain[e.ind])
  }
  if(any(t.ind)) {
     turn  = list(start = start[t.ind],
                end    = end[t.ind],
                length =(end[t.ind] - start[t.ind] + 1),
                chain = chain[t.ind])
  }
  out <- list(helix = helix, sheet=sheet, turn=turn,
              phi = phi, psi = psi)
}

