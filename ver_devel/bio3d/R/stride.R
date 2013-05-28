"stride" <-
function(pdb, exepath = "", resno=TRUE) {
  
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
  
  # DEBUG: SSE length is inconsistent with the sequence length;
  # Read ASG instead of LOC lines
  sse <- substring(raw.tor, 25,25)
  cha <- substring(raw.tor, 10,10)
  acc <- as.numeric(substring(raw.tor, 65, 69))

  h.ind <- bounds(which(sse == "H"), pre.sort=FALSE)
  g.ind <- bounds(which(sse == "G"), pre.sort=FALSE)
  e.ind <- bounds(which(sse == "E"), pre.sort=FALSE)
  t.ind <- bounds(which(sse == "T"), pre.sort=FALSE)

#  start <- as.numeric(substring(raw.loc, 23,27))
#  end   <- as.numeric(substring(raw.loc, 42,45))
#  chain <- substring(raw.loc, 29,29)
#
#  sse <- substring(raw.loc, 6,9)
#  
#  h.ind <- sse == "Alph"
#  g.ind <- sse == "310H"
#  e.ind <- sse == "Stra"
#  t.ind <- sse == "Turn"
#
#  sse.type <- sse
#  sse.type[h.ind] <- "H"
#  sse.type[g.ind] <- "G"
#  sse.type[e.ind] <- "E"
#  sse.type[t.ind] <- "T"

  h.res <- h.ind;    g.res <- g.ind
  e.res <- e.ind;    t.res <- t.ind
  if(resno) {
    res.num  <- as.numeric(substring(raw.tor, 12, 15))
    if(length(h.ind) > 0) {
      h.res[,"start"] <- res.num[h.ind[,"start"]]
      h.res[,"end"] <- res.num[h.ind[,"end"]]
    }

    if(length(g.ind) > 0) {
      g.res[,"start"] <- res.num[g.ind[,"start"]]
      g.res[,"end"] <- res.num[g.ind[,"end"]]
    }

    if(length(e.ind) > 0) {
      e.res[,"start"] <- res.num[e.ind[,"start"]]
      e.res[,"end"] <- res.num[e.ind[,"end"]]
    }

    if(length(t.ind) > 0) {
      t.res[,"start"] <- res.num[t.ind[,"start"]]
      t.res[,"end"] <- res.num[t.ind[,"end"]]
    }
  }

  sheet = list(start=NULL, end=NULL, length=NULL, chain=NULL)
  helix = list(start=NULL, end=NULL, length=NULL, chain=NULL, type=NULL)
  turn = sheet
  if(length(h.res)>1) {
    helix$start  = c(helix$start,h.res[, "start"])
    helix$end    = c(helix$end, h.res[, "end"])
    helix$length = c(helix$length, h.res[, "length"])
    helix$chain  = c(helix$chain, cha[h.ind[, "start"]])
    helix$type   = c(helix$type, sse[h.ind[, "start"]])
  }
  if(length(g.res)>1) {
    helix$start  = c(helix$start,g.res[, "start"])
    helix$end    = c(helix$end, g.res[, "end"])
    helix$length = c(helix$length, g.res[, "length"])
    helix$chain  = c(helix$chain, cha[g.ind[, "start"]])
    helix$type   = c(helix$type, sse[g.ind[, "start"]])
  }
  if(length(helix$start) > 0)
     helix <- lapply(helix, function(x) {names(x) <- 1:length(helix$start); return(x)})
  if(length(e.res)>1) {
    sheet$start  = c(sheet$start,e.res[, "start"])
    sheet$end    = c(sheet$end, e.res[, "end"])
    sheet$length = c(sheet$length, e.res[, "length"])
    sheet$chain  = c(sheet$chain, cha[e.ind[, "start"]])
  }
  if(length(sheet$start) > 0)
     sheet <- lapply(sheet, function(x) {names(x) <- 1:length(sheet$start); return(x)})
  if(length(t.res)>1) {
    turn$start  = c(turn$start,t.res[, "start"])
    turn$end    = c(turn$end, t.res[, "end"])
    turn$length = c(turn$length, t.res[, "length"])
    turn$chain  = c(turn$chain, cha[t.ind[, "start"]])
  }
  if(length(turn$start) > 0)
     turn <- lapply(turn, function(x) {names(x) <- 1:length(turn$start); return(x)})
  
#  if(any(h.ind | g.ind)) {
#     helix=list(start = c(start[h.ind], start[g.ind]),
#                end    = c(end[h.ind], end[g.ind]),
#                length = ( c(end[h.ind], end[g.ind]) -
#                          c(start[h.ind], start[g.ind]) + 1),
#                chain  = c(chain[h.ind], chain[g.ind]),
#                type   = c(sse.type[h.ind], sse.type[g.ind]))
#  }
#  if(any(e.ind)) {
#     sheet = list(start = start[e.ind],
#                end    = end[e.ind],
#                length = (end[e.ind] - start[e.ind] + 1),
#                chain = chain[e.ind])
#  }
#  if(any(t.ind)) {
#     turn  = list(start = start[t.ind],
#                end    = end[t.ind],
#                length =(end[t.ind] - start[t.ind] + 1),
#                chain = chain[t.ind])
#  }
  out <- list(helix = helix, sheet=sheet, turn=turn,
              phi = phi, psi = psi, acc = acc)
}

