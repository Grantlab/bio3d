## NOTE: 
##   We do not support old-version DSSP any longer
##   Please update your DSSP program to the newest version
"dssp" <-
function (pdb, exepath = "", resno=TRUE) {
    infile <- tempfile()
    outfile <- tempfile()
    write.pdb(pdb, file = infile)
    os1 <- .Platform$OS.type
    if(os1 == "windows") {
       shell(paste(exepath, "dssp ", infile, " ",
                 outfile, sep = ""), ignore.stderr = TRUE)
    } else {
       system(paste(exepath, "dssp ", infile, " ",
                 outfile, sep = ""), ignore.stderr = TRUE)
    }
##
## For Debug (Tue Aug  3 18:22:11 PDT 2010)
##  -- Following multi chain error report from Heiko Strathmann 
##    outfile <- "2jk2.dssp"
##    outfile <- "4q21.dssp"
##
    
    raw.lines <- readLines(outfile)
    unlink(c(infile, outfile))
    type <- substring(raw.lines, 1, 3)
    raw.lines <- raw.lines[-(1:which(type == "  #"))]
    # delete chain breaking lines
    aa <- substring(raw.lines, 14, 14)
    if(any(aa == "!"))
       raw.lines <- raw.lines[-which(aa == "!")]
    cha <- substring(raw.lines, 12, 12)
    sse <- substring(raw.lines, 17, 17)

    # column numbers of phi and psi are different between 
    # the old and new versions of DSSP 
    phi <- as.numeric(substring(raw.lines, 104, 109))
    psi <- as.numeric(substring(raw.lines, 110, 115))

    acc <- as.numeric(substring(raw.lines, 35, 38))


    h.ind <- bounds(which(sse == "H"), pre.sort=FALSE)
    g.ind <- bounds(which(sse == "G"), pre.sort=FALSE)
    e.ind <- bounds(which(sse == "E"), pre.sort=FALSE)
    t.ind <- bounds(which(sse == "T"), pre.sort=FALSE)
    
    h.res <- h.ind;    g.res <- g.ind
    e.res <- e.ind;    t.res <- t.ind
    
    if(resno) {
      res.num  <- as.numeric(substring(raw.lines, 6, 10))
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

    ## ToDo: Add "type" for turns and strands too...
    
    if(length(h.res)>1) {
#      if(is.null(nrow(h.res)))
#        h.s <- as.matrix(t(h.res))
      helix$start  = c(helix$start,h.res[, "start"])
      helix$end    = c(helix$end, h.res[, "end"])
      helix$length = c(helix$length, h.res[, "length"])
      helix$chain  = c(helix$chain, cha[h.ind[, "start"]])
      helix$type   = c(helix$type, sse[h.ind[, "start"]])
    }
    if(length(g.res)>1) {
#      if(is.null(nrow(g.res)))
#        g.s <- as.matrix(t(g.res))
      helix$start  = c(helix$start,g.res[, "start"])
      helix$end    = c(helix$end, g.res[, "end"])
      helix$length = c(helix$length, g.res[, "length"])
      helix$chain  = c(helix$chain, cha[g.ind[, "start"]])
      helix$type   = c(helix$type, sse[g.ind[, "start"]])
    }
    if(length(helix$start) > 0)
       helix <- lapply(helix, function(x) {names(x) <- 1:length(helix$start); return(x)})
    if(length(e.res)>1) {
#      if(is.null(nrow(e.res)))
#        e.s <- as.matrix(t(e.res))
      sheet$start  = c(sheet$start,e.res[, "start"])
      sheet$end    = c(sheet$end, e.res[, "end"])
      sheet$length = c(sheet$length, e.res[, "length"])
      sheet$chain  = c(sheet$chain, cha[e.ind[, "start"]])
    }
    if(length(sheet$start) > 0)
       sheet <- lapply(sheet, function(x) {names(x) <- 1:length(sheet$start); return(x)})
    if(length(t.res)>1) {
#      if(is.null(nrow(t.res)))
#        t.s <- as.matrix(t(t.res))
      turn$start  = c(turn$start,t.res[, "start"])
      turn$end    = c(turn$end, t.res[, "end"])
      turn$length = c(turn$length, t.res[, "length"])
      turn$chain  = c(turn$chain, cha[t.ind[, "start"]])
    }
    if(length(turn$start) > 0)
       turn <- lapply(turn, function(x) {names(x) <- 1:length(turn$start); return(x)})

    out <- list(helix = helix, sheet = sheet,
                turn = turn, phi = phi, psi = psi, acc = acc,
                sse = sse)
}
