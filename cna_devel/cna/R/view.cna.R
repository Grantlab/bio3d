view.cna <- function(x, layout=layout.pdb(pdb, x),
                     col.sphere=0:32, col.lines="silver",
                     weights=E(x$clustered.network)$weight,
                     radius=table(x$raw.communities$membership)/5,
                     alpha=1,
                     vmdfile="network.vmd", pdbfile="network.pdb",
                     launch=FALSE) {

  ## Draw a cna network in VMD
  
  ##-- VMD draw functions for sphere, lines and cone
  .vmd.sphere <- function(cent, radius=5, col="red", resolution=25) {
    ## .vmd.sphere( matrix(c(0,0,0, 1,1,1), ncol=3,byrow=T) )
    
    if(ncol(cent) != 3)
      stop("Input 'cent' should be a 3 col xyz matrix")

    n <- nrow(cent); scr <- rep(NA, n) 
    if(length(col) != n)
      
      col <- rep(col, length=n)
    
    if(length(radius) != n)
      radius <- rep(radius, length=n)
    
    for(i in 1:n) {
      scr[i] <- paste0("draw color ", col[i],
                       "\ndraw sphere {",
                       paste(cent[i,], collapse = " "),
                       "} radius ", radius[i],
                       " resolution ",resolution, "\n")
    }
    return(scr)
  }
  
  
  .vmd.lines <- function(start, end, radius=0.2, col="silver", resolution=25) {
    
    ## .vmd.lines( start=matrix(c(0,0,0), ncol=3,byrow=T),
    ##         end=matrix(c(1,1,1), ncol=3,byrow=T) )
    
    if(ncol(start) != 3)
      stop("Input 'start' and 'end' should be 3 col xyz matrices")
    n <- nrow(start); scr <- rep(NA, n)
    
    if(length(col) != n)
      col <- rep(col, length=n)
    
    if(length(radius) != n)
      radius <- rep(radius, length=n)
    
    for(i in 1:n) {
      scr[i] <- paste0("draw color ", col[i],
                       "\ndraw cylinder {",
                       paste(start[i,], collapse = " "),
                       "} {", paste(end[i,], collapse = " "),
                       "} radius ", radius[i],
                       " resolution ",resolution, "\n")
    }
    return(scr)
  }
  
  .vmd.cone <- function(start, end, radius=5, col="silver", resolution=25) {
    warning("not here yet")
  }

  ##- Set alpha if needed
  scr <- NULL
  if(alpha != 1)
    scr <- paste("material change opacity Transparent",
                 alpha,"\ndraw material Transparent\n")
  
  ##- Lets get drawing
  ##radius = V(x$clustered.network)$size
  ###radius = table(x$raw.communities$membership)/5
  scr <- c(scr, .vmd.sphere( layout, radius=radius, col=col.sphere))

  ## Edges
  edge.list <- unlist2(get.adjlist(x$clustered.network))
  start <- layout[as.numeric(names(edge.list)),]
  end <- layout[as.numeric((edge.list)),]
  ###weights=E(x$clustered.network)$weight ##/0.2
  scr <- c(scr, .vmd.lines( start=start, end=end,
                           radius=weights, col=col.lines))

  cat(scr, file=vmdfile, sep="")

  ## Output a PDB file with chain color
  # Use the chain field to store cluster membership data for color in VMD
  ch <- vec2resno(vec=x$raw.communities$membership, resno=pdb$atom[,"resno"])
  write.pdb(pdb, chain=LETTERS[ch], file=pdbfile)

  ## Launch option ...
  ## vmd -pdb network.pdb -e network.vmd
  if(launch) {
    cmd <- paste("vmd", pdbfile, "-e", vmdfile)

    os1 <- .Platform$OS.type
    if (os1 == "windows") {
      shell(shQuote(cmd))
    } else{
      if(Sys.info()["sysname"]=="Darwin") {
        system(paste("/Applications/VMD\\ 1.9.1.app/Contents/MacOS/startup.command",pdbfile, "-e", vmdfile))
      }
      else{
        system(cmd)
      }
    }
  }
}
