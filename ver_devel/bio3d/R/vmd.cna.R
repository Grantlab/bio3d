vmd.cna <- function(x, pdb, layout=layout.cna(x, pdb, k=3),
                     col.sphere=NULL, 
                     col.lines="silver",
                     weights=NULL,
                     radius=table(x$communities$membership)/5,
                     alpha=1,
                     vmdfile="network.vmd", pdbfile="network.pdb",
                     full=FALSE, launch=FALSE, exefile=NULL, ...) {

  ## Draw a cna network in VMD

  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }
    
  if(is.null(weights)){
    weights <- igraph::E(x$community.network)$weight
    
    if(is.null(x$call$minus.log)){
      weights <- exp(-weights)
    }
    else{
      if(x$call$minus.log){
        weights <- exp(-weights)
      }
    }
  }
  
  if(is.null(col.sphere)) {
    ## Get colors from network and convert to 0:17 VMD color index
    col.sphere <- match(igraph::V(x$community.network)$color, vmd_colors())-1
  } else {
    ## Check supplied color(s) will work in VMD
    if(!all(col.sphere %in% c(0:17))) {
      warning("Input 'col.sphere' may not work properly in VMD
               - should be 0:17 color index value")
    }
  }


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
  ##radius = V(x$community.network)$size
  ###radius = table(x$raw.communities$membership)/5
  if(!full)  {
      scr <- c(scr, .vmd.sphere( layout, radius=radius, col=col.sphere))

      ## Edges
      ##edge.list <- unlist2(get.adjlist(x$community.network))
      ##start.no <- as.numeric(names(edge.list))
      ##end.no <- as.numeric((edge.list))
      ##inds <- which(end.no > start.no)
      ##start <- layout[start.no[inds],]
      ##end <- layout[end.no[inds],]
      edge.list <- igraph::get.edges(x$community.network, 1:length(igraph::E(x$community.network)))
      start <- layout[edge.list[,1],]
      end <- layout[edge.list[,2],]
      
      ## weights=E(x$community.network)$weight ##/0.2
      scr <- c(scr, .vmd.lines(start=start, end=end,
                               radius=weights, col=col.lines))
  }
      
  ## full network
  if(full) {

      ## Calpha PDB
      if(!all(pdb$atom$elety == "CA")) {
          message("Trimming provided PDB to calpha atoms")
          pdb.ca <- trim(pdb, "calpha")
      }
      else {
          pdb.ca <- pdb
      }

      if(length(x$network[1]) != nrow(pdb.ca$atom))
          stop("Mismatch between provided network and pdb")
      
      ## Edges
      edge.list <- igraph::get.edges(x$network, 1:length(igraph::E(x$network)))
      start <- matrix(pdb.ca$xyz[, atom2xyz(edge.list[,1]) ], ncol=3, byrow=TRUE)
      end <- matrix(pdb.ca$xyz[, atom2xyz(edge.list[,2]) ], ncol=3, byrow=TRUE)

      ## Edge colors and radius
      col2 <- match(igraph::V(x$network)$color, vmd_colors())-1
      names(col2) = 1:nrow(pdb.ca$atom)
      
      col3 = apply(edge.list, 1, function(x) {
          if(col2[x[1]]==col2[x[2]])
              col2[x[1]]
          else
              16 ## black
      })
      
      rad3 = apply(edge.list, 1, function(x) {
          if(col2[x[1]]==col2[x[2]])
              0.1
          else
              0.25
      })
      
      scr2 = .vmd.lines(start=start, end=end, radius=rad3, col=col3)
      scr = c(scr, scr2)

  }
  cat(scr, file=vmdfile, sep="")

  ## Output a PDB file with chain color
  # Use the chain field to store cluster membership data for color in VMD
  ch <- vec2resno(vec=x$communities$membership, resno=pdb$atom[,"resno"])
  write.pdb(pdb, chain=LETTERS[ch], file=pdbfile)

  ## Launch option ...
  ## vmd -pdb network.pdb -e network.vmd
  if(launch) {

    ## Find default path to external program
    if(is.null(exefile)) {
      exefile <- 'vmd'
      if(nchar(Sys.which(exefile)) == 0) {
        os1 <- Sys.info()["sysname"]
        exefile <- switch(os1,
          Windows = 'vmd.exe', # to be updated
          Darwin = '/Applications/VMD\\ 1.9.*app/Contents/MacOS/startup.command',
          'vmd' )
      }
    }
    if(nchar(Sys.which(exefile)) == 0)  
      stop(paste("Launching external program failed\n",
                 "  make sure '", exefile, "' is in your search path", sep=""))
     
    cmd <- paste(exefile, pdbfile, "-e", vmdfile)

    os1 <- .Platform$OS.type
    if (os1 == "windows") {
      cmd <- paste(shQuote(exefile), pdbfile, "-e", vmdfile)
      shell(cmd)
    } else{
      system(cmd)
    }
  }
}
