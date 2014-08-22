## This script intends to do post-processing of results from WISP, 
## a software for suboptimal paths calculation in a correlation network.
## It calculates comparative length distribution of paths and node degeneracy
## for multiple networks.
##
## Arguments: 
##   inp           path to WISP results
##   state         names for multiple inputs, used to group and name outputs
##   col           colors for outputs of different states
##   pdb           reference PDB
##   outprefix     prefix of output file names
##   concise       if TRUE only show nodes on paths
##   cutoff        cutoff of node degeneracy to show in outputs
##   normalize     if TRUE node degeneracy is calculated as percentage of total number of paths 

subopt.path <- function(inp, state = 1:length(inp), col=1:length(inp), pdb = NULL, outprefix="", 
   concise = TRUE, cutoff = 0.1, normalize = TRUE) {
   require(bio3d)
   require(ggplot2)

   out <- NULL

   # read sub-optimal path lengths
   data <- lapply(inp, function(x) readLines(file.path(x, "simply_formatted_paths.txt")))
   y <- lapply(data, function(x) as.numeric(unlist(lapply(strsplit(x, " "), "[", 1))) )
  
   factor <- rep(state, unlist(lapply(y, length)))
   # plot path length distribution
   gcol = col; names(gcol) = state
   vy = unlist(y)
   bw=(max(vy) - min(vy)) * 0.02
   wd=bw/1.5
   df <- data.frame(State=as.factor(factor), Length=unlist(y))
   df <- cbind(df, bw=bw)
   pdf(file=paste(outprefix, "length_distrib.pdf", sep=""), width=3.5, height=3, pointsize=8)
   print(ggplot(df, aes(x=Length, group=State, color=NULL, fill=State, bw=bw)) + geom_histogram(aes(y=(..density..)*bw), alpha=0.8, position=position_dodge(width=wd), binwidth=bw)+scale_fill_manual(values=gcol)+ylab("Probability")+theme_bw() )
   dev.off()
   out$length <- y
   names(out$length) <- state 
   
   # read node numbers on paths
   y <- lapply(data, function(x) as.numeric(unlist(lapply(strsplit(x, " "), "[", -1))) )
   
   # correct node # (WISP return 0-based #)
   y <- lapply(y, "+", 1)
#   y[state != "GTP"] <- lapply(y[state != "GTP"], "+", 1)
   if(!is.null(pdb)) {
      resno <- as.numeric(pdb$atom[pdb$calpha, "resno"])
      y <- lapply(y, function(x) resno[x])
   }
   y0 <- y
 
   if(concise) { 
      # re-number node to get more concise plot
      ii <- sort(unique(unlist(y)))
      y <- lapply(y, match, ii)
   } 
   
   factor <- rep(state, unlist(lapply(y, length)))

   # plot node degeneracy (un-normalized)
   bw=0.9
   wd=bw/1.5
   df <- data.frame(State=as.factor(factor), Node=unlist(y))
   df <- cbind(df, bw=bw)
   pdf(file=paste(outprefix, "degeneracy_distrib.pdf", sep=""), width=5, height=3, pointsize=8)
   print(ggplot(df, aes(x=Node, group=State, color=NULL, fill=State, bw=bw)) + geom_histogram(aes(y=..count..), alpha=0.8, position=position_dodge(width=wd), binwidth=bw)+scale_fill_manual(values=gcol)+ylab("Degeneracy")+theme_bw() )
   dev.off()
#   out$path <- y0
#   names(out$path) <- state 
  
   # store node degeneracy 
   yp <- lapply(y0, table)
    
   if(normalize) {
      yp <- lapply(yp, function(x) x/max(x))
   }

   yy <- lapply(yp, function(x) x[x >= cutoff])
   nyy <- unique(names(unlist(yy)))
   i <- as.numeric(nyy)
   nyy <- nyy[order(i)]

   o <- lapply(yp, function(x) {
        x <- x[match(nyy, names(x))]
       x[is.na(x)] <- 0
       names(x) <- nyy 
       x } )
   if(!is.null(pdb)) {
      resid <- aa321(pdb$atom[pdb$calpha, "resid"])
      o <- lapply(o, function(x) {
            n <- paste(resid[match(names(x), resno)], names(x), sep="")
            names(x) <- n
            x } )
   }
   names(o) <- state 
   out$degeneracy <- do.call(rbind, o)
   if(normalize) {
#      out$degeneracy <- out$degeneracy / max(out$degeneracy) 
      out$degeneracy <- round(out$degeneracy, digits=2)
   }
   write.csv(out$degeneracy, file = paste(outprefix, "degeneracy.csv", sep=""))
   out$call <- match.call()
   class(out) <- "subopt"

   return(out)
}

print.subopt <- function(x) {
   if(!inherits(x, "subopt")) stop("Not subopt class")
   if(is.null(x$call$normalize) || is.true(x$call$normalize)) 
      print(format(x$degeneracy, nsmall=2), quote=FALSE)
   else 
      print(x$degeneracy)
}
