######## This script calculates length distribution, node degeneracy, 
######## and identifies potentially key residues based on 
######## sub-optimal paths analysis of dynamical networks

subopt.path <- function(path, state = 1:length(path), col=1:length(path), pdb = NULL, outprefix="", 
   rm.gaps = TRUE, cutoff.npath = 50, normalize = TRUE) {
   require(bio3d)
   require(ggplot2)

   out <- NULL

   # read sub-optimal path lengths
   data <- lapply(path, function(x) readLines(file.path(x, "simply_formatted_paths.txt")))
   y <- lapply(data, function(x) as.numeric(unlist(lapply(strsplit(x, " "), "[", 1))) )
  
   factor <- rep(state, unlist(lapply(y, length)))
   # plot path length distribution
   gcol = col; names(gcol) = state
   vy = unlist(y)
   bw=(max(vy) - min(vy)) * 0.04
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
 
   if(rm.gaps) { 
      # re-number node to get more concise plot
      ii <- sort(unique(unlist(y)))
      y <- lapply(y, match, ii)
   } 
   
   factor <- rep(state, unlist(lapply(y, length)))
   # plot node degeneracy
   bw=0.9
   wd=bw/1.5
   df <- data.frame(State=as.factor(factor), Node=unlist(y))
   df <- cbind(df, bw=bw)
   pdf(file=paste(outprefix, "degeneracy_distrib.pdf", sep=""), width=5, height=3, pointsize=8)
   print(ggplot(df, aes(x=Node, group=State, color=NULL, fill=State, bw=bw)) + geom_histogram(aes(y=..count..), alpha=0.8, position=position_dodge(width=wd), binwidth=bw)+scale_fill_manual(values=gcol)+ylab("Degeneracy")+theme_bw() )
   dev.off()
#   out$path <- y0
#   names(out$path) <- state 
  
   # identify potentially key residues 
   yp <- lapply(y0, table)
   yy <- lapply(yp, function(x) x[x > cutoff.npath])
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
      out$degeneracy <- out$degeneracy / max(out$degeneracy) 
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
