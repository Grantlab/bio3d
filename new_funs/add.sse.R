add.sse <- function(x, type = c("classic", "fancy"), helix.col = "gray20",
                    sheet.col = "gray80", ...) {

   if(!inherits(x, "sse"))
      stop("Input 'x' should be an object of class 'sse'")

   type = match.arg(type)
   sse <- x

   ##-- SSE annotation 
   if(type == 'classic') {
      if(is.pdb(sse)) {
         ss <- pdb2sse(sse, verbose=FALSE)
      } else {
         ss <- sse$sse
      }

      ## Add segment id to SSE element (if availabe)
      ## to take care of missing residues, multi-chain, etc
      id <- sub("^[^_]*_[^_]*_[^_]*_", "", names(ss))
      string <- paste(ss, id)
      string[!ss %in% c('H', 'E', 'G')] <- " "

      ss2no <- c(H=1, G=1, E=101, other=201)
      rl1 <- rle2(string)
      rl1$values <- ss2no[ss[rl1$inds]]
      rl1$values[is.na(rl1$values)] <- ss2no['other']

      ## check if segments of the same type were put together
      rl2 <- rle(rl1$values)
      if(any(rl2$lengths >= 2)) {
         rl2$values <- seq_along(rl2$values)
         tl <- split(rl1$values, f=as.factor(inverse.rle(rl2)))
         tl <- lapply(tl, function(x) 
            if(length(x)<=1) x 
            else x + seq_along(x) - 1
         )
         rl1$values <- as.vector(unlist(tl))
      }
   
      segments <- inverse.rle(rl1)
      stat <- sort(unique(segments))
      nstat <- seq_along(stat)
      names(nstat) <- stat
      segments <- as.vector(nstat[as.character(segments)])

      segment.col <- rep(NA, length(stat))
      segment.col[stat < 100] <- helix.col
      segment.col[stat > 100 & stat < 200] <- sheet.col
      
      dots <- list(...)
      dots$x <- segments
      dots$segment.col <- segment.col
    
      do.call(add.segment, dots)
   }
}

##- Plot colored segments in margins for e.g. sse annotation
add.segment <- function(x, pos = c("bottom", "left", "top", "right"),
                        segment.col = vmd.colors(), segment.min = 1, pad = 0.02,
                        p = NULL) { 

  if(!requireNamespace("lattice", quietly=TRUE))
     stop("Please install lattice package from CRAN")

   if(!is.numeric(x)) 
      stop("Input 'x' should be an integer vector") 
   
   margin.segments <- x 
   pos <- match.arg(pos, several.ok = TRUE)
  
   ##-- get graphic object from input or current device
   if(is.null(p))
      p <- lattice::trellis.last.object()

   ##-- check graphic object
   if(!inherits(p, "trellis"))
      stop("Input 'p' should be NULL or an object of class 'trellis'")

   xlim <- p$x.limits
   ylim <- p$y.limits
   uni <- 1/(max(xlim)-min(xlim))
   padref <- pad/uni

   ##-- Adjust margins for 'sse'/'segment'
   if("top" %in% pos) 
      ylim[2] = ylim[2] + padref
   if("bottom" %in% pos) 
      ylim[1] = ylim[1] - padref
   if("left" %in% pos) 
      xlim[1] = xlim[1] - padref
   if("right" %in% pos) 
      xlim[2] = xlim[2] + padref

   ##-- Format margin annotation object
   grps <- table(margin.segments)
   ##- Exclude small grps less than 'segment.min'
   grps = names( grps[grps > segment.min] ) 
 
   store.grps <- NULL
   for(i in 1:length(grps)) {
     store.grps <- rbind(store.grps,
       cbind( bounds(which(margin.segments == grps[i])),
             "grp"=as.numeric(grps[i])) )
   }

   ##- Margin segment colors (don't draw if color=NA)
   if(is.null(segment.col)) {
     segment.col <- store.grps[,"grp"]
   } else {
     segment.col <- segment.col[(store.grps[,"grp"])]
   }

   ##- Remove segments with color=NA
   store.grps <- store.grps[!is.na(segment.col), , drop=FALSE]
   segment.col <- segment.col[!is.na(segment.col)]

   ##-- Update trellis plot with new annotations 
   curr.args <- p$panel.args.common$draw.segment.args
   new.args <- c(curr.args, list(list(store.grps=store.grps,       
               limits = c(xlim, ylim), pad=padref, pos=pos,
               segment.col=segment.col)))
   names(new.args) <- seq_along(new.args)  ## named nested argument lists

   p <- update(p, draw.segment.args = new.args, xlim = xlim, ylim = ylim)

   inner.box <- p$panel.args.common$inner.box
   if(!is.null(inner.box) && !inner.box$draw) {
      inner.box$draw = TRUE
      p <- update(p, inner.box = inner.box)
   }
 
   ##-- Return plot
   p
}

.draw.segment <- function(start, length, limits, pad, fill.col="gray", side=1) {
  ##-- Draw Annotation On Plot Margins, used for SSE and CLUSTER members
  ##    draw.segment(store.grps[,"start"], store.grps[,"length"],
  ##                 xymin=xymin, xymax=xymax, side=1, fill.col="red")

  switch(side,
    '1' = {
    ## Bottom Margin
    lattice::lrect(x=start-0.5,
          y=limits[3], 
          col=fill.col, 
          just=c("left","bottom"),
          width=length-0.5,
          height=pad, 
          border = NA)
    },
    '2' = {
    ## Left Margin
    lattice::lrect(x=limits[1],
          y=start-0.5, 
          col=fill.col,
          just=c("left","bottom"),
          width=pad,
          height=length-0.5,
          border = NA)
    },
    '3' = {
    ## Top Margin
    lattice::lrect(x=start-0.5,
          y=limits[4] - pad,
          col=fill.col,
          just=c("left","bottom"),
          width=length-0.5,
          height=pad,
          border = NA)
    },
    '4' = {
    ## Right Margin
    lattice::lrect(x=limits[2] - pad,
          y=start-0.5,
          col=fill.col,
          just=c("left","bottom"),
          width=pad,
          height=length-0.5,
          border = NA)
    } )
}

## New panel function for drawing segment (e.g. sse) 
## and outlines surrounding plot regions
.new.lattice.panel <- function(...) {
   dots <- list(...)
   args <- dots$draw.segment.args
   lattice::panel.contourplot(...)

   ##- draw segments
   for(j in 1:length(args)) {
      for(i in match(args[[j]]$pos, c("bottom", "left", "top", "right"))) {
         .draw.segment(start = args[[j]]$store.grps[,"start"],
            length = args[[j]]$store.grps[,"length"],
            limits = args[[j]]$limits,
            pad = args[[j]]$pad,
            fill.col = args[[j]]$segment.col, side=i)
      }
   }

   ##- draw outlines of plot regions
   if(!is.null(dots$inner.box) && dots$inner.box$draw) {
      limits = dots$init.limits
      vp <- viewport(x=limits[1], y=limits[3], 
               width=diff(limits[1:2]), height=diff(limits[3:4]), 
               just=c("left","bottom"), default.units="native", 
               gp=gpar(col="black", fill=NA, lty=1, lwd=1), clip="off")
      if(dots$inner.box$type == "full") {

         grid.rect(vp=vp)

      } else if(dots$inner.box$type == "lower") {
         grid.lines(y=0, vp=vp)
         grid.lines(x=1, vp=vp)
         grid.lines(vp=vp)
      } else if(dots$inner.box$type == "upper") {
         grid.lines(y=1, vp=vp)
         grid.lines(x=0, vp=vp)
         grid.lines(vp=vp)
      }
   }
   
   ##- draw diagonal axis
   diag.scales = dots$diag.scales
   if(!is.null(diag.scales)) {
      vp <- viewport(x=limits[1], y=limits[3],
         width=diff(limits[1:2]), height=diff(limits[3:4]),
         just=c("left","bottom"), default.units="native",
         gp=gpar(col="black"), clip="off", 
         xscale=limits[1:2], yscale=limits[3:4])
      if(dots$inner.box$type == "lower") {
         grid.text(label=diag.scales$labels, x=diag.scales$at, y=diag.scales$at, 
                   just=c("right", "bottom"), vp=vp, default.units="native")
         grid.text(label=diag.scales$diag.label, x=0.45, y=0.55, 
                   just="center", rot=45, vp=vp)
      } else if(dots$inner.box$type == "upper") {
         grid.text(label=diag.scales$labels, x=diag.scales$at, y=diag.scales$at, 
                   just=c("left", "top"), vp=vp, default.units="native")
         grid.text(label=diag.scales$diag.label, x=0.55, y=0.45, 
                   just="center", rot=45, vp=vp)
      }
   }
}

