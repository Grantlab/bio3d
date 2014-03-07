# visualize molecular structures

visualize <- function(...)
  UseMethod("visualize")

visualize.default <- function(
  xyz, ele.symb = NULL, con = NULL, cell = NULL, type = "l",
  xyz.axes = TRUE, abc.axes = FALSE, pbc.box = FALSE, 
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rcov", bg.col = "#FAFAD2",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){
  
  oops <- require(rgl)
  if(!oops)
    stop("Please install the rgl package from CRAN")
  
  if(is.null(ele.symb))
    ele.symb <- rep("Xx", length(xyz)/3)
  
  if(length(ele.symb) != length(xyz)/3)
    stop("'xyz' and 'ele.symb' must have matching lengths")
  if(any(is.na(match(ele.symb, elements$symb))))
    ele.symb <- atom2ele(ele.symb)
  ele.symb[is.na(ele.symb)] <- "Xx"

  M <- match(ele.symb,elements[,"symb"])
  M[is.na(M)] <- 1 # Unrecognized elements are considered as dummy atoms    

  if(is.null(col))
    col <- do.call(rgb, elements[M, c("red","green","blue")])
  if(length(col) != length(xyz)/3){
    if(length(col)!=1) warning("'col' has been recycled")
    col <- rep(col, length = length(xyz)/3)
  }

  if(is.null(cell)) {
    if(abc.axes) {
      abc.axes <- FALSE
      warning("Cell parameters are not specified: 'abc.axes' has been set to FALSE")
    }
    if(pbc.box) {
      pdb.box <- FALSE
      warning("Cell parameters are not specified: 'pdb.box' has been set to FALSE")
    }
  }

  if(!add){
    open3d()
    par3d(windowRect = windowRect, userMatrix=userMatrix, FOV = FOV, ...)
    bg3d(color=bg.col)
  }
  ids <- rgl.ids()

  if(xyz.axes) ids <- rbind(ids, addXYZ(lwd = lwd.xyz, cex = cex.xyz))
  if(abc.axes) ids <- rbind(ids, addABC(cell, lwd = lwd.abc, cex = cex.abc))
  if(pbc.box ) ids <- rbind(ids, addPBCBox(cell, lwd = lwd.pbc.box))

  if(is.null(con)) {
    warning("Unspecifyed connectivity: Computing connectivity from coordinates...")
    con <- connectivity(x = xyz, ele.symb = ele.symb, by.block = TRUE)
  }
  
  if(is.null(con)){
    warning("'con' is empty: 'type' has been set to 'p'.")
    type <- "p"
  }
  
  if(nchar(type)>1) type <- strsplit(type, "")[[1]]
  if(!all(type %in% c("l","s","p")))
    stop("Unrecognized 'type'")

  if("l" %in% type) {
    ind <- t(con)
    seg.id <- segments3d(
    xyz[seq(1,length(xyz),3)][ind],
    xyz[seq(2,length(xyz),3)][ind],
    xyz[seq(3,length(xyz),3)][ind],
    color = col[ind], lwd=lwd, ...)

    seg.id <- data.frame(id = seg.id, type = "atom.seg")
    ids <- rbind(ids, seg.id)
  }
  if("s" %in% type) {
    if(is.character(radii[1])){
      if(! radii[1] %in% c("rcov", "rvdw") )
        stop("'radii' must be one of 'rcov', 'rvdw' or a numerical vector")
      radii <- elements[M,radii[1]]*
        ifelse(radii[1] == "rcov" & "l" %in% type, 0.5, 1)
    }
    sph.id <- spheres3d(
      xyz[seq(1,length(xyz),3)],
      xyz[seq(2,length(xyz),3)],
      xyz[seq(3,length(xyz),3)],
      color = col, radius = radii, ...)
    sph.id <- data.frame(id = sph.id, type= "atom.sph")
    ids <- rbind(ids, sph.id)
  }
  if("p" %in% type) {
    pts.id <- points3d(
      xyz[seq(1,length(xyz),3)],
      xyz[seq(2,length(xyz),3)],
      xyz[seq(3,length(xyz),3)],
      color = col, ...)
    pts.id <- data.frame(id = pts.id, type = "atom.pts")
    ids <- rbind(ids, pts.id)
  }
  invisible(ids)
}

visualize.pdb <- function(
  x, elety.custom = NULL, con = NULL, cell = NULL, type = "l",
  xyz.axes = TRUE, abc.axes = FALSE, pbc.box = FALSE, 
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rcov", bg.col = "#FAFAD2",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){

  if(!is.pdb(x)) stop("'pdb' must be an object of class pdb. See read.pdb")

  ele.symb <- atom2ele(x$atom[,"elety"], elety.custom)

  visualize.default(
    x$xyz, ele.symb = ele.symb, con, cell, type,
    xyz.axes, abc.axes, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
    cex.xyz, cex.abc, col, radii, bg.col, add, windowRect,
    userMatrix, FOV, ...)
}

# visualize.cna <- function(x, layout, weights=NULL, vertex.size = NULL, col = NULL,
#  ...){
# if(is.pdb(layout)) layout <- layout.cna(layout, layout, k=3), 
# 
#   plot3d.cna <- function(x,
# pdb = NULL,
# weights=NULL,
# vertex.size = NULL,
# layout = layout.cna(x, pdb, k=3),
# col = NULL,
# ...){
#   ##
#   ## Plot a cna network graph in 3D with RGL
#   ##    plot3d.cna(net, pdb)
#   ##  Or just
#   ##    rglplot(net$community.network)
#   ##
#   
#   ## Check if x has 'cna' class
#   if(!"cna" %in% class(x)){
#     stop("Input 'x' object must be a 'cna' class object")
#   }
#   
#   ##  oops <- require(igraph)
#   ##  if (!oops) {
#   ##    warning("igraph package missing: Please install, see: ?install.packages")
#   ##  }
#   
#   if(is.null(weights)){
#     weights <- E(x$community.network)$weight
#     
#     if(is.null(x$call$minus.log)){
#       weights <- exp(-weights)
#     }
#     else{
#       if(x$call$minus.log){
#         weights <- exp(-weights)
#       }
#     }
#     weights <- weights*10
#   }
#   
#   ## Obtain the plot coords...
#   if(!is.null(pdb) && is.null(layout)) {
#     cat("Obtaning layout from PDB structure\n")
#     layout = layout.cna(x, pdb, k=3)
#   }
#   if(is.null(pdb) && is.null(layout)) {
#     cat("Obtaning guestimated layout with fruchterman.reingold\n")
#     layout <- layout.fruchterman.reingold(x$community.network, weights=weights)
#   }
#   if(dim(layout)[2] != 3){
#     stop("Input 'layout' must be an Nx3 matrix, where N is the number of communities")
#   }
#   
#   rglplot(x$community.network,
#           edge.width = weights,
#           layout = layout,
#           vertex.size = vertex.size,
#           vertex.color <- col)
# }
# }
