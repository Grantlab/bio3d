# visualize molecular structures

visualize <- function(...)
  UseMethod("visualize")

visualize.default <- function(
  xyz, ele.symb = NULL, con = NULL, cell = NULL, type = "p",
  xyz.axes = TRUE, abc.axes = FALSE, pbc.box = FALSE, 
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rvdw", bg.col = "#FAFAD2",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){
  
  oops <- require(rgl)
  if(!oops)
    stop("Please install the rgl package from CRAN")
  
  if(is.null(ele.symb))
    ele.symb <- rep("Xx", length(xyz)/3)

  M <- match(ele.symb,elements[,"symb"])
  M[is.na(M)] <- 1 # Unrecognized elements are considered as dummy atoms    

  if(is.null(col))
    col <- do.call(rgb, elements[M, c("red","green","blue")])
  if(length(col) != length(xyz)/3)
    warning("'col' has been recycled")

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

  if(type == "l") {
    if(is.null(con)) {
      type == "p"
      warning("Undefined connectivity: 'type' has been set to 'p'.")
    }
    else {
      ind <- t(con)
      seg.id <- segments3d(
        xyz[seq(1,length(xyz),3)][ind],
        xyz[seq(2,length(xyz),3)][ind],
        xyz[seq(3,length(xyz),3)][ind],
        color = col[ind], lwd=lwd, ...)
    }
    seg.id <- data.frame(id = seg.id, type = "atom.seg")
    ids <- rbind(ids, seg.id)
  }
  else if(type == "s") {
    if(is.character(radii[1])){
      if(! radii[1] %in% c("rcov", "rvdw") )
        stop("'radii' must be one of 'rcov', 'rvdw' or a numerical vector")
      radii <- elements[M,radii[1]]
    }
    sph.id <- spheres3d(
      xyz[seq(1,length(xyz),3)],
      xyz[seq(2,length(xyz),3)],
      xyz[seq(3,length(xyz),3)],
      color = col, radius = radii, ...)
    sph.id <- data.frame(id = sph.id, type= "atom.sph")
    ids <- rbind(ids, sph.id)
  }
  else if(type == "p") {
    pts.id <- points3d(
      xyz[seq(1,length(xyz),3)],
      xyz[seq(2,length(xyz),3)],
      xyz[seq(3,length(xyz),3)],
      color = col, ...)
    pts.id <- data.frame(id = pts.id, type = "atom.pts")
    ids <- rbind(ids, pts.id)
  }
  else stop("Unrecognized 'type'")
  invisible(ids)
}

visualize.pdb <- function(
  x, elety.custom = NULL, con = NULL, cell = NULL, type = "s",
  xyz.axes = TRUE, abc.axes = FALSE, pbc.box = FALSE, 
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rvdw", bg.col = "#FAFAD2",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){

  if(!is.pdb(pdb)) stop("'pdb' must be an object of class pdb. See read.pdb")
  
  ele.symb <- atom2ele(x$atom[,"elety"], elety.custom)
  con1 <- match(con[,1], x$atom[,"eleno"])
  con2 <- match(con[,2], x$atom[,"eleno"])
  con <- data.frame(eleid.1 = con1, eleid.2 = con2)
  
  visualize.default(
    x$xyz, ele.symb = ele.symb, con, cell, type,
    xyz.axes, abc.axes, pbc.box, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
    cex.xyz, cex.abc, col, radii, bg.col, add, windowRect,
    userMatrix, FOV, ...)
}
