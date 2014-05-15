# visualize molecular structures

visualize <- function(...)
  UseMethod("visualize")

visualize.xyz <- function(
  xyz, ele.symb = NULL, con = NULL, cell = NULL, type = "l", safety = 1.2,
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rcov", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){

  if(!inherits(xyz, "xyz"))
    stop("'xyz' must be an object of class 'xyz'")
  
  if(length(xyz)==0)
    stop("No coordinates to plot. (length(xyz)==0)")
  
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
#   par.save <- par3d(skipRedraw=TRUE)
#   on.exit(par3d(par.save))

  if(!add){
    open3d()
    par3d(windowRect = windowRect, userMatrix=userMatrix, FOV = FOV, ...)
    bg3d(color=bg.col)
  }
  ids <- rgl.ids()

  if(xyz.axes) ids <- rbind(ids, addXYZ(lwd = lwd.xyz, cex = cex.xyz))
  if(abc.axes) ids <- rbind(ids, addABC(cell, lwd = lwd.abc, cex = cex.abc))
  if(pbc.box ) ids <- rbind(ids, addPBCBox(cell, lwd = lwd.pbc.box))
  
  if(nchar(type)>1) type <- strsplit(type, "")[[1]]
  if(!all(type %in% c("l","s","p")))
    stop("Unrecognized 'type'")
  
  if(centre) {
    cent <- centres(xyz)
    xyz[seq(1, length(xyz), 3)] <- xyz[seq(1, length(xyz), 3)] - cent[1]
    xyz[seq(2, length(xyz), 3)] <- xyz[seq(2, length(xyz), 3)] - cent[2]
    xyz[seq(3, length(xyz), 3)] <- xyz[seq(3, length(xyz), 3)] - cent[3]
  }

  if("l" %in% type) {
    if(is.null(con)) {
      warning("Unspecifyed connectivity: Computing connectivity from coordinates...")
      con <- connectivity.xyz(x = xyz, ele.symb = ele.symb, safety = safety, by.block = TRUE)
    }
    if(!is.null(con)){
      ind <- t(con)
      seg.id <- segments3d(
      xyz[seq(1,length(xyz),3)][ind],
      xyz[seq(2,length(xyz),3)][ind],
      xyz[seq(3,length(xyz),3)][ind],
      color = col[ind], lwd=lwd, ...)

      seg.id <- data.frame(id = seg.id, type = "atom.seg")
      ids <- rbind(ids, seg.id)
    } else{
      warning("'con' is empty: 'type' has been set to 'p'.")
      type <- "p"
    }
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
  pdb, elety.custom = NULL, atom.sel = atom.select(pdb, "notwater"),
  cell = NULL, type = "l", safety = 1.2,
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rcov", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){

  if(!is.pdb(pdb)) stop("'pdb' must be an object of class pdb. See read.pdb")

  pdb <- trim.pdb(pdb, atom.sel)  
  ele.symb <- atom2ele(pdb$atom[,"elety"], elety.custom)
  con <- connectivity.pdb(pdb)
  
  visualize.xyz(
    pdb$xyz, ele.symb = ele.symb, con, cell, type, safety,
    xyz.axes, abc.axes, pbc.box, centre, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
    cex.xyz, cex.abc, col, radii, bg.col, add, windowRect,
    userMatrix, FOV, ...)
}

# visualize.cna <- function(cna, pdb, safety = 2.7, ...){
#   if(!"cna" %in% class(cna))
#     stop("'cna' must an object of class 'cna'")
#   if(missing(pdb))
#     stop("Please specify a 'pdb' object")
#   if(!is.pdb(pdb))
#     stop("'pdb' must an object of class 'pdb'")
#   
#   ca.pdb <- trim.pdb(pdb, atom.select(pdb, "calpha", verbose = FALSE))
#   ca.con <- connectivity(ca.pdb,safety=safety)
#   net.vertex.color <- V(cna$network)$color
#   visualize(ca.pdb, con = ca.con, col = net.vertex.color)
# 
#   com.net.vertex.color <- V(cna$community.network)$color
#   radii <- V(net$community.network)$size
#   radii <- 5*radii/max(radii)
#   com.net.weight <- E(net$community.network)$weight
#   membership.centres <- centres(ca.pdb, factor = cna$communities$membership)
#   membership.centres <- matrix(membership.centres, ncol=3, byrow=TRUE)
#   spheres3d(membership.centres, col = com.net.vertex.color, radius = radii, alpha = 0.5)
#   
#   com.net.con <- apply(get.edgelist(net$community.network), 2, as.integer)
#   cyls <- apply(com.net.con, 1,
#     function(ids){
#       cyl <- cylinder3d(membership.centres[ids,], radius=com.net.weight, sides=20)
#       return(cyl)
#     })
#   cyls.ids <- lapply(cyls, shade3d, col = "white")
# }
