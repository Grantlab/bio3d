# visualize molecular structures

visualize <- function(...)
  UseMethod("visualize")

visualize.xyz <- function(
  xyz, elesy = NULL, con = NULL, cell = NULL, type = "l",
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rvdw", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...) {

  if(!inherits(xyz, "xyz"))
    stop("'xyz' must be an object of class 'xyz'")
  
  if(nrow(xyz) != 1)
    stop("'x' must be a single row 'xyz' matrix")
    
  if(is.null(elesy)) {
    warning("'elesy' is not defined. All atoms have been considered as dummy atoms")
    elesy <- rep("Xx", length(xyz)/3)
  }
    
  if(length(elesy) != length(xyz)/3)
    stop("'xyz' and 'elesy' must have matching lengths")
  
  data("elements", package = "bio3d", envir=environment()) 
  elements <- get("elements", envir=environment()) 
  
  M <- match(elesy,elements[,"symb"])
  M[is.na(M)] <- 1 # Unrecognized elements are considered as dummy atoms    
  
  if(is.null(col))
    col <- do.call(rgb, elements[M, c("red","green","blue")])
  if(length(col) != ncol(xyz)/3){
    if(length(col) != 1)
      warning("'col' has been recycled")
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

  if(requireNamespace(package = "rgl", quietly = TRUE)){
    if(!add){
      rgl::open3d()
      rgl::par3d(windowRect = windowRect, userMatrix=userMatrix, FOV = FOV, ...)
      rgl::bg3d(color=bg.col)
    }
    ids <- rgl::rgl.ids()
    
    if(xyz.axes) ids <- rbind(ids, addXYZ(lwd = lwd.xyz, cex = cex.xyz))
    if(abc.axes) ids <- rbind(ids, addABC(cell, lwd = lwd.abc, cex = cex.abc)) ## <-- MISSING FUN!!
    if(pbc.box ) ids <- rbind(ids, addPBCBox(cell, lwd = lwd.pbc.box))
    
    if(nchar(type)>1)
      type <- strsplit(type, "")[[1]]
    if(!all(type %in% c("l","s","p")))
      stop("Unrecognized 'type'")
    
    if("l" %in% type) {
      if(is.null(con)) {
        warning("Unspecifyed connectivity: Computing connectivity from coordinates...")
        con <- connectivity.xyz(x = xyz, elesy = elesy, by.block = TRUE)
      }
      if(!is.null(con)){
        ind <- t(con)
        seg.id <- rgl::segments3d(
          xyz[seq(1,ncol(xyz),3)][ind],
          xyz[seq(2,ncol(xyz),3)][ind],
          xyz[seq(3,ncol(xyz),3)][ind],
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
      sph.id <- rgl::spheres3d(
        xyz[seq(1,length(xyz),3)],
        xyz[seq(2,length(xyz),3)],
        xyz[seq(3,length(xyz),3)],
        color = col, radius = radii, ...)
      sph.id <- data.frame(id = sph.id, type= "atom.sph")
      ids <- rbind(ids, sph.id)
    }
    if("p" %in% type) {
      pts.id <- rgl::points3d(
        xyz[seq(1,length(xyz),3)],
        xyz[seq(2,length(xyz),3)],
        xyz[seq(3,length(xyz),3)],
        color = col, ...)
      pts.id <- data.frame(id = pts.id, type = "atom.pts")
      ids <- rbind(ids, pts.id)
    }
    if(centre) {
      cent <- centres(xyz, w=rep(1, length(xyz)/3))
      userMatrix <- t(rgl::translationMatrix(x = -cent[1], y = -cent[2], z = -cent[3]))
      rgl::par3d(userMatrix=userMatrix)
    }
    invisible(ids)
  } else {
    stop("Please install the rgl package from CRAN")
  }
}

visualize.pdb <- function(
  pdb, cell = pdb$cell, type = "l", con = TRUE,
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rvdw", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...){
  
  if(!is.pdb(pdb))
    stop("'pdb' must be an object of class pdb. See read.pdb")
  
#   if(any(!are.symb(pdb$atom$elesy) & is.na(pdb$atom$elesy)))
#     stop("'pdb' contains unvalid 'elesy'")
  
  if(grepl("l", type) & (is.null(pdb$con) | con)) {
    cat("Computing connectivity from coordinates...\n")
    pdb$con <- connectivity(pdb)
  }
  pdb$con$eleno.1 <- match(pdb$con$eleno.1, pdb$atom$eleno)
  pdb$con$eleno.2 <- match(pdb$con$eleno.2, pdb$atom$eleno)
  visualize.xyz(
    pdb$xyz, elesy = pdb$atom$elesy, pdb$con, cell, type,
    xyz.axes, abc.axes, pbc.box, centre, lwd, lwd.xyz, lwd.abc, lwd.pbc.box,
    cex.xyz, cex.abc, col, radii, bg.col, add, windowRect,
    userMatrix, FOV, ...)
}

visualize.character <- function(
  file, cell = NULL, type = "l", con = TRUE,
  xyz.axes = FALSE, abc.axes = FALSE, pbc.box = FALSE, centre = TRUE,
  lwd = 2, lwd.xyz = lwd, lwd.abc = lwd, lwd.pbc.box = lwd,
  cex.xyz = 2, cex.abc = 2, col = NULL, radii = "rcov", bg.col = "black",
  add = FALSE, windowRect = c(0,0,800,600), userMatrix=diag(4), FOV = 0, ...) {
  
  x <- read.pdb(file)
  if(is.null(cell))
    cell <- x$cell
  visualize.pdb(
    x, cell, type, con, xyz.axes, abc.axes, pbc.box, centre,
    lwd, lwd.xyz, lwd.abc, lwd.pbc.box, cex.xyz, cex.abc, col,
    radii, bg.col, add, windowRect, userMatrix, FOV, ...)
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
#   rgl::spheres3d(membership.centres, col = com.net.vertex.color, radius = radii, alpha = 0.5)
#   
#   com.net.con <- apply(get.edgelist(net$community.network), 2, as.integer)
#   cyls <- apply(com.net.con, 1,
#     function(ids){
#       cyl <- cylinder3d(membership.centres[ids,], radius=com.net.weight, sides=20)
#       return(cyl)
#     })
#   cyls.ids <- lapply(cyls, shade3d, col = "white")
# }
