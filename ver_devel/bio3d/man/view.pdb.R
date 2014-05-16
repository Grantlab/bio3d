sse.color <- function(x, col.coil = "gray", col.helix = "purple", col.sheet = "blue") {
  ##-- Define a color vector with helix and sheet 
  ##    annotations taken from input PDB file.
  ##      ToDo: - if no $helix and $sheet defined take 
  ## 		        annotation from dssp() or stride()
  
  col <- rep(col.coil, nrow(x$atom))
  
  if(is.null(x$helix))
    stop("'x$helix' is NULL. Please use 'dssp' or 'stride' to create it.")
  else  
    h.resno <- unbound(x$helix$start, x$helix$end)
  
  if(is.null(x$sheet))
    stop("'x$sheet' is NULL. Please use 'dssp' or 'stride' to create it.")
  else
    e.resno <- unbound(x$sheet$start, x$sheet$end)
  
  h.ind <- which(as.numeric(x$atom$ersno) %in% h.resno)
  e.ind <- which(as.numeric(x$atom$resno) %in% e.resno)
  
  col[h.ind] <- col.helix
  col[e.ind] <- col.sheet
  return(col)
}

view <- function(...)
  UseMethod("view")

view.pdb <- function(pdb, sel = "default", col = NULL, cna = NULL, ...) {
    ##-- Wrapper for visualize() to view larger PDBs the way Barry 
    ##    likes to see them most often.
    ##      To Do - This is a very quick and dirty prototype with 
    ##              no consideration of efficiency or error checking. 
    ##    For example, we use trim.pdb() rather than just using atom 
    ##    selections and we don't check the 'col' input vector length 
    ##    matches natom or convert it based on 'sel' input).
    ##
#     par.save <- par3d(skipRedraw=TRUE)
#     on.exit(par3d(par.save))
    
    if(sel == "default") {
        ca.sel <- atom.select(pdb, "calpha" , verbose = FALSE)
      prot.sel <- atom.select(pdb, "protein", verbose = FALSE)
      back.sel <- atom.select(pdb, "back"   , verbose = FALSE)
       lig.sel <- atom.select(pdb, "ligand" , verbose = FALSE)

      connectivity(pdb) <- connectivity(pdb)

      ## Ligand
      if(length(lig.sel$atom) != 0){
        pdb.lig <- trim.pdb(pdb, lig.sel)
        visualize(pdb.lig, con = FALSE, type="s") 
      }

      ## Sidechain
      tot.sel <- combine.sel(prot.sel, back.sel, op="NOT",verbose=FALSE)
      tot.sel <- combine.sel( tot.sel,   ca.sel, op="OR" ,verbose=FALSE)
      if(length(tot.sel$atom) != 0) {
        pdb <- trim.pdb(pdb, tot.sel)
        visualize(pdb, con = FALSE, add = TRUE, type = "l", col = "gray", lwd = 1)        
      }
     
      # Calpha
      ca.sel <- atom.select(pdb, "calpha" , verbose = FALSE)
      if(length(ca.sel$atom) != 0) {
        pdb <- trim.pdb(pdb, ca.sel)        
        connectivity(pdb) <- calpha.connectivity(pdb)
        if(is.null(col))
          col <- sse.color(pdb, col.coil = "#808080")
        visualize(pdb, con = FALSE, add = TRUE, type = "l", col = col)
      }
    }
    if(sel=="calpha") { ## I like this view!
      ## draw sse colored calpha trace
      pdb.ca <- trim.pdb(pdb, atom.select(pdb,"calpha",verbose=FALSE))
      connectivity(pdb.ca) <- calpha.connectivity(pdb.ca)
      if(is.null(col))
        col <- sse.color(pdb.ca)
      if(length(pdb.ca$xyz)!=0)
        visualize(pdb.ca, con=FALSE, col=col, ...)
    } 
    
    if(sel=="protein") {
      pdb.prot <- trim.pdb(pdb, atom.select(pdb,"protein",verbose=FALSE))
      connectivity(pdb.prot) <- connectivity(pdb.prot)
      if(length(pdb.prot$xyz)!=0)
        visualize(pdb.prot, con = FALSE, ...)
    }
    if(sel=="back") {
      back <- combine.sel( atom.select(pdb,"protein",verbose=FALSE), atom.select(pdb,"back",verbose=FALSE),verbose=FALSE )
      pdb.prot <- trim.pdb(pdb, back)
      connectivity(pdb.prot) <- connectivity(pdb.prot)
      if(is.null(col))
        col <- sse.color(pdb.prot)
      if(length(pdb.prot$xyz)!=0)
        visualize(pdb.prot, con = FALSE, col=col, ...)
    }
    
    #if(!is.null(cna)) { 
    #		visualize.cna(cna, pdb, ...) ## need to think more about this
    #}

}
  
view.character <- function(file, sel = "default", col = NULL, cna = NULL, ...) {
  x <- read.pdb(file)
  view.pdb(x, sel = sel, col = col, cna = cna, ...)
}
