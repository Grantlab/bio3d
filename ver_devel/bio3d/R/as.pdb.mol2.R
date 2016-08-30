
as.pdb.mol2 <- function(mol2, ...) {

    cl <- match.call()
    natoms <- nrow(mol2$atom)
    xyz <- mol2$xyz
    tmp.pdb <- list()
    
    if(!is.null(mol2$substructure)) {
        rownames(mol2$substructure) <- mol2$substructure$name
        resid <- mol2$substructure[mol2$atom$resid, "sub_type"]
        chain <- mol2$substructure[mol2$atom$resid, "chain"]
    }
    else {
        resid <- "UNK"
        chain <- NA
    }
    
    tmp.pdb$atom <- data.frame(cbind(type=rep("ATOM", natoms),
                                     eleno=seq(1, natoms),
                                     elety=mol2$atom$elena,
                                     alt=NA,
                                     resid=resid,
                                     chain=chain, 
                                     resno=mol2$atom$resno,
                                     insert=NA,
                                     x=mol2$atom$x, y=mol2$atom$y, z=mol2$atom$z,
                                     o=NA, b=NA, segid=NA,
                                     elesy=unlist(lapply(strsplit(mol2$atom$elety, split="[.]"),
                                                         function(x) x[1])),
                                     charge=mol2$atom$charge),
                               stringsAsFactors=FALSE)
    
    tmp.pdb$xyz <- xyz
    class(tmp.pdb) <- "pdb"
    tmp.pdb$calpha <- rep(FALSE, nrow(tmp.pdb$atom))
    tmp.pdb$call <- cl
    
    return(tmp.pdb)
}
