
as.pdb.mol2 <- function(mol, ...) {

    cl <- match.call()
    natoms <- nrow(mol$atom)
    pdbn <- list()

    allset <- FALSE
    resid <- strtrim(mol$atom$resid, 3)
    chain <- NULL
    
    if(!is.null(mol$substructure)) {
        if(all(c("root_atom", "subst_type", "chain") %in% colnames(mol$substructure))) {
            key1 <- paste(mol$atom$resno, mol$atom$resid, sep="-")
            rownames(mol$substructure) <- key1[ mol$subs$root_atom ]
            resid <- mol$substructure[key1, "sub_type"]
            chain <- mol$substructure[key1, "chain"]
            allset <- TRUE
        }
    }

    if(!allset & length(unique(mol$resid)) > 1) {
        warning("insuffient data in SUBSTRUCTURE to set residue and chain identifiers in PDB")
    }

    pdb <- as.pdb.default(pdb    = NULL,
                          xyz    = mol$xyz,
                          type   = rep("ATOM", natoms),
                          resno  = mol$atom$resno,
                          resid  = resid,
                          eleno  =  mol$eleno,
                          elety  = mol$atom$elena,
                          chain  = chain,
                          insert = NULL,
                          alt    = NULL,
                          o=NULL, b=NULL, segid=NULL,
                          elesy=unlist(lapply(strsplit(mol$atom$elety, split="[.]"),
                                                function(x) x[1])),
                          charge=mol$atom$charge)
    
    pdb$call <- cl
    return(pdb)
}
