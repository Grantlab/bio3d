write.mol2 <- function(mol, file="R.mol2", append=FALSE) {

    if(!is.mol2(mol))
        stop("input should be of class 'mol2' as obtained by 'read.mol2'")

    if(any(is.na(mol$atom))) {
        mol$atom[ is.na(mol$atom) ] = ""
    }
    if(!is.null(mol$bond)) {
        if(any(is.na(mol$bond)))
            mol$bond[ is.na(mol$bond) ] = ""
    }
    if(!is.null(mol$substructure)) {
        if(any(is.na(mol$substructure)))
            mol$substructure[ is.na(mol$substructure) ] = ""
    }
  
    raw.lines <- c()
    
    raw.lines <- c(raw.lines, "@<TRIPOS>MOLECULE")
    raw.lines <- c(raw.lines, mol$name)

    if(length(mol$info) < 5) mol$info <- c(mol$info, rep(NA, 5-length(mol$info)))
    mol$info[is.na(mol$info)] = ""
    fmt <- paste(rep("%7s", length(mol$info)), collapse=" ")
    raw.lines <- c(raw.lines,
                   sprintf(fmt, mol$info[1], mol$info[2], mol$info[3], mol$info[4], mol$info[5])
                   )
    
    
    raw.lines <- c(raw.lines, "@<TRIPOS>ATOM")
    fmt <- "%7s %-8s %9.4f %9.4f %9.4f %-5s %5s %-9s %8.4f "
    for ( i in 1:nrow(mol$atom) ) {
        raw.lines <- c(raw.lines,
                       paste0(
                           sprintf(fmt, mol$atom[i, 1], mol$atom[i, 2],
                                   mol$atom[i, 3], mol$atom[i, 4],
                                   mol$atom[i, 5], mol$atom[i, 6],
                                   mol$atom[i, 7], mol$atom[i, 8],
                                   mol$atom[i, 9] 
                                   ),  mol$atom[i, 10])
                       )
    }

    if(!is.null(mol$bond)) {
        raw.lines <- c(raw.lines, "@<TRIPOS>BOND")
        fmt <- "%7s %5s %5s %-5s "
        for ( i in 1:nrow(mol$bond) ) {
            raw.lines <- c(raw.lines, 
                           paste0(
                               sprintf(fmt,
                                       mol$bond[i, 1], mol$bond[i, 2],
                                       mol$bond[i, 3], mol$bond[i, 4]
                                       ), mol$bond[i, 5])
                           )
        }
    }

    if(!is.null(mol$substructure)) {
        raw.lines <- c(raw.lines, "@<TRIPOS>SUBSTRUCTURE")
        fmt <- "%7s %-8s %-10s %8s %4s %-6s %4s %4s "
        for ( i in 1:nrow(mol$substructure) ) {
            raw.lines <- c(raw.lines, 
                           paste0(
                               sprintf(fmt,
                                       mol$substructure[i, 1], mol$substructure[i, 2],
                                       mol$substructure[i, 3], mol$substructure[i, 4],
                                       mol$substructure[i, 5], mol$substructure[i, 6],
                                   mol$substructure[i, 7], mol$substructure[i, 8]
                                   ), mol$substructure[i, 9])
                           )
        }
    }
     
    
    write.table(raw.lines, file = file, quote = FALSE, na="",
                row.names = FALSE, col.names = FALSE, append = append)
    

}
