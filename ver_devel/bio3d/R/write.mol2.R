write.mol2 <- function(mol, file="R.mol2", append=TRUE) {

    if(!is.mol2(mol))
        stop("input should be of class 'mol2' as obtained by 'read.mol2'")

    if(any(is.na(mol$atom))) {
        mol$atom[ is.na(mol$atom) ] = ""
    }
    if(any(is.na(mol$bond))) {
        mol$bond[ is.na(mol$bond) ] = ""
    }
    if(any(is.na(mol$substructure))) {
        mol$substructure[ is.na(mol$substructure) ] = ""
    }
  
    raw.lines <- c()
    
    raw.lines <- c(raw.lines, "@<TRIPOS>MOLECULE")
    raw.lines <- c(raw.lines, mol$name)

    fmt <- "%7s %7s %7s"
    mol$info[is.na(mol$info)]=" "
    raw.lines <- c(raw.lines, 
                   sprintf(fmt, mol$info[1], mol$info[2], mol$info[3]),
                   "")

    
    
    raw.lines <- c(raw.lines, "@<TRIPOS>ATOM")
    fmt <- "%7s %-8s %9.4f %9.4f %9.4f %-5s %5s %-9s %8.4f "
    for ( i in 1:nrow(mol$atom) ) {
        raw.lines <- c(raw.lines,
                       paste0(
                           sprintf(fmt, mol$atom[i, 1], mol$atom[i, 2],
                                   mol$atom[i, 3], mol$atom[i, 4],
                                   mol$atom[i, 5], mol$atom[i, 6],
                                   mol$atom[i, 7], mol$atom[i, 8],
                                   mol$atom[i, 9], mol$atom[i, 10]
                                   ),  mol$atom[1, 10])
                       )
    }


    raw.lines <- c(raw.lines, "@<TRIPOS>BOND")
    fmt <- "%7s %5s %5s %-5s"
    for ( i in 1:nrow(mol$bond) ) {
        raw.lines <- c(raw.lines, 
                       paste0(
                           sprintf(fmt,
                                   mol$bond[i, 1], mol$bond[i, 2],
                                   mol$bond[i, 3], mol$bond[i, 4]
                                   ), mol$bond[1, 5])
                           )
    }
    
    raw.lines <- c(raw.lines, "@<TRIPOS>SUBSTRUCTURE")
    fmt <- "%7s %-8s %-10s %6s %4s %-6s %4s"
    for ( i in 1:nrow(mol$substructure) ) {
        raw.lines <- c(raw.lines, 
                       paste0(
                           sprintf(fmt,
                                   mol$substructure[i, 1], mol$substructure[i, 2],
                                   mol$substructure[i, 3], mol$substructure[i, 4],
                                   mol$substructure[i, 5], mol$substructure[i, 6],
                                   mol$substructure[i, 7], mol$substructure[i, 8]
                                   ), mol$substructure[1, 9])
                           )
    }

 
    
    
    
    write.table(raw.lines, file = file, quote = FALSE, na="",
                row.names = FALSE, col.names = FALSE, append = append)
    

}
