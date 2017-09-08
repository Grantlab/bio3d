.get.exepath <- function(prg) {
    
    paths <- list(

        pymol = list(
            
            Linux = c("/usr/bin/pymol",
                      "/usr/local/bin/pymol"),
            
            Darwin = c("/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL",
                       "/usr/bin/pymol",
                       "/usr/local/bin/pymol"),
            
            Windows = c("C:/python27/PyMOL/pymol.exe",
                        "C:/Program Files/PyMOL/PyMOL/PymolWin.exe",
                        "C:/Program Files/PyMOL/PymolWin.exe"),

            ver = "-cq"
        ),


        muscle = list(

            Linux = c("/usr/bin/muscle",
                      "/usr/local/bin/muscle"),
            
            Darwin = c("/usr/bin/muscle",
                       "/usr/local/bin/muscle"),
            
            Windows = c("C:/Program Files/muscle.exe",
                        "C:/Program Files/muscle3.8.31_i86win32.exe",
                        "C:/Program Files/muscle/muscle.exe",
                        "C:/Program Files/Muscle/muscle.exe",
                        "C:/Program Files/seaview/muscle.exe",
                        "C:/Program Files/seaview4/muscle.exe"),

            ver = "-version"
        ),

        clustalo = list(

            Linux = c("/usr/bin/clustalo",
                      "/usr/local/bin/clustalo"),
            
            Darwin = c("/usr/bin/clustalo",
                       "/usr/local/bin/clustalo"),
            
            Windows = c("C:/Program Files/clustalo.exe", 
                        "C:/Program Files/clustalo/clustalo.exe",
                        "C:/Program Files/Clustalo/clustalo.exe",
                        "C:/Program Files/seaview/clustalo.exe",
                        "C:/Program Files/seaview4/clustalo.exe"),

            ver = "--version"
        ), 
        
        dssp = list(

            Linux = c("/usr/bin/dssp",
                      "/usr/local/bin/dssp"),
            
            Darwin = c("/usr/bin/dssp",
                       "/usr/local/bin/dssp",
                       "/usr/bin/mkdssp",
                       "/usr/local/bin/mkdssp"),
            
            Windows = c("C:/Program Files/dssp.exe",
                        "C:/Program Files/dssp-2.0.4-win32.exe",
                        "C:/Program Files/dssp/dssp.exe",
                        "C:/Program Files/Dssp/dssp.exe"),
            
            ver = "--version"
        )
 
    )

    ## user provided full path
    if(file.exists(prg) & !dir.exists(prg)) {
        return(prg)
    }
    
    ## try to automatically determine path
    exefile <- Sys.which(prg)

    if(nchar(exefile) == 0) {

        if(prg %in% c("pymol", "muscle", "clustalo", "dssp")) {
            ## determine os
            os1 <- Sys.info()["sysname"]
            
            ## use guess-paths defined above
            exefiles <- paths[[prg]][[os1]]
            fe <- file.exists(exefiles)
            
            if(any(fe)) {
                exefile <- exefiles[which(fe)[1]]
            }
            else {
                exefile <- NULL
            }
        }
        else {
            exefile <- NULL
        }
    }
    
    if(is.null(exefile)) {
        stop(paste0("could not determine path to '", prg, "'"))
    }
    return(exefile)
}

.test.exefile <- function(exefile) {
    prg <- tolower(basename(exefile))

    if(grepl("muscle", prg)) {
        ver <- "-version"
    }
    if(grepl("pymol", prg)) {
        ver <- "-cq"
    }
    if(grepl("clustalo", prg)) {
        ver <- "--version"
    }
    if(grepl("dssp", prg)) {
        ver <- "--version"
    }

    
    os1 <- Sys.info()["sysname"]
    if (os1 == "Windows") {
        success <- shell(paste(shQuote(exefile), ver))
    }
    else {
        success <- system(paste(exefile, ver),
                          ignore.stderr = TRUE, ignore.stdout = TRUE)
    }
    
    if(!(success %in% c(0,1))) {
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}
    
