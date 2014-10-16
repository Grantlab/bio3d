"write.dcd" <-
function(xyz, header = NULL, file = "trj.dcd", verbose = TRUE, ...) {

   # Description:
   #  Writes a CHARMM or X-PLOR/NAMD binary
   #  trajectory file with either big- or
   #  little-endian storage formats
   #

   #===DCD=FORMAT==============================================
   #HDR  NSET  ISTRT NSAVC 5-ZEROS NATOM-NFREAT DELTA   9-ZEROS
   #CORD files step1 step  zeroes  (zero)      timestep  zeroes
   #C*4  INT   INT   INT    5INT    INT         DOUBLE   9INT
   #  [CHARACTER*20]
   #===========================================================
   #NTITLE   TITLE
   #INT      C*MAXTITL
   #C*2      C*80
   #===========================================================
   #NATOM
   #INT
   #===========================================================
   #CELL(I), I=1,6          (DOUBLE)
   #===========================================================
   #X(I), I=1,NATOM         (SINGLE)
   #Y(I), I=1,NATOM         
   #Z(I), I=1,NATOM         
   #===========================================================


   # Open file conection
   fout <- file(file, "wb")
 
   if( is.null(header) ) { 

     # Default header 
     header <- list(84L, 
         hdr = "CORD",
         nframe = as.integer(nrow(xyz)),
         first = 0L,
         step = 1L,
         nstep = as.integer(nrow(xyz)),
         integer(3),
         ndegf = 0L,
         nfixed = 0L,
         delta = 0,
         cryst = 0L,
         block = 0L,
         integer(7),
         vers = 24L,
         84L,
         84L,
         ntitle = 1L,
         title = sprintf("%-80s", "Automatically generated header by Bio3D 2.0"),
         84L,
         4L,
         natom = as.integer(ncol(xyz) / 3),
         4L
     )
   } 
  
   # Write the header 
   writeBin(header[[1]], fout, ...)
   writeChar(header$hdr, fout, eos=NULL)
   for(i in 3:9) writeBin(header[[i]], fout, ...)
   writeBin(header$delta, fout, size = 4, ...)
   for(i in 11:17) writeBin(header[[i]], fout, ...)
   writeChar(header$title, fout, eos=NULL)
   for(i in 19:22) writeBin(header[[i]], fout, ...) 
   
   # Write XYZ
   apply(xyz, 1, function(x) {
      writeBin(header$natom * 4L, fout, ...)
      writeBin(x[seq(1, length(x), 3)], fout, size = 4, ...)
      writeBin(header$natom * 4L, fout, ...)
 
      writeBin(header$natom * 4L, fout, ...)
      writeBin(x[seq(2, length(x), 3)], fout, size = 4, ...)
      writeBin(header$natom * 4L, fout, ...)
 
      writeBin(header$natom * 4L, fout, ...)
      writeBin(x[seq(3, length(x), 3)], fout, size = 4, ...)
      writeBin(header$natom * 4L, fout, ...)
   } )

   invisible(close(fout))
}

