
## UNREFINED functions for Bio3D incoperation at some point
##  Most of these are useful for protein structure modeling
##
## These include (list to be updaed):
##  # fit.pdbs    - Quick Fit Fitter for PDBs
##  pdbname     - Extract PDB identifier from filename
##  srxn.bd     - Find BD trajectories that complete a give reaction
##  srxn2trj.bd - Make XML BD trajectorys given srxn output
##  motif.find  - Return indices of a motif within a sequence
##  getArgs     - Parse command line options when using Rscript
##  txt2num     - Convert a character string to numeric
##  read.apbs   - Read elec binding energy from APBS log files
##  chain.pdb   - Find possible chian breaks
##  pdbaln      - Quick and dirty alignment of PDB sequences
##  alitrim     - Trim cols from alignment data structure
##   ncmap      - *see bio3d 'cmap'
##   ndm        - *see bio3d 'dm.xyz'
##  seq2aln     - Add a sequence to an existing alignment
##  write.pir   - write alignment for modeler
##  write.sge   - write a series of Sun Grid Engine scripts
##                for running AMBER on chemcca2, 31 or oolite
##  interp      - interpolate between two vectors, useful for
##                PCA z-score interpolation in combination
##                with "pca.z2xyz"
##  fit         - A simple wrapper for fit.xyz when using
##                pdbs style objects
##  ide.group   - Return the indices of the largest group of
##                sequences that have identity values above
##                a particular 'cutoff'
##  rama.inds   - Return xyz indices for PHI-PSI Ramachendran atoms
##  plot.rama   - Ramachendron plot (basic)
##  aln2aln     - Add one alignment to another that contains
##                at least one similar entry
##  get.pdb     - download PDB files from a list of ids
##  get.uniprot - download FASTA sequence files from a list of
##                swissprot or uniprot ids
##  seq.pdb     - Return basic 1-letter calpha ATOM sequence from a
##                pdb object
## read.propka  - Read the output of PropKa produced by pdb2pqr
## write.crdbox - Write AMBER CRD format trajectory files
## blast.n      - Blast with nucleotide sequences
## gb2fasta     - Convert GENEBANK to FASTA format
## free2fasta   - Convert FREE (SCA) to FASTA format
## fasta2free   - Convert FASTA to FREE (SCA) format. 
## tlsq & ulsq  - Fitting rotation and translation matrices
## renumber.pdb - Renumber resno and eleno records
## vec2seq      - Match a vector via matching sequence to alignment
## vec2resno    - Replicate a vector based on concetive resno entries
## bgr.colors   - blue-gray-red color range
## lsos         - a better list of current objects sorted by size
## read.hmmer.tbl - read HMMER3 hmmsearch log files
## cons.aln     - score residue conservation in an alignment
## dis.ftmap    - Process FTMAP results (min dist to residue of probes)
## bootstrap.rmsf - Bootstrap sampling of frames for RMSF determination.
## mustang        -  Structural alignment with mustang
## add.dccm.grid  -  Add a grid or colored boxes to a plot.dccm() plot.
## col.wheel    - useful for picking plot colors (e.g. col.wheel("dark") ) 
##
##
## See also:
##  ~/work/scop/scop.sf  - Have access to the full SCOP database in R
##  ~/tmpwork/Rpackage/new_funs/HB.back.R
##  ~/tmpwork/Rpackage/new_funs/Hbond.R
##  ~/tmpwork/Rpackage/new_funs/io.R
##  ~/tmpwork/Rpackage/new_funs/plot.cij.R
##  ~/tmpwork/Rpackage/new_funs/read.prmtop.r
##  ~/tmpwork/Rpackage/new_funs/read.inpcrd.r
##  ~/tmpwork/Rpackage/new_funs/read.psf.R
##  ~/tmpwork/Rpackage/new_funs/read.pqr.R
##  ~/tmpwork/Rpackage/new_funs/read.pdb.R (with resolution etc)
##  ~/tmpwork/Rpackage/new_funs/rot_trans/rot.trans.R



### --- See bio3d pdbfit()
###fit.pdbs <- function(pdbs, inds=NULL, outpath=NULL, ...) {
###  ##
###  ## Quick Fit Fitter for PDBs
###  ##
###  if(class(pdbs)!="3dalign") {
###    stop("Input 'pdbs' should be of class '3dalign'")
###  }
###  full <- ifelse(is.null(outpath), FALSE, TRUE)
###  if(is.null(inds)) {    inds <-gap.inspect(pdbs$xyz)$f.inds  }
###  if(is.list(inds)){ inds=inds$xyz }
###  return( fit.xyz( fixed=pdbs$xyz[1,], mobile=pdbs, fixed.inds=inds,
###                  mobile.inds=inds, outpath=outpath,
###                  full.pdbs=full, ... ))
###}


###pdbname <- function(x, mk4=TRUE) {
###  ##
###  ## Extract PDB identifer from filenames
###  ## like "basename()" for PDB files
###  ##
###  x <- sub(".pdb$","", basename(x))
###  if(mk4) {
###    return(substr(x,1,4))
###  } else { return(x) }
###}

pdbname <- function(x, pdbid=FALSE) {
  ## Extract PDB name or ID from file names
  ## like basename() but for PDB files
  y <- sub("\\.pdb$", "", basename(x))
  if(pdbid)
    y <- substr(y, 1, 4)

  names(y) = x
  return(y)
}


srxn.bd <- function(traj.file, reaction, outdir="trjout/") {
  ##
  ## Find BD trajectories that complete a give reaction
  ##

  ## srx <- srxn.bd("traj0", reaction=8)
  ## xml <- srxn2trj.bd(srx, dcd=FALSE)

  dir.create(outdir, FALSE)

  trj.name <- basename(traj.file)
  log.file  <- paste(outdir,"/tlog_rxn_",reaction,
                     "_", trj.name, ".xml", sep="")

  cmd1 <- paste("~bgrant/software/browndye/bin/process_trajectories -traj ",
                traj.file,".xml -index ",traj.file,
                ".index.xml -srxn ", reaction,
                " > ", log.file, sep="")


  srxlog <- function(logfile) {

    htrj <- scan(logfile, what="", sep="\n", quiet=TRUE)
    hind <- cbind( grep("<trajectory>", htrj),
                  grep("<number>", htrj),
                  grep("<subtrajectory>", htrj) )
    hval <- cbind(as.numeric( gsub("..number>","", htrj[hind[,2]]) ),
                  as.numeric( gsub("..subtrajectory>","", htrj[hind[,3]]) ))
    colnames(hval) <- c("num","sub")
    return(hval)
  }

  cat(paste("\n Creating SRXN file:\n\t",log.file,"\n\n"))
  system(cmd1)
  hval <- srxlog(log.file)
  return(list(numb=hval[,"num"], nsub=hval[,"sub"],
              logfile=log.file, cmd=cmd1,
              traj.file=traj.file, reaction=reaction) )
}




srxn2trj.bd <- function(srx,
                        Aatoms.xml="../m-atoms.xml",
                        Batoms.xml="../k-atoms.xml",
                        xyz=TRUE, dcd=TRUE,
                        outdir="trjout/") {

  ##
  ## Make XML BD trajectorys given srxn output
  ##

  ## srx <- srxn.bd("traj0", reaction=8)
  ## xml <- srxn2trj.bd(srx, dcd=FALSE)

  dir.create(outdir, FALSE)
  
  xml.file <- paste(outdir,"/trxn_", srx$reaction,
                    "_", basename(srx$traj.file),
                    ".", srx$numb, ".",
                     srx$nsub, ".xml", sep="")

  xyz.file <- sub(".xml$", ".xyz", xml.file)

  dcd.file <- sub(".xml$", ".dcd", xml.file)

  cmd2 <- paste(##"~/software/browndye/bin/",
                "process_trajectories -traj ",
                srx$traj.file,".xml -index ",srx$traj.file,
                ".index.xml -n ", srx$numb,
                " -sn ", srx$nsub, " > ",
                xml.file, sep="")

  cat(paste("\n Creating XML TRJ file(s):\n\t",
            paste(xml.file, collapse="\n\t "), "\n\n"))

  for(cmd in cmd2) { system(cmd) }

  
  ## xyz_trajectory
  cmd3 <- paste(##"~/software/tmp/browndye/bin/",
                "xyz_trajectory ",
                "-mol0 ", Aatoms.xml," -mol1 ", Batoms.xml,
                " -trajf ", xml.file," > ", xyz.file, sep="")
              
  if(xyz){
    cat(paste("\n Creating XYZ TRJ file(s):\n\t",
              paste(xyz.file, collapse="\n\t "), "\n\n"))
  
    for(cmd in cmd3) { system(cmd) }
  }

  
  ## catdcd
  cmd4 <- paste("catdcd -o ", dcd.file,
                " -otype dcd -xyz ", xyz.file, sep="")
    
  if(dcd){
    cat(paste("\n Creating DCD TRJ file(s):\n\t",
              paste(dcd.file, collapse="\n\t "), "\n\n"))

    for(cmd in cmd4) { system(cmd) }
    
  }

  
  return( list(xml.file=xml.file,
               xyz.file=xyz.file,
               dcd.file=dcd.file,
               cmd.xml=cmd2,
               cmd.xyz=cmd3,
               cmd.dcd=cmd4,
               numb=srx$numb, nsub=srx$nsub,
               traj.file=srx$traj.file, reaction=srx$reaction) )
}


### --- See bio3d motif.find()
###motif.find <- function(motif, sequence) {
###  ## return indices of motif within sequence
###  position <- regexpr( paste(motif, collapse=""), paste(sequence,collapse=""))
###  inds <- c(position):c(position+attr(position, "match.length")-1)
###  return(inds)
###}

getArgs <- function(set) {

  ## Parse command line options (for Rscript IO handeling)
  ## If there is an equals sign (=) in the input args then it
  ## is parsed and the first part (before the = sign) treated
  ## as the variable name, with the second (after the = sign )
  ## taken as the value of the argument.
  ## If there is no = sign then the arg value is assigned
  ## logical TRUE
  ##
  ## Note that all values will be characters (i.e. strings).
  ## so YOU must sort them out afterwords
  ##
  ## E.g.
  ##
  ##  #!/usr/bin/env Rscript
  ##  got <- getArgs( c("dir","id","nfile") )
  ##
  ##  if(!got["dir"])
  ##    stop("Command line argument 'dir=?' is required")
  ##  if(!got["id"]) {
  ##    warning("Assuming 'id' is the same as basename(dir)")
  ##    id <- basename(dir) }
  ##  if(!got["nfile"])
  ##    nfile <- NULL
  ##  ...


  passed.args <- commandArgs(TRUE)
  cat( paste("   Called with", length(passed.args),"arguments:\n\t-  ",
             paste(passed.args, collapse=",\n\t-   ")),"\n\n")

  for (e in passed.args) {
    ta = strsplit(e,"=",fixed=TRUE)
    if(ta[[1]][1] %in% set) {
      if(! is.na(ta[[1]][2])) {
        assign(ta[[1]][1], ta[[1]][2], env = .GlobalEnv)
        cat("   Assigned ",ta[[1]][1]," the value of:",ta[[1]][2],"\n")
      } else {
        assign(ta[[1]][1],TRUE, env = .GlobalEnv)
        cat("   Assigned ",ta[[1]][1]," the value of: TRUE\n")
      }
    } else {
      warning(paste(" Command line argument:", ta[[1]][1],
                    "is undefined and will be ignored"))
    }
  }
  sapply(set, exists)
}




txt2num <- function(x) {
  ## Convert a character string to numeric
  ## txt2num( "1:4,6,lkk,1:10,fdsfsd,99a" )

  if(length(x)==0)
    return(x)
  
  num1 <- unlist(strsplit(x, split = ","))
  num2 <- suppressWarnings(as.numeric(num1))

  nsplit <- grep(":",num1)
  num3 <- unlist(strsplit( num1[ grep(":",num1) ], split = ":"))
  num4 <- apply( matrix(as.numeric(num3),nrow=2),
                2, function(x){ x[1]:x[2] })
  
  num <- NULL; j <- 1
  for(i in 1:length(num1)) {
    if(!is.na(num2[i])) {
      num <- c(num, num2[i])
    } else {
      if( i %in% nsplit) {
        if( ncol(num4)==1 ) {
          num <- c(num, num4[,j])
        } else {
          num <- c(num, num4[[j]])
        }
          j <- j+1
      } else {
        cat(paste(" Ignoring character component: ",num1[i],"\n")) 
      }
    }
  }
  return(num)
}


### --- See bio3d::chain.pdb()
###chain.break <- function(pdb, ca.dist=4, blank="X", rtn.vec=TRUE) {
###  ##- Find possible chian breaks
###  ##  (i.e. Caplpa's that are too far apart).
###  ##
###  ## chn <- chain.break(pdb)
###  x <- diff( as.numeric(pdb$atom[pdb$calpha,"x"]) )^2
###  y <- diff( as.numeric(pdb$atom[pdb$calpha,"y"]) )^2
###  z <- diff( as.numeric(pdb$atom[pdb$calpha,"z"]) )^2
###  d <- sqrt((y+y+z)/3)
###
###  ind <- which(d > ca.dist)
###  len <- diff( c(1,ind,length(d)) )
###  ## vec <- rep( LETTERS[1:length(len)], len)
###
###  cat(paste("Found",length(ind), " possible chain breaks\n"))
###  cat(paste("After resno(s)",
###        paste( pdb$atom[pdb$calpha,"resno"][(ind)], collapse=", " ),"\n" ))
###  cat(paste("chain length(s)",
###        paste(len, collapse=", " ),"\n" ))
###
###  if(rtn.vec) {
###    ## Make a chain id vector
###    resno.ind <- as.numeric(c(1, sort(as.numeric(c(ind,(ind+1)))), length(d)))
###    resno.val <- pdb$atom[pdb$calpha,"resno"][resno.ind]
###    resno.val <- matrix(as.numeric(resno.val),nrow=2)
###
###    vec <- rep(blank, nrow(pdb$atom))
###    for(i in 1:(length(resno.val)/2)) {
###      sel.ind <- atom.select(pdb,
###                             resno=c(resno.val[1,i]:resno.val[2,i]),
###                             verbose=FALSE)
###      ##rep(LETTERS[i], length(sel.ind$atom))
###      vec[sel.ind$atom]=LETTERS[i]
###    }
###    return(vec)
###  }
###}
###
### --- See bio3d::chain.pdb()
###chain.pdb <- function(pdb, ca.dist=4, blank="X", rtn.vec=TRUE) {
###  ##- Find possible chian breaks
###  ##  i.e. Concetive Caplpa's that are further than 'ca.dist' apart,
###  ##        print basic chain info and rtn a vector of chain ids
###  ##        consisting of the 26 upper-case letters of the Roman
###  ##        alphabet
###  ##
###  ## chn <- chain.pdb(pdb)
###  ## pdb$atom[,"chain"] <- chain.pdb(pdb)
###  ##
###
###  ## Distance between concetive C-alphas
###  ca <- atom.select(pdb, "calpha", verbose=FALSE)
###  xyz <- matrix(pdb$xyz[ca$xyz], nrow=3)
###  d <- sqrt( rowSums( apply(xyz , 1, diff)^2 )/3 )
###
###  ## Chain break distance check
###  ind <- which(d > ca.dist)
###  len <- diff( c(1,ind,length(d)) )
###
###  cat(paste("Found",length(ind), " possible chain breaks\n"))
###  cat(paste("After resno(s)",
###        paste( pdb$atom[ca$atom,"resno"][(ind)], collapse=", " ),"\n" ))
###  cat(paste("chain length(s)",
###        paste(len+1, collapse=", " ),"\n" ))
###
###  ## Make a chain id vector
###  if(rtn.vec) {
###    resno.ind <- as.numeric(c(1, sort(as.numeric(c(ind,(ind+1)))), (length(d)+1)))
###    resno.val <- pdb$atom[ca$atom,"resno"][resno.ind]
###    resno.val <- matrix(as.numeric(resno.val),nrow=2)
###
###    vec <- rep(blank, nrow(pdb$atom))
###    for(i in 1:(length(resno.val)/2)) {
###      sel.ind <- atom.select(pdb,
###                             resno=c(resno.val[1,i]:resno.val[2,i]),
###                             verbose=FALSE)
###      vec[sel.ind$atom]=LETTERS[i]
###    }
###    return(vec)
###  }
###}
###
  
## pdbaln.R ==> see bio3d()


alitrim <- function(aln, cols=NULL, rows=NULL, atom=TRUE) {
  ##- Trim cols from alignment data structure
  ##  l <- alitrim(pdbs, cols=c(1:351))
  if(!atom) {
    stop("Only atom col inds supported")
  }
  naln <- NULL
  if(!is.null(cols)) {
    naln$id <- aln$id
    naln$b  <- aln$b[,-cols]
    naln$ali <- aln$ali[,-cols]
    naln$resno <- aln$resno[,-cols]
    naln$chain <- aln$chain[,-cols]
    naln$xyz <- aln$xyz[,-atom2xyz(cols)]
  } else { naln <- aln }
  if(!is.null(rows)) {
    naln$id <- naln$id[-rows]
    naln$b  <- naln$b[-rows,]
    naln$ali <- naln$ali[-rows,]
    naln$resno <- naln$resno[-rows,]
    naln$chain <- naln$chain[-rows,]
    naln$xyz <- naln$xyz[-rows,]
  }

  return(naln)
}


### --- See bio3d::seq2aln()
seq2aln.old <- function(seq, aln, aln.ind=NULL, id="seq") {
  ##- Add a sequence 'seq' to an existing alignment 'aln'
  ##  aln.ind is the row indice of the closest guy in aln to seq
  ##  l <- seq2aln(c(seq$ali), aln, which(aln$id %in% "d1i6ia_") )

  if(is.null(aln.ind)) {
    ## no aln.ind provided so lets try one from each aln group
    stop("Need to provide an 'aln.ind' as I havent botherd doing this yet")
  }

  ##- Align seq to masked template from alignment
  tmp.msk <- aln$ali[aln.ind,]
  tmp.msk[is.gap(tmp.msk)] <- "X"
  seq2tmp <- seqaln.pair(seqbind(seq, tmp.msk), id=c("seq","tmp"))


  ##- Insert gaps to adjust alignment
  ins <- which(is.gap( seq2tmp$ali[2,] ))
  if( length(ins)==0 ) {
    ntmp <- aln$ali
  } else {
    ins <- which(is.gap( seq2tmp$ali[2,] ))
    ntmp <- matrix(".", nrow=nrow(aln$ali), ncol=(ncol(aln$ali)+length(ins)))
    ntmp[,-ins] <- aln$ali
  }

  ## Add seq to top of adjusted alignment
  naln <- seqbind( seq2tmp$ali[1,], ntmp)
  rownames(naln) <- c(id, aln$id)
  return( list(id=c(id, aln$id), ali=naln) )
}


"write.pir" <- function (pdbs=NULL, ##ids=NULL, seqs=alignment$ali,
                         file, append = FALSE) {

#  if (is.null(seqs))
#    stop("write.fasta: please provide a 'seqs' or 'alignment' input object")

#  if (!is.null(alignment)) {
#    if (is.null(alignment$id) | is.null(alignment$ali)) {
#      stop("write.fasta: 'alignment' should be a list with '$id' and '$ali'components")
#    }
#    if (is.null(ids)) {
#      ids=alignment$id
#    }
#  } else {
#    if (is.null(ids)) {
#      n.ids <- nrow(seqs)
#      if(is.null(n.ids)) { n.ids=1 }
#      ids=seq( 1, length=n.ids )
#    }
#  }

  ## RESIDUE NUMBERING MUST BE CONSECITIVE! (bounds will have multiple rows)

  line2 <- rep("sequence::::::::0.00: 0.00", nrow(pdbs$resno))
  stru <- which(rowSums(is.na(pdbs$resno)) != ncol(pdbs$resno))
  ## structure:1ii6:1:A:347:A:bio3d:::
  for(i in 1:length(stru)) {
    rnum <- bounds( as.numeric(pdbs$resno[stru[i],]) )
    line2[stru[i]] <- paste("structure:", pdbs$id[ stru[i] ],
                            ":",rnum[,"start"],
                            ":",unique(na.omit(pdbs$chain[ stru[i] ,])),
                            ":", rnum[,"end"],":",
                            unique(na.omit(pdbs$chain[ stru[i] ,])),
                            ":bio3d:::",sep="")
  }

  if (!append)
    file.remove(file, showWarnings = FALSE)
  nseqs <- length(pdbs$id)
  if (nseqs == 1) {
    cat(">P1;", pdbs$id, "\n", line2, "\n", pdbs$ali, "*\n", file = file,
        append = TRUE, sep = "")
  }
  else {
    for (i in 1:nseqs) {
      cat(">P1;", pdbs$id[i], "\n", line2[i], "\n", pdbs$ali[i,],
          "*\n", file = file, append = TRUE, sep = "")
    }
  }
}

##======================================##
## All atom contact Functions
##======================================##

### --- See bio3d::cmap()
###ndm <- function(xyz, grpby=NULL, scut=NULL) {
###  ##-- New distance matrix function with 'grpby' option
###  ##  ndm(pdb$xyz, grpby=pdb$atom[,"resno"], scut=3)
###
###  
###  ##- Full Distance matrix (could use 'dm' or 'dist.xyz')
###  dmat <- as.matrix(dist(matrix(xyz, ncol = 3, byrow = TRUE)))
###  
###  if(is.null(grpby)) {
###    ##- Mask lower.tri  
###    dmat[lower.tri(dmat)] = NA
###    ##- Mask concetive residues
###    if (!is.null(scut))
###      dmat[diag.ind(dmat, n = scut)] = NA
###
###    return(dmat)
###    
###  } else {
###    ##- Group by concetive numbers in 'grpby'
###    if( length(xyz) != (length(grpby)*3) )
###      stop("dimension miss-match in 'xyz' and 'grpby', check lengths")
###
###    ##- Bounds of 'grpby' numbers
###    inds <- bounds(grpby, dup=TRUE)
###    nres <- nrow(inds)
###    
###    ##- Per-residue matrix
###    m <- matrix(, ncol=nres, nrow=nres)
###    ij <- pairwise(nres)
###    
###    ##  Ignore concetive residues
###    if (!is.null(scut))
###      ij <- ij[ij[,2]-ij[,1] > (scut-1),]
###    
###    ##- Min per residue
###    for(k in 1 : nrow(ij) ) {
###      m[ij[k,1],ij[k,2]] <- min( dmat[ (inds[ij[k,1],"start"]:inds[ij[k,1],"end"]),
###                                      (inds[ij[k,2],"start"]:inds[ij[k,2],"end"])],
###                                na.rm=TRUE )
###    }
###    return(m)
###  }
###}
###
###
###ncmap <- function(xyz, grpby,  dcut=4, scut=3) {
###
###  ## Distance matrix (all-atom)
###  dmat <- ndm( xyz, grpby, scut)
###  ## Contact map
###  return(matrix(as.numeric(dmat < dcut),
###                ncol = ncol(dmat),
###                nrow = nrow(dmat)))
###
###}
###






write.sge <- function(prefix = "r",
                      runtime = "72:00:00",
                      ncpu = 72,
                      iter = 1:10,
                      cluster="flux") {

  ## Write out a series of PBS or SGE
  ## shell scripts for running AMBER on flux or chemcca31
  ## > write.sge(prefix="r4q21")
  ## > write.sge("amd", cluster="amd")


  if (nchar(prefix) > 8) {
    warning(paste("Your filename 'prefix' is over 8 characters in length,\n\t",
                  "consider shortning so iteration digits are visible with qstat"))
  }

  if(cluster=="flux") {
    ## Run pmemd cMD on FLUX
    pbsfiles <- paste(prefix,".",sprintf("%02.0f", iter),".pbs",sep="")
    submitfile <- paste("submit_",prefix,".sh",sep="")
    
    for(i in 1:length(iter)) {

      cat(paste("#!/bin/bash
#PBS -S /bin/bash
#PBS -A bjgrant_flux
###PBS -N test
#PBS -q flux
#PBS -M bjgrant@umich.edu
#PBS -m abe
#PBS -j oe
#PBS -o out.$PBS_JOBNAME.log
#PBS -e err.$PBS_JOBNAME.log
#PBS -V
#PBS -l procs=",ncpu,",qos=bjgrant_flux,walltime=",runtime,"

echo Running job name $PBS_JOBNAME with ID $PBS_JOBID on host $PBS_O_HOST
echo Nodes for run:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cores
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Time is `date`
echo Directory is `pwd`

module delete openmpi
module load openmpi/1.4.2-intel
module load amber

##prv=07
##cur=08
prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

MPIRUN=/home/software/rhel5/openmpi-1.4.2/intel-11.0/bin/mpirun
AMBER=/home/software/rhel5/amber/11/exe/pmemd.MPI

AMBER_ARGS=\"-O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst -o dyna.$cur.out  -r dyna.$cur.rst -x dyna.$cur.traj.nc -inf dyna.$cur.inf\"

echo Running cmd: $MPIRUN -np $NPROCS  $AMBER $AMBER_ARGS

$MPIRUN -np $NPROCS $AMBER $AMBER_ARGS

echo Finished at time: `date`\n", sep=""),
        file=pbsfiles[i])
    }

    ##-- Write a master submission shell script for dependent jobs
    ## pbsfiles <- paste("r",".",sprintf("%02.0f", iter),".pbs",sep="")
    pbsids <- paste("r", iter, sep="") 
    head <- paste("#!/bin/bash \n## qsub -W depend=afterok:6664397 test1.pbs\n\n",
                  pbsids[1],"=`qsub ",pbsfiles[1],"`\necho $r1\n", sep="")

    middle <- paste(paste(pbsids[-1],"=`qsub -W depend=afterok:$",
                    pbsids[-length(pbsids)]," ",pbsfiles[-1],"`\necho $",
                    pbsids[-1], "\n", sep=""), collapse="", sep="")
    
    cat(head, middle, file=submitfile)
  }
    

#######

  if(cluster=="amd") {
    ## Run pmemd aMD on FLUX
    pbsfiles <- paste(prefix,".",sprintf("%02.0f", iter),".pbs",sep="")
    submitfile <- paste("submit_",prefix,".sh",sep="")
    
    for(i in 1:length(iter)) {

      cat(paste("#!/bin/bash
#PBS -S /bin/bash
#PBS -A bjgrant_flux
###PBS -N test
#PBS -q flux
#PBS -M bjgrant@umich.edu
#PBS -m abe
#PBS -j oe
#PBS -o out.$PBS_JOBNAME.log
#PBS -e err.$PBS_JOBNAME.log
#PBS -V
#PBS -l procs=",ncpu,",qos=bjgrant_flux,walltime=",runtime,"

echo Running job name $PBS_JOBNAME with ID $PBS_JOBID on host $PBS_O_HOST
echo Nodes for run:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cores
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Time is `date`
echo Directory is `pwd`

module delete intel-comp/11.0
module load intel-comp/12.0
module delete openmpi
module load openmpi/1.4.3/intel/12.0
module add med
module add amber-bjg

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

MPIRUN=/home/software/rhel5/openmpi-1.4.3/intel-12.0/bin/mpirun
AMBER=/home/software/rhel5/med/amber-bjg/11.0/exe/pmemd.MPI

AMBER_ARGS=\"-O -i dyna_amd.sander -p sys_box.prmtop -c dyna.$prv.rst -o dyna.amd.$cur.out -r dyna.amd.$cur.rst -x dyna.amd.$cur.traj.nc -inf dyna.amd.$cur.inf\"

echo Running cmd: $MPIRUN -np $NPROCS  $AMBER $AMBER_ARGS

$MPIRUN -np $NPROCS $AMBER $AMBER_ARGS

echo Finished at time: `date`\n", sep=""),
        file=pbsfiles[i])
    }

    ##-- Write a master submission shell script for dependent jobs
    ## pbsfiles <- paste("r",".",sprintf("%02.0f", iter),".pbs",sep="")
    pbsids <- paste("r", iter, sep="") 
    head <- paste("#!/bin/bash \n## qsub -W depend=afterok:6664397 test1.pbs\n\n",
                  pbsids[1],"=`qsub ",pbsfiles[1],"`\necho $r1\n", sep="")

    middle <- paste(paste(pbsids[-1],"=`qsub -W depend=afterok:$",
                    pbsids[-length(pbsids)]," ",pbsfiles[-1],"`\necho $",
                    pbsids[-1], "\n", sep=""), collapse="", sep="")
    
    cat(head, middle, file=submitfile)

  } else {
      cat("Have not added the code for other clusters yet\n")
    }

  ## Running instructions
  cat(paste(" *  SCP files to Cluster ",
            cluster, "\n\t(including:\t", submitfile, " dyna.",
            sprintf("%02.0f", (iter[1]-1)),
            ".rst sys_box.prmtop and \n\t\t\tdyna_prod.sander or dyna_amd.sander)\n\t",
            "(possibly cp dyna_equil.rst  to dyna.00.rst)\n\t",
            "> cp dyna_equil.rst dyna.00.rst\n\t",
            "> scp dyna.00.rst sys_box.prmtop dyna_prod.sander dyna_amd.sander ",
            prefix,".*.pbs submit_",prefix,".sh bgrant@",cluster,":somepath/.\n\n",

            " *  Submit all dependent jobs with run script:\n\t> sh ",
            submitfile,"\n\n",

              " *  Or submit individual jobs with:\n\t> qsub ",
            paste(prefix,".",sprintf("%02.0f", iter[1]),
                  ".sge",sep=""), "\n\t> qsub -hold_jid <JOBID> ",
            paste(prefix,".",sprintf("%02.0f", iter[2]),
                  ".sge",sep=""),"\n\t ...etc...\n\n",sep=""))
}







write.sge.OLD <- function(prefix = "r",
                      runtime = "72:00:00",
                      ncpu = 72,
                      iter = 1:10,
                      cluster="flux") {

  ## Write out a series of Sun Grid Engine
  ## shell scripts for running AMBER on flux or chemcca31
  ## > write.sge(prefix="4q21")


  if (nchar(prefix) > 8) {
    warning(paste("Your filename 'prefix' is over 8 characters in length,\n\t",
                  "consider shortning so iteration digits are visible with qstat"))
  }

  if(cluster=="flux") {
    ## Run pmemd on FLUX
    pbsfiles <- paste(prefix,".",sprintf("%02.0f", iter),".pbs",sep="")
    submitfile <- paste("submit_",prefix,".sh",sep="")
    
    for(i in 1:length(iter)) {

      cat(paste("#!/bin/bash
#PBS -S /bin/bash
#PBS -A bjgrant_flux
###PBS -N test
#PBS -q flux
#PBS -M bjgrant@umich.edu
#PBS -m abe
#PBS -j oe
#PBS -o out.$PBS_JOBNAME.log
#PBS -e err.$PBS_JOBNAME.log
#PBS -V
#PBS -l procs=",ncpu,",qos=bjgrant_flux,walltime=",runtime,"

echo Running job name $PBS_JOBNAME with ID $PBS_JOBID on host $PBS_O_HOST
echo Nodes for run:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cores
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Time is `date`
echo Directory is `pwd`

module delete openmpi
module load openmpi/1.4.2-intel
module load amber

##prv=07
##cur=08
prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

MPIRUN=/home/software/rhel5/openmpi-1.4.2/intel-11.0/bin/mpirun
AMBER=/home/software/rhel5/amber/11/exe/pmemd.MPI

AMBER_ARGS=\"-O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst -o dyna.$cur.out  -r dyna.$cur.rst -x dyna.$cur.traj.nc -inf dyna.$cur.inf\"

echo Running cmd: $MPIRUN -np $NPROCS  $AMBER $AMBER_ARGS

$MPIRUN -np $NPROCS $AMBER $AMBER_ARGS

echo Finished at time: `date`\n", sep=""),
        file=pbsfiles[i])
    }
    ##-- Write a master submission shell script for dependent jobs
    head <- paste("#!/bin/bash
## qsub -W depend=afterok:6664397 test1.pbs
FIRST=`qsub ",pbsfiles[1],"`
echo $FIRST\n")
    
    middle <- paste(paste(letters[iter[-1]],"=`qsub -W depend=afterok:$FIRST ",pbsfiles[-1],"` 
echo $",letters[iter[-1]],"\n", sep=""), collapse="", sep="")
    cat(paste(head, middle), file=submitfile)

  }


  
#########
  if(cluster=="c2") {
    ## Run pmemd on new cluster2
    for(i in 1:length(iter)) {

      cat(paste("#!/bin/bash
#$ -cwd
#$ -l h_rt=",runtime,"
#$ -V
#$ -S /bin/bash
#$ -pe orte ",ncpu,"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

MPIRUN=/soft/linux/pkg/openmpi/bin/mpirun
AMBER=/soft/linux/pkg/amber10/exe/pmemd

AMBER_ARGS=\"-O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst -o dyna.$cur.out  -r dyna.$cur.rst -x dyna.$cur.traj.crd -inf dyna.$cur.inf\"

$MPIRUN -np $NSLOTS  $AMBER $AMBER_ARGS

echo Finished at time: `date`\n", sep=""),
        file=paste(prefix,".",sprintf("%02.0f", iter[i]),".sge",sep=""))
    }
  }


  if(cluster=="c2AMD") {
    ## Run AMD on new cluster2
    for(i in 1:length(iter)) {

      cat(paste("#!/bin/bash
#$ -cwd
#$ -l h_rt=",runtime,"
#$ -V
#$ -S /bin/bash
#$ -pe orte ",ncpu,"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

MPIRUN=/soft/linux/pkg/openmpi/bin/mpirun
##AMBER=/soft/linux/pkg/amber10/exe/pmemd
AMBER=/gpfs/cesar/amber10_cesar/exe/sander.MPI

AMBER_ARGS=\"-O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst -o dyna.$cur.out  -r dyna.$cur.rst -x dyna.$cur.traj.crd -inf dyna.$cur.inf\"

$MPIRUN -np $NSLOTS  $AMBER $AMBER_ARGS

echo Finished at time: `date`\n", sep=""),
        file=paste(prefix,".",sprintf("%02.0f", iter[i]),".sge",sep=""))
    }
  }



  
  if(cluster=="c31") {
    ## Run on pmemd chemcca31
    for(i in 1:length(iter)) {

      cat(paste("#!/bin/bash
#$ -cwd
#$ -l h_rt=",runtime,"
#$ -V
#$ -S /bin/bash
#$ -pe mpich ",ncpu,"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
cat $TMPDIR/machines

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

/soft/linux/pkg/openmpi-1.2.7-intel-mx2g/bin/mpirun -v -np $NSLOTS -machinefile $TMPDIR/machines /soft/linux/pkg/amber10/exe/pmemd -O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst  -o dyna.$cur.out  -r dyna.$cur.rst  -x dyna.$cur.traj.crd  -inf dyna.$cur.inf\n", sep=""),
        file=paste(prefix,".",sprintf("%02.0f", iter[i]),".sge",sep=""))
    }
  } else {
    
    if(cluster=="oolite") {
      ## Run on pmemd oolite
      for(i in 1:length(iter)) {
        cat(paste("#!/bin/tcsh
#$ -cwd
#$ -l h_rt=",runtime,"
#$ -e sge.err
#$ -o sge.out
#$ -pe mpi ",ncpu,"
#$ -S /bin/tcsh

setenv MPI_MAX_CLUSTER_SIZE 4
setenv P4_GLOBMEMSIZE 32000000
#$ -v P4_GLOBMEMSIZE
#$ -V

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
cat $TMPDIR/machines

echo This job has allocated $NSLOTS processors

set mpirun=/usr/local/topspin/mpi/mpich/bin/mpirun
set pmemd=/share/apps/amber9/exe.pmemd.ts/pmemd

set prv=",sprintf("%02.0f", (iter[i]-1)),"
set cur=",sprintf("%02.0f",  iter[i]),"

  $mpirun -machinefile $TMPDIR/machines -np $NSLOTS $pmemd -O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst  -o dyna.$cur.out  -r dyna.$cur.rst  -x dyna.$cur.traj.crd  -inf dyna.$cur.inf\n", sep=""),
        file=paste(prefix,".",sprintf("%02.0f", iter[i]),".sge",sep=""))
      }

    }
############
    ## source("~/tmpwork/Rpackage/new_funs/model.R")
    ## write.sge("a3_apo", ncpu=32, iter=1:20, cluster="ctbp1m")
    
    if(cluster=="ctbp1m") {
      ## Run on pmemd ctbp1m
      for(i in 1:length(iter)) {

        cat(paste("#!/bin/bash
#$ -cwd
#$ -l h_rt=",runtime,"
#$ -V
#$ -S /bin/bash
#$ -pe mpich ",ncpu,"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
cat $TMPDIR/machines

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

/soft/linux/openmpi-1.2.7-intel-mx2g/bin/mpirun -v -np $NSLOTS -machinefile $TMPDIR/machines /soft/linux/amber10/exe/pmemd -O -i dyna_prod.sander -p sys_box.prmtop -c dyna.$prv.rst  -o dyna.$cur.out  -r dyna.$cur.rst  -x dyna.$cur.traj.crd  -inf dyna.$cur.inf\n", sep=""),
            file=paste(prefix,".",sprintf("%02.0f", iter[i]),".sge",sep=""))
      }
    }
###########    
    if (cluster=="amd") {
      ## Run aMD on chemcca31
      for(i in 1:length(iter)) {
        
        cat(paste("#!/bin/bash
#$ -cwd
#$ -l h_rt=",runtime,"
#$ -V
#$ -S /bin/bash
#$ -pe mpich ",ncpu,"

unsetenv LD_LIBRARY_PATH
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
cat $TMPDIR/machines

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

/opt/mpich/myrinet/intel/bin/mpirun -v -np $NSLOTS -machinefile $TMPDIR/machines /home/dhamelbe/amber8/exe.torsion.whole.fix/sander -O -i production.in -p sys_box.prmtop -c production.$prv.restrt -o production.$cur.out -r production.$cur.restrt -x production.$cur.mdcrd -inf mdinfo.$cur.inf\n", sep=""),
            file=paste(prefix,".",sprintf("%02.0f", iter[i]),".sge",sep=""))
      }
    } else {
      cat("Have not added the code for other clusters yet\n")
    }
  }

  ## Running instructions
  cat(paste(" *  SCP files to ",
            cluster, "\n\t(including dyna.",
            sprintf("%02.0f", (iter[1]-1)),
            ".rst sys_box.prmtop and dyna_prod.sander)\n\t",
            "(possibly cp dyna_equil.rst  to dyna.00.rst)\n\t",
            "> cp dyna_equil.rst dyna.00.rst\n\t",
            "> scp dyna.00.rst sys_box.prmtop dyna_prod.sander ",
            prefix,".*.pbs submit_",prefix,".sh bgrant@",cluster,":somepath/.\n\n",
            " *  Submit jobs with:\n\t> qsub ",
            paste(prefix,".",sprintf("%02.0f", iter[1]),
                  ".sge",sep=""), "\n\t> qsub -hold_jid <JOBID> ",
            paste(prefix,".",sprintf("%02.0f", iter[2]),
                  ".sge",sep=""),"\n\t ...etc...\n",sep=""))
}





##---- some PCA utility functions ----##
## z.coord: a numeric vector or row-wise matrix of z scores (i.e. centered and rotated data projected onto PCs) for conversion to xyz coordinates
##
## xyz.coord: a numeric vector or row-wise matrix of xyz coordinates for projection onto PCs 

## see pca.z2xyz and pca.xyz2z

interp <- function(start, end, nbeads=25) {
  ## interpolate between two vectors
  incby <- (end-start)/nbeads
  curr <- start; store <- curr
  for (i in 1:nbeads) {
    curr <- curr+incby
    store <- rbind(store, curr)
  }
  return(store)
}


### ---- See bio3d::pdbfit()
#####-- Define a new functions (an easy fit.xyz wrapper)
###fit <- function(pdbs, inds, outpath=NULL, prefix = "", pdbext = "") {
###  full <- ifelse(is.null(outpath), FALSE, TRUE)
###  return( fit.xyz( fixed=pdbs$xyz[1,], mobile=pdbs,
###                  fixed.inds  = inds, mobile.inds = inds,
###                  prefix = prefix, pdbext = pdbext,
###                  outpath = outpath, full.pdbs = full,het=TRUE ))
###}
###
###
###fit2 <- function(pdbs, inds=NULL, outpath=NULL, ...) {
###  full <- ifelse(is.null(outpath), FALSE, TRUE)
###  if(is.null(inds)) {
###    inds <- gap.inspect(pdbs$xyz)$f.inds
###  }
###  return( fit.xyz( fixed=pdbs$xyz[1,], mobile=pdbs,
###                  fixed.inds  = inds, mobile.inds = inds,
###                  outpath = outpath, full.pdbs = full, ... ))
###}
###

### --- See bio3d::ide.filter()
ide.group <- function(aln=NULL, ide=NULL, cutoff=0.6) {

  ##
  ## Return the indices of the largest group of sequences
  ## that have identity values above a particular 'cutoff'
  ##  (see ide.filter for different functionality)
  
  if(is.null(ide)) {
    if(is.null(aln)) 
      stop("Must provide either an alignment 'aln' or identity matrix 'ide'")
    ide  <- seqidentity(aln)
  }
  tree <- hclust( as.dist(1-ide) )
  grps <- cutree(tree, k = NULL, h = (1-cutoff))
  freq <- table(grps)
  if( max(freq)==1 ) {
    warning("No pair of sequences are above cutoff\n... returing a single index of the most representative sequence")
    return( which.max(colSums(ide)) )
  } else {
    ord  <- rev( order(freq) )
    return( which(names(freq[ord])[1]==grps) )
  }
}




##plot.blast <- function(b, cutoff=110) {
##  ##- Plot alignment stats
##  par(mfcol=c(4,1), mar=c(4, 4, 1, 2), cex.lab=1.5)
##  plot(b$mlog.evalue, xlab="Hit No", ylab="-log(Evalue)")
##  if(!is.null(cutoff))
##    abline(h=cutoff, col="red", lty=3)
##  plot(b$bitscore, xlab="Hit No", ylab="Bitscore")
##  plot(b$hit.tbl[,"identity"], xlab="Hit No", ylab="Identity")
##  plot(b$hit.tbl[,"alignmentlength"], xlab="Hit No", ylab="Length")
##
##  i <- b$mlog.evalue >= cutoff
##  
##  out <- cbind(b$pdb.id[i], b$gi.id[i])
##  rownames(out) <- which(i)
##  colnames(out) <- c("pdb.id", "gi.id")
##  return(out)
##}



##get.pdb <- function(ids, path="./") {
##  
##  if(any(nchar(ids) != 4)) {
##    warning("ids should be standard 4 character PDB formart: trying first 4 char...")
##    ids <- substr(basename(ids),1,4)
##  }
##  ids <- unique(ids)
##  
##  pdb.files <- paste(ids, ".pdb", sep="")
##  get.files <- file.path("http://www.rcsb.org/pdb/files", pdb.files)
##  put.files <- file.path( path, pdb.files)
##
##  dir.create(path)
##  rtn <- rep(NA, length(pdb.files))
##  for(k in 1:length(pdb.files)) {
##    rtn[k] <- download.file( get.files[k], put.files[k] )
##  }
##
##  names(rtn) <- paste( path, "/", ids, ".pdb", sep="")
##  if(any(rtn==1)) {
##    warning("Some files could not be downloaded, check returned value")
##    return(rtn)
##  } else {
##    return(names(rtn))
##  }
##}

##get.uniprot <- function(ids, path="./", onefile=FALSE) {
##  
##  if(any(nchar(ids) != 6)) {
##    warning("ids should be standard 6 character SWISSPROT/UNIPROT formart: trying first 6 char...")
##    ids <- substr(basename(ids),1,6)
##  }
##  ids <- unique(ids)
##  
##  pdb.files <- paste(ids, ".fasta", sep="")
##  get.files <- file.path("http://www.uniprot.org/uniprot", pdb.files)
##  put.files <- file.path( path, pdb.files)
##  mode="w"
##  
##  if(onefile) {
##    mode="a"
##    put.files <- rep(file.path( path, "combined_seqs.fasta"), length(pdb.files))
##  }
##  dir.create(path)
##  rtn <- rep(NA, length(pdb.files))
##  for(k in 1:length(pdb.files)) {
##    rtn[k] <- download.file( get.files[k], put.files[k], mode=mode )
##  }
##
##  names(rtn) <- ids
##  rtn
##}


##get.nr <- function(gi, path="./", onefile=FALSE) {
##
##  ids <- unique(gi)
##  pdb.files <- paste(ids, ".fasta", sep="")
##  get.files <- paste("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=",
##                     ids, "&dopt=fasta&sendto=t", sep="")
##  put.files <- file.path( path, pdb.files)
##
##  mode="w"
##  if(onefile) {
##    mode="a"
##    put.files <- rep(file.path( path, "combined_nr_seqs.fasta"), length(pdb.files))
##  }
##  dir.create(path)
##  rtn <- rep(NA, length(pdb.files))
##  for(k in 1:length(pdb.files)) {
##    rtn[k] <- download.file( get.files[k], put.files[k], mode=mode )
##  }
##
##  names(rtn) <- ids
##  rtn
##}



##fetch.seqs <- function(gi, path="./", onefile=TRUE) {
##  ## download FASTA format sequences from the NR database via
##  ## there gi number (see also get.pdb, get.uniprot and get.nr)
##  
##  ids <- unique(gi)
##  pdb.files <- paste(ids, ".fasta", sep="")
##  get.files <- paste("http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=",
##                     ids, "&dopt=fasta&sendto=t", sep="")
##  put.files <- file.path( path, pdb.files)
##
##  mode="w"
##  if(onefile) {
##    mode="a"
##    put.files <- rep(file.path( path, "nr_seqs.fasta"), length(pdb.files))
##  }
##  dir.create(path)
##  rtn <- rep(NA, length(pdb.files))
##  for(k in 1:length(pdb.files)) {
##    rtn[k] <- download.file( get.files[k], put.files[k], mode=mode )
##  }
##
##  names(rtn) <- ids
##  if(all(!rtn)) {
##    if(onefile) {
##      return(read.fasta( file.path( path, "nr_seqs.fasta") ))
##    } else {
##      return(rtn)
##    }
##  } else {
##    warning("Not all downloads were sucesfull, see returned values")
##    return(rtn)
##  }
##}



##fetch.pdbs <- function(ids, path="./") {
##  
##  if(any(nchar(ids) != 4)) {
##    warning("ids should be standard 4 character PDB formart: trying first 4 char...")
##    ids <- substr(basename(ids),1,4)
##  }
##  ids <- unique(ids)
##  
##  pdb.files <- paste(ids, ".pdb", sep="")
##  get.files <- file.path("http://www.rcsb.org/pdb/files", pdb.files)
##  put.files <- file.path( path, pdb.files)
##
##  dir.create(path)
##  rtn <- rep(NA, length(pdb.files))
##  for(k in 1:length(pdb.files)) {
##    rtn[k] <- download.file( get.files[k], put.files[k] )
##  }
##
##  names(rtn) <- ids
##  rtn
##}




##split.pdb <- function(raw.pdbs, path="split_chain/") {
##
##  dir.create(path)
##  for(i in 1:length(raw.pdbs)) {
##    pdb <- read.pdb(raw.pdbs[i], het2atom=TRUE)
##    chains <- unique(pdb$atom[,"chain"])
##    if(length(chains) > 0) {
##      for(j in 1:length(chains)) {
##        if( !is.na(chains[j]) ) {
##          sel <- paste("//",chains[j],"/////")
##          new <- atom.select(pdb,sel)
##          new.pdb <- NULL
##          new.pdb$atom <- pdb$atom[new$atom,]
##          new.pdb$xyz <-  pdb$xyz[new$xyz]
##
##          new.name <- paste(path,
##                            substr(basename(raw.pdbs[i]),1,4),
##                            "_", chains[j], ".pdb", sep="")
##
##          write.pdb(new.pdb, file=new.name)
##        }
##      }
##    }
##  }
##}



### --- See bio3d::pdbaln()
###pdb.seq.aln <- function(pdb.files, alnfile="pdb_aln.fa", ...) {
###
###  raw <- NULL
###  for(i in 1:length(pdb.files)) {
###    pdb <- read.pdb(pdb.files[i], het2atom=TRUE)
###    seq <- aa321(pdb$atom[pdb$calpha,"resid"])
###    if (is.null(seq)) {
###      raw <- seqbind(raw, c("-","-"))
###    } else {   
###      raw <- seqbind(raw, seq)
###    }
###  }
###  return(seqaln(raw, id=pdb.files, file=alnfile, ...))
###}
###

### --- See bio3d::torsion.pdb()
rama.inds <- function( pdb, zone=c(33:40) ) {
  
  ##- Return xyz indices for PHI-PSI Ramachendran atoms 

  phi <- matrix(NA,ncol=4*3, nrow=length(zone))
  psi <- phi

  for(i in 1:length(zone)) {
    b <- zone[i]; a <- b-1; c <- b+1
    
    phi[i,1:3]   <- atom.select(pdb, paste("///", a ,"///C/"),verbose=FALSE)$xyz
    phi[i,4:6]   <- atom.select(pdb, paste("///", b ,"///N/"),verbose=FALSE)$xyz
    phi[i,7:9]   <- atom.select(pdb, paste("///", b ,"///CA/"),verbose=FALSE)$xyz
    phi[i,10:12] <- atom.select(pdb, paste("///", b ,"///C/"),verbose=FALSE)$xyz

    psi[i,1:3]   <- atom.select(pdb, paste("///", b ,"///N/"),verbose=FALSE)$xyz
    psi[i,4:6]   <- atom.select(pdb, paste("///", b ,"///CA/"),verbose=FALSE)$xyz
    psi[i,7:9]   <- atom.select(pdb, paste("///", b ,"///C/"),verbose=FALSE)$xyz
    psi[i,10:12] <- atom.select(pdb, paste("///", c ,"///N/"),verbose=FALSE)$xyz
  }  
  return( list(psi = psi, phi=phi) )
}


plot.rama <- function(phi, psi, col=densCols( cbind(phi, psi)), ...) {
  library("geneplotter")
  require("RColorBrewer")
  plot(phi, psi, col=col,
       xlim=c(-180,180), ylim=c(-180,180),pch=20,
       xlab=expression(phi),
       ylab=expression(psi), ...)
}


### --- See bio3d::seq2aln()
aln2aln <- function(aln.mv, aln.mv.ind=1, aln.rf, aln.rf.ind=1) {

  ## saln <- read.fasta("coords/saln.15.afasta")
  ## moved <- saln
  ## for (i in 1:length(rep.str)) {
  ##  galn  <- read.fasta(paste("coords/grp_",i,".fa",sep=""))
  ##  ## mask X to Z
  ##  galn$ali[galn$ali=="X"] <- "Z"
  ## 
  ##  moved <- aln2aln(aln.mv = galn,
  ##                   aln.mv.ind = which(basename(galn$id) %in% moved$id),
  ##                   aln.rf = moved,
  ##                   aln.rf.ind = which(moved$id %in% basename(galn$id)) )
  ## }
  ## write.fasta(moved, file="tmp.fa") 
  ## moved$ali[which(moved$ali == "X",arr.ind=T)]="-"
  ## write.fasta(moved, file="moved.fa")


  
  mask <- function(x, char="X") {
    if(is.list(x))
      x <- x$ali
    x[which(x == "-",arr.ind=T)]=char
    x[which(x == ".",arr.ind=T)]=char
    return(x)
  }


  rf <- mask(aln.rf)
  mv <- mask(aln.mv)

  rf.seq <- rf[aln.rf.ind,]
  mv.seq <- mv[aln.mv.ind,]
  
  rm.aln <- seqaln.pair(seqbind(mv.seq, rf.seq), id=c("mv","rf") )$ali

  ## Insure ncol consistancy
  n.col  <- max( ncol(rm.aln), length(mv.seq), length(rf.seq) )
  rm.aln <- seqbind(rm.aln, rep("-",n.col))

  rf.ind <- which(!is.gap(rm.aln["rf",]))
  mv.ind <- which(!is.gap(rm.aln["mv",]))


  blank <- matrix("-", ncol=n.col, nrow=sum(nrow(rf),nrow(mv)))
  blank[1:nrow(rf), rf.ind] <- rf
  blank[(nrow(rf)+1):(nrow(rf)+nrow(mv)), mv.ind] <- mv

  naln <- NULL
  naln$ali <- blank
  naln$id  <- c(aln.rf$id, aln.mv$id)
  return(naln)
}


### --- See bio3d::pdbseq()
##seq.pdb <- function(pdb, inds=NULL) {
##  ## b.inds <- atom.select(pdb, "//B////CA/")
##  ## seq.pdb(pdb, b.inds)
##
##  if(is.null(inds))
##    inds <- atom.select(pdb, "//////CA/")$atom
##
##  if(is.list(inds))
##    inds <- inds$atom
##  
##  return(aa321(pdb$atom[inds,"resid"]))
##}


read.propka <- function(file) {

  ##-- examine out.propka
  ## ~/software/pdb2pqr/pdb2pqr.py --with-ph=7 --chain --salt --ffout=amber --ff=amber 1BE9.pdb out
  ## pka <- read.propka("out.propka")
  ## wiki.tbl(pka$his)


  
  ## lineformat
  first <- c(1, 4,  8, 16, 25, 48, 64, 80)
  last  <- c(3, 7, 15, 24, 47, 63, 79, 95)

  split.string <- function(x) {
    x <- substring(x, first, last)
    x[nchar(x) == 0] <- as.character(NA)
    x
  }
  trim <- function (s) {
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-NA
    s
  }
  
  raw.lines <- readLines(file)
  
  type <- substring(raw.lines,1,6)
  blank <- type=="  "
  type <- type[!blank]; raw.lines <- raw.lines[!blank]
  col <- which(type=="RESIDU")
  start <- col+2
  end <- which(type=="SUMMAR")-2
  
  ##trim(split.string(raw.lines[start]))
  
  mat <- matrix(trim(sapply(raw.lines[start:end], split.string)),byrow=TRUE, ncol=length(first))
  colnames(mat) <- c("resid","resno","pka","location","desolv","hside","hback","coulombic")
  return(list(pka=mat, his=mat[mat[,"resid"]=="HIS",]))
}




write.crdbox <- function(x, trjfile = "R.mdcrd") {
  ##- Output an AMBER crdbox format FORMAT(10F8.3) trajectory file

  cat(paste("TITLE : Created by Bio3D with",ncol(x)/3,"atoms"),
      fill=80, file=trjfile)
  for(j in 1:nrow(traj)) {
    cat( sprintf("%8.3f",x[j,]), sep="", fill=80,
        file=trjfile, append=TRUE )
    cat("   0.000    0.000    0.000", sep="",
        fill=80, file=trjfile, append=TRUE )
  }
}




## Basic utility functions 'cat0' and 'paste0', which are 
## like 'cat' and 'paste' except that their default code{sep}
## value is '""' (i.e. no space)

### --- See paste0() [now in base R]
###paste0 <- function(..., sep = "") {
###  paste(..., sep = sep)
##}


cat0 <- function (..., sep = "") {
  cat(..., sep = sep)
}


odd <-  function(x) {
  x != as.integer(x/2) * 2
}

even <- function(x) {
  x == as.integer(x/2) * 2
}


less <- function(x) {
  txt <- tempfile()
  sink(txt)
  print(x)
  sink()
  system(paste("less", txt))
  unlink(txt)
}

## makes all letters except first in word lower case
##gsubfn("\\B.", tolower, "I WOULD LIKE A BANANA SPLIT", perl = TRUE)



g03.min <- function(pdbfile, charge, multiplicity=1, comfile="out.com") {

  ## Run Geometry optimization with Gaussian (generate input for g03)
  ##
  ## N.B. the pdbfile and total charge values NEED to be set below
  ##
  ##  You can get charge and  multiplicity values (and your input
  ##  pdbfile with hydrogens) from mastro (or some other such source).
  ##  See the and wizard setup option and Jaguar molecule tab in Mastro.
  ##  Could also check http://ligand-expo.rcsb.org 

  if(is.null(charge))
    stop("Need know the net charge, charge=")

  if(is.null(pdbfile))
    stop("Need an input pdbfile, pdbfile=")

      
  pdb <- read.pdb(pdbfile, het2atom=TRUE)

  atoms <- substr(pdb$atom[,"elety"],1,1)

  if(sum(atoms=="H")==0)
    stop("Looks like your molecule has _NO_HYDROGENS_???")
  
  x <- as.numeric(pdb$atom[,"x"])
  y <- as.numeric(pdb$atom[,"y"])
  z <- as.numeric(pdb$atom[,"z"])
  format <- "%1s        %10.6f      %10.6f      %10.6f"

  lines <- NULL
  for (i in 1:length(atoms)) {
    lines <- rbind(lines, sprintf(format, atoms[i],
                                x[i], y[i], z[i]) )
  }

  ##- TO Do. add charge determination section!!
  
  header <- paste("%Chk=Job1-gaussian.chk
%Mem=256MB
%NProc=1

#P hf/6-31G* Opt=(Tight,CalcFC) Freq SCF(Conver=8) Test

Gaussian optimization input, to be used by R.E.D.

",charge, multiplicity,"\n")

  ## Could use some other options, e.g. 
  ##P AM1 Opt=Tight GFInput GFPrint SCF(Conver=8) Test\n


  cat(header, file = comfile)
  cat(lines, file = comfile, sep = "\n", append = TRUE)
  cat("\n", file = comfile, append = TRUE)

  cat(paste("   To run Gaussian optimization use the output script (",
            comfile,") with the cmd:",
            "\n\t g03 < ",comfile," > min.log\n", sep=""))

}




read.xyz <- function(xyzfile) {
  raw <- as.matrix(read.fwf(xyzfile,
                            widths=c(2,15,14,15),
                            skip = 2, strip.white=TRUE,
                            col.names=c("elety","x","y","z")) )
  return(raw)
}


g03.epot <- function(xyzfile=NULL, pdbfile=NULL, charge, multiplicity=1,
                     comfile="epot.com") {

  if(is.null(charge))
    stop("Need know the net charge, charge=")
    
  if(is.null(xyzfile)) {
    if(is.null(pdbfile)) {
      stop("Need an input coordfile, either xyzfile= or pdbfile=")
    } else {
      xyz <- read.pdb(pdbfile, het2atom=TRUE)$atom
    }
  } else {
    xyz <- read.xyz(xyzfile)
  }

  elety <- xyz[,"elety"]
  x <- as.numeric(xyz[,"x"])
  y <- as.numeric(xyz[,"y"])
  z <- as.numeric(xyz[,"z"])

  format <- "%1s        %10.6f      %10.6f      %10.6f"

  lines <- NULL
  for (i in 1:length(elety)) {
    lines <- rbind(lines, sprintf(format, elety[i],
                                  x[i], y[i], z[i]) )
  }


  header <- paste("%Mem=64MB
%NProc=1\n
#P HF/6-31G* SCF(Conver=6) NoSymm Test
   Pop=mk IOp(6/33=2) GFInput GFPrint\n
Molecular Electrostatic Potential min.xyz2epot.com\n
",charge, multiplicity,"\n")

  
  ## Could use some other options, e.g. 
  ##P AM1 Opt=Tight GFInput GFPrint SCF(Conver=8) Test\n


  cat(header, file = comfile)
  cat(lines, file = comfile, sep = "\n", append = TRUE)
  cat("\n", file = comfile, append = TRUE)

  cat(paste("   To run Gaussian potential calculation use the output script (",
            comfile,") with the cmd:",
            "\n\t g03 < ",comfile," > epot.log\n", sep=""))

}



elety2num <- function(atoms) {

  ele <- c("H", "He", "Li", "Be", "B", "C", "N", "O", "F",
           "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
           "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
           "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
           "As", "Se", "Br")

  convert <- function(x) {
    if (all(x != ele)) {
      return("X")
    } else {
      return(which(x == ele))
    }
  }
  return(as.vector(unlist(sapply(atoms, convert))))
}



xyz2resp <- function(xyzfile, charge, rest=NULL,
                     respfile="resp_input.out",
                     ioutopt=0, iqopt=1, nmol=1,
                     ihfree=1, irstrnt=1, qwt=0.0005) {

  
  if(is.null(charge))
    stop("Need know the net charge, charge=")
    
  if(is.null(xyzfile)) 
    stop("Need an input coordfile, xyzfile=")
    
  raw <- read.xyz(xyzfile)
  charge <- as.character(charge) ## <--- !****! --->
  natom  <- as.character(nrow(raw))
  nums   <- as.character(elety2num(raw[,"elety"]))

  if(is.null(rest)) ## <--- !****! --->
    rest   <- as.character(rep(0,nrow(raw))) 


  format <- "%5s%5s"
  
  lines <- sprintf(format,charge,natom)
  for (i in 1:nrow(raw)) {
    lines <- rbind(lines, sprintf(format, nums[i], rest[i]) )
  }
  
  header <- paste("  title\n",
                  " &cntrl\n",
                  "  ioutopt=", ioutopt,
                  ", iqopt=", iqopt,
                  ", nmol=", nmol,
                  ", ihfree=", ihfree,
                  ", irstrnt=",irstrnt,
                  ", qwt=", qwt,
                  "\n &end\n  1.0\n title\n", sep="")

  cat(header, file = respfile)
  cat(lines, file = respfile, sep = "\n", append = TRUE)
  cat("\n", file = respfile, append = TRUE)
}



### --- See bio3d::blast.pdb()
`blast.n` <-
function(seq, database="nr") {

  ## Run NCBI blastn on a given 'seq' sequence against a given 'database'
  if(!is.vector(seq)) {
    stop("Input 'seq' should be a single sequence as a single or multi element character vector")
  }
  seq <- paste(seq, collapse="")
  
  if( !(database %in% c("pdb", "nr", "swissprot")) )
    stop("Option database should be one of pdb, nr or swissprot")

  ##- Submit
  urlput <- paste("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?CMD=Put&DATABASE=",
                  database,"&HITLIST_SIZE=20000&PROGRAM=blastn&CLIENT=web&QUERY=",
                  paste(seq,collapse=""),
                  sep="")


  txt <- scan(urlput, what="raw", sep="\n", quiet=TRUE)
  rid <- sub("^.*RID = " ,"",txt[ grep("RID =",txt) ])

  cat(paste(" Searching ... please wait (updates every 5 seconds) RID =",rid,"\n "))

  
  ##- Retrive results via RID code (note 'Sys.sleep()')
  urlget <- paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get",
                  "&FORMAT_OBJECT=Alignment",
                  "&ALIGNMENT_VIEW=Tabular",
                  "&RESULTS_FILE=on",
                  "&FORMAT_TYPE=CSV",
                  "&RID=",rid, sep="")
  

  raw  <- read.csv(urlget,
                   header = FALSE, sep = ",", quote="\"", dec=".",
                   fill = TRUE, comment.char="")

  
  ## Check for job completion (retrive html or cvs?)
  html <- 1
  while(length(html) == 1) {
    cat("."); Sys.sleep(5)
    
    raw  <- read.csv(urlget,
                     header = FALSE, sep = ",", quote="\"", dec=".",
                     fill = TRUE, comment.char="")

    html <- grep("DOCTYPE", raw[1,])
  }

  
  colnames(raw) <- c("queryid", "subjectids", "identity", ##"positives",
                     "alignmentlength", "mismatches", "gapopens",
                     "q.start", "q.end", "s.start", "s.end",
                     "evalue", "bitscore")
  
  ## expand 'raw' for each hit in 'subjectids' (i.e. split on ";")
  rawm <- as.matrix(raw)
  
  eachsubject <- strsplit(rawm[,"subjectids"],";")
  subjectids  <- unlist(eachsubject)
  n.subjects  <- sapply(eachsubject, length)
  
  rawm <- apply(rawm, 2, rep, times=n.subjects)
  rawm[,"subjectids"] <- subjectids
  
  ## parse ids
  all.ids <- strsplit(subjectids, "\\|")
  gi.id  <- sapply(all.ids, '[', 2)
  pdb.id <- paste(sapply(all.ids, '[', 4),"_",sapply(all.ids, '[', 5),sep="")

  ## N.B. hack: zero evalues to arbirtrly high value!!
  mlog.evalue <- -log(as.numeric(rawm[,"evalue"]))
  mlog.evalue[is.infinite(mlog.evalue)] <- -log(1e-308)

  
  cat(paste("\n Reporting",length(pdb.id),"hits\n"))
  return(  list(bitscore=  as.numeric(rawm[,"bitscore"]),
                evalue =  as.numeric(rawm[,"evalue"]),
                mlog.evalue = mlog.evalue,
                gi.id = gi.id,
                pdb.id = pdb.id,
                hit.tbl = rawm,
                raw = raw) )
}



gb2fasta <- function(source.file,
                     destination.file = paste(basename(source.file),".fasta",sep="")) {
  input <- readLines(source.file)
  head <- input[1]
  head <- unlist(strsplit(head, split = " "))
  head <- head[nchar(head) > 0]
  seqname <- head[2]
  seqsize <- as.integer(head[3])
  outheader <- sprintf(">%s %d bp", seqname, seqsize)
  confile <- file(destination.file, open = "w")
  writeLines(outheader, confile)
  debut <- which(substring(input, 1, 6) == "ORIGIN") + 1
  if (length(debut) > 1) 
    stop("Multiple entries not yet implemented !")
  fin <- which(substring(input, 1, 2) == "//") - 1
  if (length(fin) > 1) 
    stop("Multiple entries not yet implemented !")
  input <- input[debut:fin]
  input <- sapply(input, function(x) {
    return(paste(substr(x, 11, 20), substr(x, 22, 31), substr(x, 
            33, 42), substr(x, 44, 53), substr(x, 55, 64), substr(x, 
            66, 75), sep = "", collapse = ""))
    })
  names(input) <- NULL
  writeLines(input, confile)
  close(confile)
  return(read.fasta(destination.file))
}


free2fasta <- function(infile, outfile) {
  
  ## Convert FREE (SCA) format to FASTA format
  ## fasta2free("PDZ_final.fasta", "poo")
  ## free2fasta("PDZ_final.free", "PDZ_final.fasta")

  raw <- read.table(infile, stringsAsFactors=FALSE)

  suppressWarnings( file.remove(outfile) )
  for(i in 1:nrow(raw)) {
    cat(paste(">",raw[i,1],sep=""), append=TRUE, sep ="\n", file=outfile)
    cat(raw[i,2], append=TRUE, sep ="\n",file=outfile)
  }
  
}


fasta2free <- function(infile, outfile) {

  ## Convert FATA format to FREE (SCA) format
  ## fasta2free("PDZ_final.fasta", "poo")
  ## free2fasta("PDZ_final.free", "PDZ_final.fasta")

  raw <- read.fasta(infile)

  suppressWarnings( file.remove(outfile) )
  for(i in 1:length(raw$id)) {
    space <- paste( rep(" ", (10 - nchar(raw$id[i]))), collapse="")
    cat(paste( raw$id[i], space, 
              paste(raw$ali[i,], collapse=""), sep=""),
        append=TRUE, sep ="\n",file=outfile)
  }
}



read.apbs <- function(f) {
  ## Read one or more APBS log files
  ## use regular UNIX * eild card for multi file reads e.g.
  ##  read.apbs("tub_mutant_pqrs/apbs_C.*.log")
  
  cmd <- paste("grep 'Global net ELEC energy ='",f)
  raw <- read.table( pipe(cmd) )
  ind <- which(raw[1,]=="energy")+2
  return( as.numeric(raw[1,ind]) )
  
##  if(length(grep("\\*",f))>0 ) {
##    raw <- c(as.matrix(read.table( pipe(cmd),
##                colClasses=c(rep('NULL',6), "numeric", "NULL")) ))
##  } else {
###    raw <- c(as.matrix( read.table( pipe(cmd),
###                colClasses=c(rep('NULL',5),"numeric", "NULL")) ))
##    raw <- c(as.matrix( read.table( pipe(cmd),
##                colClasses=c(rep('NULL',6),"numeric", "NULL")) ))
##  }
##  return(unique(raw))
}

"ulsq" <-function(xx,
                    yy,
                    xfit=rep(TRUE,length(xx)),
                    yfit=xfit,
                    verbose=FALSE) {
  
  ## Return the rotation and translation matrix ... (see tlsq)
  ##
  ## Coordinate superposition with the Kabsch algorithm
  ## from Acta Cryst (1978) A34 pp827-828 (to which equation no. refer)
  ## yy is the target (i.e. fixed)

  xx <- matrix(xx,nrow=3, ) 
  x <- matrix(xx[xfit],nrow=3, ) 
  y <- matrix(yy[yfit],nrow=3, )
  
  if(length(x) != length(y)) stop("dimension mismatch in x and y")

  ## mean positions 
  xbar <- apply(x,1,mean) ; ybar <- apply(y,1,mean)

  ## center both sets
  xx <- sweep(xx,1,xbar) # NB xx centred on xbar 
  x <- sweep(x,1,xbar) ; y <- sweep(y,1,ybar)



  
  ##irmsd <- sqrt(sum((x-y)^2)/dim(y)[2])
  ##cat("#irmsd= ",round(irmsd,6),"\n")

  # generate the 3x3 moment matrix: R (Equation 3)
  R <- y %*% t(x)

  # form R'R
  RR <- t(R) %*% R
  
  # diagonalize R'R
  prj <- eigen(RR)
  prj$values[prj$values < 0 & prj$values >= -1.0E-12]<-1.0E-12
  
  # form A
  A <- prj$vectors

  # make explicitly rh system
  #	A[,3] <- v3cross(A[,1],A[,2])
  # inline the cross-product function call.
  b<-A[,1]; c <- A[,2]
  A[1,3] <- (b[2] * c[3]) - (b[3] * c[2])
  A[2,3] <- (b[3] * c[1]) - (b[1] * c[3])
  A[3,3] <- (b[1] * c[2]) - (b[2] * c[1])

  # form B (==RA) (Equation 8)
  B <- R %*% A

  # normalize B
  # B <- sweep(B,2,sqrt(apply(B^2,2,sum)),"/")
  B <- sweep(B,2,sqrt(prj$values),"/")

  # make explicitly rh system
  #	B[,3] <- v3cross(B[,1],B[,2])
  # inline the cross-product function call.
  b<-B[,1]; c <- B[,2]
  B[1,3] <- (b[2] * c[3]) - (b[3] * c[2])
  B[2,3] <- (b[3] * c[1]) - (b[1] * c[3])
  B[3,3] <- (b[1] * c[2]) - (b[2] * c[1])

  # form U (==Ba) (Equation 7)
  # U is the rotation matrix
  U <-  B %*% t(A)

  ####====###
  return(list(U=U, ybar=ybar, xbar=xbar))
  ####====###
  
  ## here we apply transformation matrix to *all* elements of xx
  ## rotate xx (Uxx)
###  xx <- U %*% xx

  ## return xx centred on y
###  xx <- sweep(xx,1,ybar,"+")
###  as.vector(xx)
}

tlsq <- function(xyz, u) {
  ## transform coords bu rotation and translation from ulsq
  return( as.vector(sweep( (u$U %*%
                            sweep( matrix(xyz, nrow=3),1,u$xbar) ),
                          1, u$ybar,"+" )) )
}


`read.maciek` <-
function (file, maxlines=50000, multi=FALSE,
                      rm.insert=FALSE, rm.alt=TRUE, verbose=TRUE) {

  if(missing(file)) {
    stop("read.pqr: please specify a PQR 'file' for reading")
  }
  if(!is.numeric(maxlines)) {
    stop("read.pqr: 'maxlines' must be numeric")
  }
  if(!is.logical(multi)) {
    stop("read.pqr: 'multi' must be logical TRUE/FALSE")
  }
  
  # PDB FORMAT v2.0:    colpos,  datatype,    name,      description
  atom.format <- matrix(c(-6,     NA,          NA,       # (ATOM)
                          5,     'numeric',   "eleno",   # atom_no
                         -1,     NA,          NA,        # (blank)
                          3,     'character', "elety",   # atom_ty
                          1,     'character', "alt",     # alt_loc
                          4,     'character', "resid",   # res_na 
                          1,     'character', "chain",   # chain_id 
                          5,     'numeric',   "resno",   # res_no
                          1,     'character', "insert",  # ins_code
                         -3,     NA,           NA,       # (blank)
                          10,     'numeric',   "x",       # x
                          10,     'numeric',   "y",       # y
                          10,     'numeric',   "z",       # z
                          8,     'numeric',   "o",       # o  ### 6 for pdb
                          8,     'numeric',   "b",       # b  ### 6 for pdb
                         -6,     NA,           NA,       # (blank)
                          4,     'character', "segid"    # seg_id
                         ), ncol=3, byrow=TRUE,
                       dimnames = list(c(1:17), c("widths","what","name")) )


  ## split.string(raw.atom[1])
  ## raw.atom[1]
  split.string <- function(x) {
    # split a string 'x'
    x <- substring(x, first, last)
    x[nchar(x) == 0] <- as.character(NA)
    x
  }
  is.character0 <- function(x){length(x)==0 & is.character(x)}
  
  trim <- function (s) {
    # Remove leading and traling
    # spaces from character strings
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-NA
    s
  }


  # finds first and last (substr positions)
  widths <-  as.numeric(atom.format[,"widths"]) # fixed-width spec  
  drop.ind <- (widths < 0) # cols to ignore (i.e. -ve)
  widths <- abs(widths)    # absolute vales for later  
  st <- c(1, 1 + cumsum( widths ))
  first <- st[-length(st)][!drop.ind] # substr start
  last <- cumsum( widths )[!drop.ind] # substr end

  # read n lines of PDB file
  raw.lines  <- readLines(file, n = maxlines)
  type <- substring(raw.lines,1,6)

  # check number of END/ENDMDL records
  raw.end <- sort(c(which(type == "END"),
                    which(type == "ENDMDL")))
  
  if (length(raw.end) > 1) {
    print("PDB has multiple END/ENDMDL records")
    if (!multi) {
      print("multi=FALSE: taking first record only")
      raw.lines <- raw.lines[ (1:raw.end[1]) ]
      type <- type[ (1:raw.end[1]) ]
    } else {
      print("multi=TRUE: 'read.dcd' will be quicker!")
    }
  }
  if ( length(raw.end) !=1 ) {
    if (length(raw.lines) == maxlines) {
      # have not yet read all the file
      print("You may need to increase 'maxlines'")
      print("check you have all data in $atom")
    }
  }

  # split by record type
  raw.header <- raw.lines[type == "HEADER"]
  raw.seqres <- raw.lines[type == "SEQRES"]
  raw.helix  <- raw.lines[type == "HELIX "]
  raw.sheet  <- raw.lines[type == "SHEET "]
  raw.atom   <- raw.lines[type == "ATOM  "]
  het.atom   <- raw.lines[type == "HETATM"]
  # also look for "TER" records
  rm(raw.lines)
  
  if (verbose) {
    if (!is.character0(raw.header)) { cat(" ", raw.header, "\n") }
  }
  seqres <- unlist(strsplit( trim(substring(raw.seqres,19,80))," "))

  helix  <- list(start = as.numeric(substring(raw.helix,22,25)),
                 end   = as.numeric(substring(raw.helix,34,37)),
                 chain = trim(substring(raw.helix,20,20)),
                 type  = trim(substring(raw.helix,39,40)))

  sheet  <- list(start = as.numeric(substring(raw.sheet,23,26)),
                 end   = as.numeric(substring(raw.sheet,34,37)),
                 chain = trim(substring(raw.sheet,22,22)),
                 sense = trim(substring(raw.sheet,39,40)))
  
  # format ATOM records as a character matrix
  atom <- matrix(trim(sapply(raw.atom, split.string)), byrow=TRUE,
                 ncol=nrow(atom.format[ !drop.ind,]), 
                 dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )

  # Alt records with m[,"alt"] != NA
  if (rm.alt) {
    if ( sum( !is.na(atom[,"alt"]) ) > 0 ) {
      cat("   PDB has ALT records, taking A only, rm.alt=TRUE\n")
      alt.inds <- which( (atom[,"alt"] != "A") ) # take first alt only
      if(length(alt.inds)>0)
        atom <- atom[-alt.inds,]
    }
  }
  # Insert records with m[,"insert"] != NA
  if (rm.insert) {
    if ( sum( !is.na(atom[,"insert"]) ) > 0 ) {
      cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
      insert.inds <- which(!is.na(atom[,"insert"])) # rm insert positions
      atom <- atom[-insert.inds,]
    }
  }
  het <- matrix(trim(sapply(het.atom, split.string)), byrow=TRUE,
                ncol=nrow(atom.format[ !drop.ind,]), 
                dimnames = list(NULL, atom.format[ !drop.ind,"name"]) )

  output<-list(atom=atom,
               het=het,
               helix=helix,
               sheet=sheet,
               seqres=seqres,
               xyz=as.numeric(t(atom[,c("x","y","z")])),
               calpha = as.logical(atom[,"elety"]=="CA"))

  class(output) <- "pdb"
  return(output)
  
}


occur <- function(v) {
  ## The number of occurances of something in a vector at that point
  ## v <- c("HSPA1", "HSPA2", "HSPA8", "HSPA8", "HSPA2", "HSPA8")
  ## occur(v)
  ## see also make.unique(v) and table
  ##
  
  n <- max(table(v))
  ape <- rep(1, length(v))
  tmp <- v
  for(i in 1:n) {
    dup <- duplicated(tmp, incomparables=NA)
    ape[dup] = ape[dup] + 1
    tmp[!dup] = NA
  }
  names(ape) = v
  return(ape)
}

cat.pdb <- function(pdb1, pdb2) {
  ## Concatenate two PDB files
  ## n <- cat.pdb( read.pdb( "4q21" ), read.pdb( "5p21" ))
  ## write.pdb(n, file="mycatstructure.pdb")

  npdb <- pdb1
  npdb$atom <- rbind(npdb$atom, pdb2$atom)
  npdb$xyz <- c(npdb$xyz, pdb2$xyz)
  npdb$calpha <- c(npdb$calpha, pdb2$calpha)
  return(npdb)
}


###--- See bio3d::convert.pdb()
###"renumber.pdb" <-
###function(pdb, first.resno=1, first.eleno=1, bychain=FALSE) {
###
###  ## Renumber 
###  ori.num <- as.numeric(pdb$atom[,"resno"])
###  ##ori.res <- pdb$atom[,"resid"]
###
###  if(bychain) {
###    stop("No code here yet")
###    ##ch <- chain.pdb(pdb)
###    
###  } else {
###    ## concetive residue numbers
###    ##tbl <- table(ori.num)
###    ##new.nums <- first.resno:(first.resno+length(tbl)-1)
###    ##pdb$atom[,"resno"] <- rep(new.nums, tbl)
###
###    ## Everytime old resno changed
###    new.res.len <- rle(ori.num)$lengths
###    nres <- length(new.res.len)
###    new.num <- first.resno:(first.resno+nres-1)
###    new.num <- rep(new.num,  new.res.len)
###  }
###
###  ## Eleno
###  npdb <- pdb
###  npdb$atom[,"eleno"] <- seq(first.eleno, length=nrow(pdb$atom))
###  npdb$atom[,"resno"] <- new.num
###  
###  return(npdb)
###}
###


## Match vector via matching sequence to alignment
vec2seq <- function(vec, aa.seq, aln) {
  ## Match a vector of values 'vec' for a given 
  ## alignment 'aln' to a new seq 'aa.seq'
  Naln <- seq2aln(aa.seq, aln, 1)
  Nmatch <- Naln$ali[1,]
  if(length(Nmatch) == length(vec)) {
    Nvec <- vec[!is.gap(Nmatch)]
  } else {
    Nvec <- vec[!is.gap(Nmatch)]
    Nvec[is.na(Nvec)] = 0
    ##Ngap <- (Naln$ali[2,]==".")
    ##Nvec[Ngap] = 0
    ##print(Nmatch)
    warning("Lengths don't match")
  }
  return(Nvec)
}

### --- See bio3d::vec2resno()
###vec2resno <- function(vec, resno) {
###  ## replicate vec based on concetive
###  ## similar resno entries
###
###  if(class(resno)=="pdb")
###    resno <- pdb$atom[,"resno"]
###
###  res.len <- rle(resno)$lengths
###  if(length(vec) != length(res.len))
###    stop("Length miss-match of 'vec' and concetive 'resno'")
###
###  if( sum(res.len) != length(resno) )
###    stop("Replicated length Miss-match")
###  
###  return( rep(vec,  times=res.len))
###}
###


### --- See bio3d::vec2resno()
vec2dimer <- function(vec, aln, pdb, add.chain=TRUE) {
  ## Function to make a numeric vector from
  ## conservation data suitable for the Bfactor
  ## filed of a dimeric PDB file with different
  ## atom numbers

  if(add.chain) {
    ## Boundarys of chain A and B
    ch <- chain.pdb(pdb)
    pdb$atom[,"chain"] <- ch
  }
    
  ## Extract sequence
  Aseq <- seq.pdb(pdb, atom.select(pdb, chain="A", elety="CA"))
  Bseq <- seq.pdb(pdb, atom.select(pdb, chain="B", elety="CA"))

  ## Vector equivalence from adding sequence to alignment
  Avec <- vec2seq(vec=vec, aa.seq=Aseq, aln=aln)
  Bvec <- vec2seq(vec=vec, aa.seq=Bseq, aln=aln)

  ## Replicate values for each atom in all residues
  val <- vec2resno( c(Avec,Bvec), resno=pdb$atom[,"resno"])
  return(val)
}



### --- See bio3d::bwr.colors()
"bgr.colors" <-
function (n) {

  # Derived from the function colorpanel
  # by Gregory R. Warnes

  if(n<3)
    warning("not sensible to ask for less than 3 colors")
  
  odd = FALSE
  
  if (n != as.integer(n/2) *2) {
    n   <- n + 1
    odd = TRUE
  }
  low  <- col2rgb("blue")
  mid  <- col2rgb("gray")
  high <- col2rgb("red")
  
  lower <- floor(n/2)
  upper <- n - lower
  
  red   <- c(seq(low[1, 1], mid[1, 1], length = lower),
             seq(mid[1,1], high[1, 1], length = upper))/255
  
  green <- c(seq(low[3, 1], mid[3, 1], length = lower),
             seq(mid[3, 1], high[3, 1], length = upper))/255
  
  blue  <- c(seq(low[2, 1], mid[2, 1], length = lower),
             seq(mid[2, 1], high[2, 1], length = upper))/255
    
  if (odd) {
    red   <- red[-(lower + 1)]
    green <- green[-(lower + 1)]
    blue  <- blue[-(lower + 1)]
  }
  rgb(red, blue, green)
}



## improved list of objects
.ls.objects <- function(pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {

  napply <- function(names, fn) sapply(names, function(x)
                                       fn(get(x, pos = pos)))

  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
                      as.numeric(dim(x))[1:2]))

  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")

  obj.dim[vec, 1] <- napply(names, length)[vec]

  out <- data.frame(obj.type, obj.size, obj.dim)

  names(out) <- c("Type", "Size", "Rows", "Columns")

  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]

  if (head)
    out <- head(out, n)

  out
}

# shorthand

lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}




read.hmmer.tbl <- function(infile) {
  ##- Function to Read HMMER3 log files produced with
  ##  the --tblout option to hmmsearch

  col.vals <- c("target", "t.accession",
                "query", "q.accession",
                "evalue", "score", "bias",
                "D.evalue", "D.score", "D.bias",
                "exp", "reg", "clu",  "ov",
                "env", "dom", "rep", "inc",
                "description")

  ## Find width of variable first field from third line
  col.widths <- c(11,21,10, 10,7,6,10,7, 7,4,4,4,4,4,4,4,4, 200)
  raw <- readLines(infile, n=3)
  first.field <- (motif.find("- -", raw[3])[1]+1)
  col.widths <- c(first.field, col.widths)

  return(read.fwf(infile, widths=col.widths, col.names=col.vals))
}


### --- See bio3d::conserv()
###cons.aln <- function(aln, gap.col=0.6) {
###  ## Score residue conservation
###  ## Exclude positions with >= gap.col percent gaps
###  ##
###  ##aln <- read.fasta("~/work/eb1/pfam_aln2.fa")
###  ##con <- cons.aln(aln)
###  ##plot(con$sim, typ="l")
###  ##plot(con$ent.10, typ="h", col=con$col.ent.10)
###  
###  sim <- conserv(x=aln$ali, method="similarity", sub.matrix="bio3d")
###  ent.10 <- conserv(x=aln$ali, method="entropy10")
###  ent.22 <- conserv(x=aln$ali, method="entropy22")
###  ide <- conserv(x=aln$ali, method="identity")
###
###  ## excluding positions >=60 percent gaps)
###  H <- entropy(aln$ali)
###  gap60.inds <- which(apply(H$freq[21:22,],2,sum)>=gap.col) ##0.6)
###
###  ent.10[gap60.inds] = 0
###  ent.22[gap60.inds] = 0
###
###  col.sim <- rep("black", length(sim))
###  col.ent.10 <- col.sim
###  col.ent.22 <- col.sim
###  col.ide <- col.sim
###
###  col.sim[sim >= 0.7] = "red"
###  col.ent.10[ent.10 >= 0.7] = "red"
###  col.ent.22[ent.22 >= 0.7] = "red"
###  col.ide[ide >= 0.7] = "red"
###
###  return(list(sim=sim, ent.10=ent.10, ent.22=ent.22, ide=ide,
###              gap60.inds=gap60.inds, col.sim=col.sim,
###              col.ent.10=col.ent.10, col.ent.22=col.ent.22,
###              col.ide=col.ide) )
###}
###

pathInterp <- function(z, npoint=12) {
  ## Given PCA Z coords for two conformers (rows) return
  ## 'npoint' interpolated Z coords btwn these two
  ## then pass this to pcaz2xyz()
  ##  E.g.
  ##   pca.z2xyz( pathInterp(pc.xray$z[inds,1:3]), pc.xray)
  ##
  ## z=z.pth; npoint=12
  
  if(nrow(z)!=2)
    stop("Input should have 2 rows only: Start and End")
  ndims <- ncol(z)
  blank=matrix(NA, nrow=npoint, ncol=ndims)
  for(i in 1:ndims) { 
    blank[,i] <- seq(from=z[1,i], to=z[2,i], len=npoint)
  }
  return(blank)
}


##
## Description:  Functions for FTmap analysis of trajectorys and aligned PDBs
## Date:         Wed Feb  8 10:42:39 EST 2012 
## Author:       Barry Grant bjgrant@umich.edu
##

probe.trj <- function(pdbfiles, max1=TRUE) {
  ##
  ## Analyse FTMAP results for trajectory PDB files
  ##  i.e. pdbfiles that have the same number of resdues/atoms
  ##
  ##  For input files with different sequences and/or different
  ##  residue counts use diferent "probe.pdbs()" function
  ##
  ##  E.G.
  ##     pth <- "cMD_ftmap/"
  ##     pdbfiles <- list.files(pth, pattern="ftmap.pdb$",full.names=TRUE)
  ##     probes <- probe.trj(pdbfiles)
  ##     apply(probes,2,mean, na.rm=TRUE)
  ##
  
  pdb <- read.pdb( pdbfiles[1] )
  chains <- unique(pdb$atom[,"chain"])

  ca.inds <- atom.select(pdb, chain=chains[1], elety="CA")
  nres <- length(ca.inds$atom)
  nfiles <- length(pdbfiles)

  ##- Create some blank data structure for keeping results
  probe.count <- matrix(NA, ncol=nres, nrow=nfiles)

  ##- Loop through input processing each ftmap pdb file.
  for(j in 1:nfiles) {

    ## read an ftmap pdb
    pdb <- read.pdb(pdbfiles[j])

    ## What chains do we have
    chains <- unique(pdb$atom[,"chain"])

    ## Select Calpha chain 1 indices.              ### <- change per chains
    ca.inds <- atom.select(pdb, chain=chains[1], elety="CA") 

    ## Coords for protein
    ca.xyz <- pdb$xyz[ca.inds$xyz]

    ## Keep top 5 chains (probe clusters only)
    chn=length(chains)
    if(chn > 6) chn=6
  
    ## blank matrix to store min distance to probe atoms 
    dmin.allres <- matrix(NA, nrow=length(ca.inds$atom), ncol=(chn-1))

  
    for(i in 2:chn) {
      ## Probe chains
      probe.inds <- atom.select(pdb, chain=chains[i])
      probe.xyz <- pdb$xyz[probe.inds$xyz]
    
      d <- dist.xyz(ca.xyz, probe.xyz, all.pairs=TRUE)
      d.min <- apply(d, 1, min)
      dmin.allres[,(i-1)] <- d.min
    }

    ##rownames(dmin.allres) <- pdb$atom[ca.inds$atom,"resno"]
    colnames(dmin.allres) <- chains[2:chn]

    ##-- Store frequencey 'count' of chain probes per residue
    probe.count[j,] <- rowSums(dmin.allres <= 5)
  }
  colnames(probe.count) = pdb$atom[ca.inds$atom,"resno"]
  if(max1==TRUE) {
    probe.count[probe.count>1]=1
  }
  return(probe.count)
}




probe.val <- function(p.norm) {
  return( colSums(p.norm, na.rm=T)/nrow(p.norm) )
}



probe.pdbs <- function(pdbfiles, pdbs, max1=TRUE) {
  ##
  ##  Analyse FTMAP results for aligned PDB files
  ##   "pdbfiles" is the filenames of ftmap output PDBS
  ##   "pdbs" is a bio3d class structure alignment object
  ##     as returned from "read.fasta.pdb()" and "pdbaln()"
  ##  NOTE:
  ##    The order of "pdbfiles" and the rows in "pdbs" must
  ##    be the SAME!
  ##    E.G.
  ##       ids <- substr(basename(pdbs$id),1,6)
  ##       pth <- "raw_FTMap_results/"
  ##       pdbfiles <- paste(pth, ids,".fftmap.pdb",sep="")
  ##       out <- probe.pdbs(pdbfiles, pdbs)
  ##       p.all <- probe.val(out)
  ##
  
  probe.count <- matrix(NA, ncol=ncol(pdbs$ali), nrow=nrow(pdbs$ali))

  for(j in 1:length(pdbfiles)) {
    pdb <- read.pdb(pdbfiles[j])

    ## What chains do we have
    chains <- unique(pdb$atom[,"chain"])

    ## Select Calpha chain 1 indices abd coords
    ca.inds <- atom.select(pdb, chain=chains[1], elety="CA") 
    ca.xyz <- pdb$xyz[ca.inds$xyz]

    ## Keep only top 5 clusters
    chn=length(chains)
    if(chn > 6) chn=6
  
    ## blank matrix to store distance to 
    dmin.allres <- matrix(NA, nrow=length(ca.inds$atom), ncol=(chn-1)) 
  
    for(i in 2:chn) {                  
      ## Probe chains
      probe.inds <- atom.select(pdb, chain=chains[i])
      probe.xyz <- pdb$xyz[probe.inds$xyz]
    
      d <- dist.xyz(ca.xyz, probe.xyz, all.pairs=TRUE)
      d.min <- apply(d, 1, min)
      dmin.allres[,(i-1)] <- d.min
    }

    ##rownames(dmin.allres) <- pdb$atom[ca.inds$atom,"resno"]
    colnames(dmin.allres) <- chains[2:chn]

    ##-- Store frequencey 'count' of chain probes per residue
    probe.count[j,!is.na(pdbs$resno[j,]) ] <- rowSums(dmin.allres <= 5)
  }
  ## Norm to max 1 value
  if(max1==TRUE) {
    probe.count[probe.count>1]=1
  }
  return(probe.count)
}


ttest.probe <- function(Aprobe.norm, Bprobe.norm, cutval=0.05) {
  ## T-test for probe occupancy differences between cols (i.e. positions)
  ## of two probe norm occupancy matrices
  ## Note there should be no NAs in input - usually => 0
  ##  E.G.
  ##   t.tp.dp <- ttest.probe(val.tp, val.dp)
  ##
  if(ncol(Aprobe.norm) != ncol(Bprobe.norm)) {
    stop("Two input matrices must have same number of cols/positions")
  }
  p <- rep(NA, ncol(val.tp))
  for(i in 1:ncol(val.tp)) {
    p[i] <- try(t.test(x=Aprobe.norm[,i], y=Bprobe.norm[,i],
                       alternative = "two.sided")$p.value, silent=TRUE)
  }
  p <- as.numeric(p)
  return( which(p<cutval) )
}



dis.ftmap <- function(pdbfile, prot.chains=1) {
  ##
  ## Return a minimal distance matrix of the min
  ## dist of a probe chain atom (clustered/grouped
  ## probes from ftmap) to any atom of a particular
  ## protein residue (i.e. min dist between probe
  ## and residue.
  ##  e.g.   m <- dis.ftmap( 1NF3_B.fftmap.pdb )
  ##
  
  ## Read FTMAP pdb format result file
  pdb <- read.pdb(pdbfile)

  ## What chains do we have
  chains <- unique(pdb$atom[,"chain"])

  ## Select Calpha chain 1 indices.              ### <- change per chains
  ca.inds <- atom.select(pdb, chain=chains[prot.chains], elety="CA") 

  ## Coords for protein
  ca.xyz <- pdb$xyz[ca.inds$xyz]

  ## blank matrix to store distance to 
  dmin.allres <- matrix(NA, nrow=length(ca.inds$atom),
                        ncol=(length(chains)-1)) ### <- change per chains

  nc <- max(prot.chains)+1
  for(i in nc:length(chains)) {                   ### <- change per chains

    ## Probe chains
    probe.inds <- atom.select(pdb, chain=chains[i])
    probe.xyz <- pdb$xyz[probe.inds$xyz]
    
    d <- dist.xyz(ca.xyz, probe.xyz, all.pairs=TRUE)
    d.min <- apply(d, 1, min)
    dmin.allres[,(i-1)] <- d.min
  }

  ##rownames(dmin.allres) <- pdb$atom[ca.inds$atom,"resno"]
  colnames(dmin.allres) <- chains[-prot.chains]
  return(dmin.allres)
}




bootstrap.rmsf <- function(xyz, N=5, n.repeat=10) {
  ## Bootstrap sampling of 'N' frames/rows from
  ## 'xyz' for RMSF determination.
  ##
  ## Here we take a sample of size 'N' from 'xyz' for
  ## RMSF calculation and repeat 'n.repeat' times.
  ## e.g. rf <- bootstrap.rmsf(trj)
  
  ## Create a blank RMSF storage matrix
  n.frames <- nrow(xyz)
  n.atoms  <- ncol(xyz)/3
  rf.store <- matrix(NA, ncol=n.atoms, nrow=n.repeat)

  for(i in 1:n.repeat) {
    ## bootstrap sampling of frames/rows
    row.inds <- sample(1:n.frames, size=N, replace=TRUE)
    rf.store[i,] <- rmsf(xyz[row.inds,])
  }

  return( list(ave=apply(rf.store, 2, mean),
               stdev=apply(rf.store, 2, sd),
               rf=rf.store) )
}



mustang <- function(pdbfiles, alnfile="mustangaln") {
  ##
  ## Structural alignment with mustang
  ## pdbfiles should be a vector of PDB files names
  ## alnfile is the output fasta format alignment
  ##
  if(length(pdbfiles) > 6) 
    warning("Multiple structural alignment will be VERY slow for this many files")
    
  cmd <- paste("~/bin/mustang -o", alnfile,
               "-F fasta -s OFF -r ON -i",
               paste(pdbfiles, collapse=" "))
  cat("\nRunning:\n",cmd,"\n\n")
  system(cmd)
  aln <- read.fasta( paste(alnfile,".afasta",sep="") )
  aln$ali[aln$ali=="Z"] <- "-"
  gaps <- gap.inspect(aln$ali)
  aln$ali <- aln$ali[,gaps$col!=nrow(aln$ali)]
  return(aln)
}


### --- See add.dccm.grid() below
add.grid.old <- function(sse, col="gray", lty="dashed") {
  ## Add a grid to a plot.dccm() plot
  ##   sse <- dssp(pdb)
  ##   plot.dccm(cij, sse=see)
  ##   add.grid(sse)

  if(is.list(sse)) {
    gridpoints <- c(sse$sheet$start[sse$sheet$length>4],
                    sse$helix$start[sse$helix$length>4] )
  } else {
    gridpoints <- sse
  }

  grid.rect(x=unit(gridpoints, "native"),
            y=unit(gridpoints, "native"),
            gp = gpar(fill=NA, col=col, lty=lty),
            just=c("right","top"),
            width=1, height=1,
            vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))

  grid.rect(x=unit(gridpoints, "native"),
            y=unit(gridpoints, "native"),
            gp = gpar(fill=NA, col=col, lty=lty),
            just=c("left","bottom"),
            width=1, height=1,
            vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))

  grid.lines(x=unit(c(0,max(gridpoints)), "native"),y=c(0,0),
             gp = gpar(fill=NA, col="black"),
             vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
  grid.lines(x=c(0,0),y=unit(c(max(gridpoints),0), "native"),
             gp = gpar(fill=NA, col="black"),
             vp=vpPath("plot_01.toplevel.vp","plot_01.panel.1.1.vp"))
}


add.dccm.grid <- function(x, fill.col=NA, helix.col="purple", sheet.col="yellow",
                          line.col="black", lty=1, lwd=1, segment.min=1,
                          alpha=0.3, side=c(1,2)) {
  
  ## Add a grid or colored boxes to a plot.dccm() plot
  ##   sse <- dssp(pdb)
  ##   plot.dccm(cij, sse=see)
  ##   add.dccm.grid(sse)
  ##   add.dccm.grid(list("start"=c(50, 100, 150), "length"=c(10,25,50)))
  ##   add.dccm.grid(list("start"=c(50, 100, 150), "length"=c(10,25,50)),
  ##                 fill.col=c("blue","red","green"))
  ##   plot.dccm2(cij, margin.segments=net$membership, segment.min=15)
  ##   add.dccm.grid( net$membership, segment.min=15)
  ##
  ##   add.dccm.grid( net$membership, segment.min=25)
  ##   add.dccm.grid( net$membership, segment.min=25, fill.col="gray")

  ##-- Function to draw box on plot
  draw.box <- function(start, length, xymin=0, xymax=1,
                       fill.col="gray", alpha=0.3, line.col="black", lty=1, lwd=1, 
                       side=c(1,2) ) {
    
    ##-- Draw Annotation Blocks On DCCM Plots
    ##    draw.box(150,20) 
    ##    draw.box(50,100, side=1, fill.col=NA, lty=2)

    ## Grid graphics paramaters
    gp <- gpar(fill=fill.col, col=line.col,
               lty=lty, lwd=lwd, alpha=alpha)

    vp <- vpPath("plot_01.toplevel.vp",
                 "plot_01.panel.1.1.vp")

    ##- Side 1: From Bottom Margin
    if( (side==1) || (side=="both") || all(side==c(1,2)) ) {
      grid.rect(x=unit(start-0.5, "native"),
                y=xymin,
                width=unit(length-0.5, "native"),
                height=xymax,
                gp=gp, just=c("left","bottom"), vp=vp) 
    }
    
    ##- Side 2: From Left Margin
    if( (side==2) || (side=="both") || all(side==c(1,2)) ) {
      grid.rect(x=xymin, 
                y=unit(start-0.5, "native"),
                width=xymax,
                height=unit(length-0.5, "native"),
                gp=gp, just=c("left","bottom"), vp=vp)      
    }
  }


  ##NOTE: For all 'x' objects that are not vectors we will exclude
  ##      segments that are under 'segment.min' in length
  segment.min.exclusion = TRUE 
  


  ##-- Parse 'x' Input --##
  
  ##- For vector input objects - e.g. $membership from cutree()
  if( (is.vector(x)) && (!is.list(x)) ) {
    grps <- table(x)

    ## Exclude small grps less than 'segment.min'
    ##  But do not do more filtering below 
    grps = names( grps[grps > segment.min] )
    segment.min.exclusion=FALSE ##<--- good idea but plots are too crowded!!

    store.grps <- NULL; 
    for(i in 1:length(grps)) {
      store.grps <- rbind(store.grps,
          cbind( bounds(which(x == grps[i])),
                "grp"=as.numeric(grps[i])) )
    }
    ## convert to matrix for use below
    x=store.grps

    ## Dont do any more filtering
    
  }

  ##- For SSE objects
  if(class(x) == "dssp") {
    start <- c(x$helix$start, x$sheet$start)
    length <- c(x$helix$length, x$sheet$length)
    
    ## If no 'fill.col' is provided use helix and sheet specific fill.col
    if(is.na(fill.col)) {
      fill.col <- c(rep(helix.col, length(x$helix$start)),
                    rep(sheet.col, length(x$sheet$start)) )
    }
  }  

  ##- For other list objects
  if(class(x) == "list") {
    start <- x$start
    length <- x$length
  }

  ##- For matrix objects - e.g. from bounds()
  if(class(x) =="matrix") {
    start  <- x[,"start"]
    length <- x[,"length"]
  }


  ##-- Filter out short segments based on input 'segment.min'
  if(segment.min.exclusion) {
    inds <- !(length < segment.min)
    start <- start[inds]
    length <- length[inds]
    if(length(fill.col > 1)) { fill.col <- fill.col[inds] }
  }
  
  ## Draw
  draw.box(start, length, fill.col=fill.col)
}

  





col.wheel <- function(str, cex=0.75) {
  ## Produce a simple chart for picking color names
  ##  e.g. col.wheel("slate")
  ##       col.wheel("dark")
  ##       col.wheel("red") # etc...
  cols <- colors()[grep(str, colors())]
  pie(rep(1, length(cols)), labels=cols, col=cols, cex=cex)
  cols
}
