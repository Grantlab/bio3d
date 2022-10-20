# v2.4 (Version 2.4, released November 2019)

## Overview of new features and enhancements

* New facilities for sequence alignment using the 'msa' package from Bioconductor or alternatively using the online server of the European Bioinformatics Institute (EMBL-EBI) (and so a local MUSCLE program is no longer a mandate).

* A more robust sequence fetching function supporting both the EMBL-EBI and NCBI (National Center for Biotechnology Information) servers.

* A new function for "core" detection using contact maps, new and improved functions for processing protein structures and doing structural network analysis.

* More advanced system environment checking for automatically locating essential external programs.  

## Other updates

* We have also updated online vignettes and other documentations. For a fine-grained list of changes, or to report a bug, please consult:

    * [The issues log](https://bitbucket.org/Grantlab/bio3d/issues)
    * [The commit log](https://bitbucket.org/Grantlab/bio3d/commits/all)

* For full install instructions see: 
    * <http://thegrantlab.org/bio3d/tutorials/installing-bio3d>

## Major new/enhanced functions

* seqaln:  Supports using 'msa' and the EMBL-EBI server to do sequence alignment 
* get.seq:  More robust to fetch a large number of sequences and supports EMBL-EBI and NCBI servers
* atom.select.pdbs:  New function for atom selection of a 'pdbs' object
* core.cmap:  New function for "core" residues detection using contact maps
* chain.pdb:  Supports inspection on peptide bond (C-N) length
* community.aln:  Comparison of networks with different numbers of nodes
* dccm.pca, dccm.xyz:  Tidier interface by merging "lmi()" and providing a new argument "method"
* core.find, dccm.nma, dm.xyz, nma.pdbs, pdbaln, pdbsplit, read.fasta.pdb:  Improved progress bar
* dssp.pdb, pdbaln, pymol.dccm, pymol.modes, pymol.pdbs, seqaln:  Exefile argument with OS dependent defaults
* rmsd, seqidentity:  More sensible output by fetching row and column names from input

### Happy Bio3Ding!


# v2.3 (Version 2.3, released September 2016)

## Overview of new features and enhancements

* New facilities for ensemble normal mode analysis (NMA) with all-atom elastic network model (ENM) and Gaussian network model (GNM).

* Enhanced NMA calculations with the rotation-translation block (RTB) method and the new "4-bead" coarse-grained ENM.

* More efficient reading of large PDB files using Rcpp.

* PDB annotation from the PFAM database.

* More supported I/O file formats.  

## Other updates

* We have also updated online vignettes and other documentations. For a fine-grained list of changes, or to report a bug, please consult:

    * [The issues log](https://bitbucket.org/Grantlab/bio3d/issues)
    * [The commit log](https://bitbucket.org/Grantlab/bio3d/commits/all)

* For full install instructions see: 
    * <http://thegrantlab.org/bio3d/tutorials/installing-bio3d>

## Major new/enhanced functions

* aanma:  All-atom ENM normal mode analysis (with RTB and 4-bead ENM supported)
* aanma.pdbs:  Ensemble NMA with all-atom ENM
* gnm:  Gaussian network model (GNM) calculations
* gnm.pdbs:  Ensemble NMA with GNM
* dccm.gnm:  Dynamical cross-correlation for GNM
* pdbs2sse: Retrieve SSE from pdbs object with appropriate residue numbers for plotting
* mask.dccm:  Produce a new DCCM object with selected atoms masked
* pdb.pfam:  Function for PFAM annotation of PDB IDs
* pymol.pdbs:  Builds a pymol session from a 'pdbs' object
* read.cif:  Read a Protein Data Bank (mmCIF) coordinate file
* read.dssp:  For reading existing DSSP output files
* read.stride:  For reading existing STRIDE output files
* read.crd:  Read a CHARMM CARD (CRD) or AMBER coordinate file
* read.prmtop:  Read parameter and topology data from an AMBER PrmTop file
* read.pdb:  Use Rcpp to (more rapidly) read and parse PDB files
* read.pdb2: Renamed old read.pdb function
* plot.matrix.loadings:  For plotting loadings obtained from pca.array()
* community.aln:  To align communities from two or more related networks
* atom.select: Supports 'insert' identifier
* vmd.cna and vmd.cnapath: Renamed view.cna and view.cnapath
* pymol.dccm, pymol.modes, pymol.nma, and pymol.pca: Renamed view.xxx functions
* plot.fasta: Improved plotting function for multiple sequence alignment
* read.mol2, write.mol2, atom.select, trim, as.pdb:  Read, write and manipulate mol2 files with functions  

### Happy Bio3Ding!


# v2.2 (Version 2.2, released in Feb 2015) 

## Added new facilities for:

* sub-optimal path analysis of biomolecular correlation networks

* constructing biological units

* identification and tidying of malformed PDB files

* improved secondary structure annotation of 'pdbs' objects and various plots.

## Other updates

* We have also updated and enhanced atom selection functionality and developed a new vignette detailing PDB structure manipulation and analysis facilities. For a fine-grained list of changes, or to report a bug, please consult:

    * [The issues log](https://bitbucket.org/Grantlab/bio3d/issues)
    * [The commit log](https://bitbucket.org/Grantlab/bio3d/commits/all)

## Major new functions

* cnapath: Suboptimal Path Analysis for Correlation Networks
* biounit: Biological Unit Construction
* clean.pdb: Inspect And Clean Up A PDB Object
* cat.pdb: Concatenate Multiple PDB Objects
* pdb2sse: Obtain An SSE Sequence Vector From A PDB Object
* bounds.sse: Obtain A SSE Object From An SSE Sequence Vector
* aa.table: Updated amino acid reference data that replaces older 'aa.mass'
* as.fasta: Convert alignment/sequence in matrix/vector format to a FASTA object.
* as.pdb: Convert coordinate data to PDB format
* as.select: Convert atomic indices to a atom.select object
* as.xyz: Convert vectors and matrices to 'xyz' class objects
* atom.select.pdb: Atom selection from PDB objects has been extensively updated
* basename.pdb: Utility for manipulation of PDB file names
* check.utility: Check and Report on Missing Bio3D Utility Programs
* cmapt: Update contact map methods for pdb and xyz objects
* cna: Update correlation network analysis methods for dccm and ensmb objects"
* cnapath: Suboptimal Path Analysis for Correlation Networks
* com: Updated center of mass methods for pdb and xyz objects
* combine.select: Combine atom.select objects, renamed from previous 'combine.sel'
* cov.enma: New method to Calculate Covariance Matrix from Ensemble Normal Modes"
* cov2dccm: Calculates the N-by-N cross-correlation matrix from a 3N-by-3N covariance matrix
* covsoverlap: New methods for nma and enma objects
* dm: Distance matrix gets new methods for pdb and xyz class objects
* dssp: Secondary Structure Analysis with DSSP gets new methods for pdb, xyz and pdbs class objects
* geostas: Geometrically stable domain finder gets new methods for nma, enma, pdb, pdbs and xyz objects.
* is.pdbs: Is an Object of Class pdbs
* mono.colors: New color palette
* pdb2sse: Obtain An SSE Sequence Vector From A PDB Object
* pdbfit: Coordinate superposition gets new methods for multi-model pdb objects and pdbs objects.
* read.crd: Can Now Read Coordinate Data from Amber or CHARMM
* read.prmtop: Read AMBER Parameter/Topology files
* var.pdbs: Pairwise Distance Variance in Cartesian Coordinates
* plot:  New or updated plot methods for 'cmap', 'geostas', and 'pca' class objects as well as a new plot.fluct() function that expands on plot.bio3d() for plotting atomic fluctuations from MD and NMA results.
* print:  New print methods for cnapath, enma, geostas, mol2, nma, pca, pdb, prmtop, rle2,  select and sse objects.


# v2.1 (Version 2.1, released in Sep 2014)

## Overview of major changes

* Added new facilities for correlation Network Analysis (cna) and Geometrically Stable Domain finding (geostas).

* We have also changed 'PDB object data' storage from a matrix to a data.frame
format. 

* Improved methods and functionality for ensemble NMA are now also included along with extensive improvements to package vignettes and function documentation. For a fine-grained list of changes, or to report a bug,
please consult:

    * [The issues log](https://bitbucket.org/Grantlab/bio3d/issues)
    * [The commit log](https://bitbucket.org/Grantlab/bio3d/commits/all)

## Major new functions

* cna: Protein Dynamic Correlation Network Construction and Community Analysis.
* plot.cna: Protein Structure Network Plots in 2D and 3D.
* print.cna: Summarize and Print Features of a cna Network Graph
* identify.cna: Identify Points in a CNA Protein Structure Network Plot
* layout.cna: Protein Structure Network Layout
* view.cna: View CNA Protein Structure Network Community Output in VMD
* prune.cna: Prune A cna Network Object
* community.tree: Reconstruction of the Girvan-Newman Community Tree for a CNA Class Object.
* network.amendment: Amendment of a CNA Network According To A Input Community Membership Vector.
* lmi: Linear Mutual Information Matrix
* dccm.pca: Dynamic Cross-Correlation from Principal Component Analysis
* filter.dccm: Filter for Cross-correlation Matrices (Cij)
* cmap.filter: Contact Map Consensus Filtering
* geostas (amsm.xyz): GeoStaS Domain Finder
* bhattacharyya Bhattacharyya Coefficient
* covsoverlap: Covariance Overlap
* sip: Square Inner Product
* cov.nma: Calculate Covariance Matrix from Normal Modes
* mktrj.enma: Ensemble NMA Atomic Displacement Trajectory
* pca.array: Principal Component Analysis of an array of matrices
* hmmer: HMMER Sequence Search
* plot.hmmer: Plot a Summary of HMMER Hit Statistics.
* uniprot: Fetch UniProt Entry Data.
* pfam: Download Pfam FASTA Sequence Alignment
* hclustplot: Dendrogram with Clustering Annotation
* write.pir: Write PIR Formated Sequences
* mustang: Structure-based Sequence Alignment with MUSTANG
* pdbs.filter: Filter or Trim a pdbs PDBs Object
* dssp.pdbs: Secondary Structure Analysis of Aligned PDB Structures with DSSP
* plot.fasta: Plot a Multiple Sequence Alignment
* print.fasta: Printing Sequence Alignments
* inspect.connectivity: Check the Connectivity of Protein Structures
* var.xyz: Pairwise Distance Variance in Cartesian Coordinates
* is.xyz(as.xyz, print.xyz): Is an Object of Class
* setup.ncore: Setup for Running Bio3D Functions using Multiple CPU Cores

# v2.0 (Version 2.0, released in Nov 2013)

* Contains over 30 new functions including enhanced Normal Mode Analysis facilities as well extensive improvements to existing code and documentation. For a fine-grained list of changes or to report a bug, please consult:

    * [The issues log](https://bitbucket.org/Grantlab/bio3d/issues)
    * [The commit log](https://bitbucket.org/Grantlab/bio3d/commits/all)

## Major new functions

* aa2mass: Amino Acid Residues to Mass Converter
* atom.index: Index of Atomic Masses
* atom2mass(atom2ele, formula2mass): Atom Names to Mass Converter
* binding.site: Binding Site Residues
* com(com.xyz): Center of Mass
* combine.sel: Combine Atom Selections From PDB Structure
* dccm.enma: Cross-Correlation for Ensemble NMA (eNMA)
* dccm.mean: Filter DCCM matrices
* dccm.nma: Dynamic Cross-Correlation from Normal Modes Analysis
* dccm.xyz: DCCM: Dynamical Cross-Correlation Matrix
* deformation.nma: Deformation Analysis
* dssp.trj: Secondary Structure Analysis of Trajectories with DSSP
* fluct.nma: NMA Fluctuations
* inner.prod: Mass-weighted Inner Product
* is.pdb: Is an Object of Class pdb
* is.select: Is an Object of Class atom.select
* load.enmff(ff.calpha, ff.calphax, ff.anm, ff.pfanm, ff.sdenm, ff.reach): ENM Force Field Loader
* mktrj.nma: NMA Atomic Displacement Trajectory
* nma(build.hessian, print.nma): Normal Mode Analysis
* nma.pdbs(print.enma): Ensemble Normal Mode Analysis
* normalize.vector: Mass-Weighted Normalized Vector
* pdb.annotate: Get Customizable Annotations From PDB
* pdb2aln: Align a PDB structure to an existing alignment
* pdb2aln.ind: Mapping between PDB atomic indices and alignment positions
* pdbfit: PDB File Coordinate Superposition
* pdbs2pdb: PDBs to PDB Converter
* plot.enma: Plot eNMA Results
* plot.nma: Plot NMA Results
* plot.rmsip: Plot RMSIP Results
* read.mol2: Read MOL2 File
* sdENM: Index for the sdENM ff
* sse.bridges: SSE Backbone Hydrogen Bonding
* struct.aln: Structure Alignment Of Two PDB Files
* view.dccm: Visualization of Dynamic Cross-Correlation
* view.modes: Vector Field Visualization of Modes
* vmd.colors: Color as in VMD Molecular Viewer


## Versioning

* Releases will be numbered with the following semantic versioning format:

    * \<major\>.\<minor\>-\<patch\>

    * E.g.: 2.0-1

* And constructed with the following guidelines:

    * Breaking backward compatibility bumps the major (and resets the minor and patch)
    * New additions without breaking backward compatibility bumps the minor (and resets the patch)
    * Bug fixes and misc changes bumps the patch

### For more information on SemVer, please visit http://semver.org/.


# For changes prior to v1.1-6 (Apr 2013) please see the bio3d wki:
* [Whats new wki page](http://bio3d.pbworks.com/w/page/7824486/WhatsNew)
