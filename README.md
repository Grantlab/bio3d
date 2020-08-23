# Bio3D Package for Biological Structure Analysis #

Utilities to analyze, process, organize and explore biomolecular structure, sequence and dynamics data.


## Features ##

Features include the ability to read and write structure, sequence and dynamic trajectory data, perform database searches, atom summaries, atom selection, re-orientation, superposition, rigid core identification, clustering, torsion analysis, distance matrix analysis, structure and sequence conservation analysis, normal mode analysis (NMA), correlation network analysis (CNA) and principal component analysis (PCA).  

In addition, various utility functions are provided to enable the statistical and graphical power of the R environment to work with biological sequence and structural data.  Please refer to the main [Bio3D website](http://thegrantlab.org/bio3d/) for more background information.


## Installing Bio3D ##

For the majority of users we recommend the use of the last stable release available from [CRAN](http://cran.r-project.org/web/packages/bio3d/) and the main [Bio3D website](http://thegrantlab.org/bio3d/). To install from within R issue the command:

```
#!r
install.packages("bio3d", dependencies=TRUE)
```

The development version is available from our bitbucket repository and typically contains new functions and bug fixes that have not yet been incorporated into the latest stable release. The simplest method for development version installation is to use the R function install_bitbucket() from the devtools package:


```
#!r
install.packages("devtools")
library(devtools)
install_bitbucket("Grantlab/bio3d", subdir = "ver_devel/bio3d/")
```

Alternative installation methods and additional instructions are posted to the [wiki section](https://bitbucket.org/Grantlab/bio3d/wiki/Home) of our bitbucket repository. 


## Other packages in the Bio3D familiy ##

With growing functionality in Bio3D we have decided to implement larger Bio3D modules into separate packages:

* [Bio3D-web](https://bitbucket.org/Grantlab/bio3d-web/) for online interactive analysis of protein structure ensembles
* [Bio3D-nma](https://bitbucket.org/Grantlab/bio3d-nma/) for Normal Mode Analysis (NMA) of protein structures and ensembles
* [Bio3D-cna](https://bitbucket.org/Grantlab/bio3d-cna/) for Correlation Network Analysis (CNA) for protein structure ensembles
* [Bio3D-eddm](https://bitbucket.org/Grantlab/bio3d-eddm/) for Distances matrix analysis for the identification of significant conformational changes underling functional processes
* [Bio3D-gesostas](https://bitbucket.org/Grantlab/bio3d-geostas/) for the identification of geometrically stable domains in biomolecules
* [Bio3D-view](https://bitbucket.org/Grantlab/bio3d-view/) for interactive visualization of structures in R

Current ongoing projects entails:

* [Bio3D-mmtf](https://bitbucket.org/Grantlab/bio3d-mmtf/) [In development]
* [Bio3D-amber](https://bitbucket.org/Grantlab/bio3d-amber/) [In development]
* [Bio3D-muscle](https://bitbucket.org/Grantlab/bio3d-muscle/) [In development]
* [Bio3D-cheminf](https://bitbucket.org/larsss/cheminf/src/master/) 
* [Bio3D-dcna](https://bitbucket.org/xinqyao/dcna)

## Contributing to Bio3D ##

We are always interested in adding additional functionality to Bio3D. If you have ideas, suggestions or code that you would like to distribute as part of this package, please contact us (see below). You are also encouraged to contribute your code or issues directly to [this repository](https://bitbucket.org/Grantlab/bio3d/) for incorporation into the development version of the package. For details on how to do this please see the [developer wiki](https://bitbucket.org/Grantlab/bio3d/wiki/Home).  

  
## Contact ##

You are welcome to:

* Submit suggestions and bug-reports at: https://bitbucket.org/Grantlab/bio3d/issues
* Send a pull request on: https://bitbucket.org/Grantlab/bio3d/pull-requests
* Compose a friendly e-mail to: bjgrant@ucsd.edu
