This application demonstrates normal modes analysis (NMA) with Bio3D. It uses function `read.pdb` to fetch and parse the requested PDB structure; `trim.pdb` to trim the PDB to specified chain ID; `nma` to calculate the normal modes; and `dccm` to calculate the cross correlation matrix.

The analysis can be performed with the following code:

```
pdb = read.pdb("1hel")
pdb = trim.pdb(pdb, chain="A")

modes = nma(pdb)
plot(modes)

cij = dccm(modes)
plot(cij)


```