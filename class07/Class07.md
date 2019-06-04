class07
================
Luan Tran
4/18/2019

``` r
library(bio3d)
```

``` r
#PlotPDB function reads any Protein PDB data (input) and outputs a plot for the specified protein residue against their b factor values
PlotPDB <- function(x) {
 PDB <- read.pdb(x) #Read a Protein Data Bank (PDB) coordinate file. Outputs protein information including the sequence
 chainA <- trim.pdb(PDB, Chain="A", elety ="CA") #Produces a new smaller PDB object from the original input PDB, In this case it returns the Chain A and C Alpha atoms
 bfactor <- chainA$atom$b #Specifies subset of PDB file, in this case the atoms in Chain A and extracting the b factor value
 plotb3(bfactor, sse=chainA, typ="l", ylab="Bfactor") #Outputs a plot of the bfactor with secondary structure elements (residue)
}
```

``` r
PlotPDB("4AKE")
```

    ##   Note: Accessing on-line PDB file

![](Luan_Tran_Hw6_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
PlotPDB("1AKE")
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

![](Luan_Tran_Hw6_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
PlotPDB("1E4Y")
```

    ##   Note: Accessing on-line PDB file

![](Luan_Tran_Hw6_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->
