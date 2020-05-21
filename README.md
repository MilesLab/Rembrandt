# Rembrandt

This is the Github repository for the Rembrandt R package that contains functions for processing Rembrant data. 

# Installation

The package can be installed from GitHub as follows.

```
library(devtools)
install_github("MilesLab/Rembrandt")
```

# Required Tools

R version 3.5.3 or higher will be needed to perform the analyses. The package will also require *Shortread* and *Rsubread* Bioconductor packages. This should be done automatically when installing *Rembrandt* however we provide instructions on obtaining said packages. 

The *ShortRead* package can be installed using the following code.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ShortRead")
```

The *RSubread* package can be installed using the following code.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
```

# Current File Requirements

**Given that we are using paired-end data, it is a requirement that both pairs must be ordered in exactly the same way or these functions will not work correctly**

# Vignette and Data

The package includes a 'fasta' file and SAF annotation file for amplicon sequences as well as test data.

They can be obtained from the 'inst/ext' folder of the package or downloaded from the link below:

[https://mileslab.github.io/RembrandtData.zip](https://mileslab.github.io/RembrandtData.zip)

A vignette is provided for running this package. Due to space constraints, we host it on another repository.  The link for the 'html' vignette is as below:

[https://mileslab.github.io/UsingRembrandt.html](https://mileslab.github.io/UsingRembrandt.html)
 
The .Rmd file for the vignette is available at this link

[https://mileslab.github.io/vignette_rmd.zip](https://mileslab.github.io/vignette_rmd.zip)

The processed example output can be downloaded from here.

[https://mileslab.github.io/temp_vignette.zip](https://mileslab.github.io/temp_vignette.zip)

# Original Protocol

If you want to refer to our original protocol, refer here:

[https://mileslab.github.io/analysisprotocol.html](https://mileslab.github.io/analysisprotocol.html)

# Citation

If you use any of our codes, please cite as below:

**Palmieri, Dario, et al. "REMBRANDT: A high-throughput barcoded sequencing approach for COVID-19 screening." bioRxiv (2020).**

# Contact

If you have any questions, comments, and suggestions for our software, please contact <jalal.siddiqui@osumc.edu>

