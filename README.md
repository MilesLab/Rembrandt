# Rembrandt

The following repository will contain R codes for analyses of paired end COVID-19 sequencing data from Rembrandt. At the moment, this is still under construction, however a protocol for the analysis is provided.  

# Required Tools

R version 4.0 will be needed to perform the analyses. Additionally, the ShortRead and the RSubread packages will need to be obtained from Bioconductor.

The ShortRead package can be installed using the following lines.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ShortRead")

```


