# Rembrandt


The following repository will contain R codes for analyses of paired end COVID-19 sequencing data from Rembrandt. At the moment, this is still under construction, however a protocol for the analysis is provided. It is our goal to create an R package to assist with analysis. 

# Required Tools

R version 4.0 will be needed to perform the analyses. Additionally, the ShortRead and the RSubread packages will need to be obtained from Bioconductor.

The ShortRead package can be installed using the following code.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ShortRead")
```

The RSubread package can be installed using the following code.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
```

# Protocol for Analyses

## De-multiplexing Reads

Prior to aligning the ‘fastq’ files, we need to demultiplex RNASeq reads using the ShortRead Bioconductor R package. This will separate the initial multiplexed ‘fastq’ file into smaller files for each patient sample. Of note, given that we are using paired-end data, both files must be ordered in the same way and care must be ensured to ensure that this order is preserved. The steps for de-multiplexing reads are as follows:

1. Make CSV files containing lists of Illumina adaptor barcodes, plate barcodes, and well or patient barcodes.  These can be labeled as ‘illumina_barcodes.csv’, ‘plate_barcodes’ and ‘patient_barcodes’.

2. Make a patient mapping CSV file containing 4 columns linking patient ID to each of the barcodes in the format as below.
```
Patient_ID		Illumina_barcode		Plate_barcode		Patient_barcode
```
3. Demultiplexing based on Illumina Barcode Sequences
•	Use FastqStreamer() function to build a connection that will allow us to import a set number of reads (nreads) from each paired end file as below:
```
f1 <- FastqStreamer({pair1_filepath}, nreads)
f2 <- FastqStreamer({pair2_filepath}, nreads)
```

•	The set number of reads specified by nreads will be extracted using the yield() function as below in steps until there are no more reads available
```
fq1 <- yield(f1))
fq2 <- yield(f2))
```
•	For each set of extracted reads, a for-loop will run over the set of Illumina barcodes
•	In the for-loop, the grepl() string searching function will search for the selected barcode. If a barcode is present, the reads will be written to a ‘fastq’ file for that particular Illumina barcode.  Also, the corresponding reads from the other mate will be written to the other ‘fastq’ file. 
```
selected_reads = grepl(pattern, sread(fq1))
fqsub1 <- fq1[selected_reads]
fqsub2 <- fq2[selected_reads]
writeFastq(fqsub1,{fq_illumina_barcodeX_pair1_filepath},mode="a",compress=FALSE)
writeFastq(fqsub2,{fq_illumina_barcodeX_pair2_filepath},mode="a",compress=FALSE)
```
•	Final ‘fastq’ files from Illumina demultiplexing will be used for demultiplexing by plate-specific barcodes

4. Demultiplexing based on plate barcodes
•	Similar to demultiplexing based on Illumina barcodes
•	Procedures repeated for each resulting file from Illumina demultiplexing
•	Use FastqStreamer() function to build a connection that will allow us to import a set number of reads (nreads) from each paired end file as below:
```
f1 <- FastqStreamer({fq_illumina_barcodeX_pair1_filepath}, nreads)
f2 <- FastqStreamer({fq_illumina_barcodeX_pair2_filepath}, nreads)
```
•	The set number of reads specified by nreads will be extracted using the yield() function as below in steps until there are no more reads available
```
fq1 <- yield(f1))
fq2 <- yield(f2))
```
•	For each set of extracted reads, a for-loop will run over sets of plate barcodes
•	In the for-loop, the grepl() string searching function will search for the selected barcode. If a barcode is present, the reads will be written to a ‘fastq’ file for selected plate barcode.  Also, the corresponding reads from the other mate will be written to the other ‘fastq’ file. 
```
selected_reads = grepl(pattern, sread(fq1))
fqsub1 <- fq1[selected_reads]
fqsub2 <- fq2[selected_reads]
writeFastq(fqsub1,{fq_plate_barcodeY_illumina_barcodeX_pair1_filepath},mode="a",
compress=FALSE)
writeFastq(fqsub2,{fq_plate_barcodeY_illumina_barcodeX_pair2_filepath},mode="a",
compress=FALSE)
```
•Final ‘fastq’ files from plate demultiplexing will be used for patient-specific barcodes

5. Demultiplexing based on patient
•	Similar to demultiplexing based on Illumina and plate-specific barcode
•	Procedures repeated for each resulting file from plate demultiplexing
•	Use FastqStreamer() function to build a connection that will allow us to import a set number of reads (nreads) from each paired end file as below:
```
f1 <- FastqStreamer{(fq_plate_barcodeY_illumina_barcodeX_pair1_filepath}, nreads)
f2 <- FastqStreamer({fq_plate_barcodeY_illumina_barcodeX_pair2_filepath}, nreads)
```
•	The set number of reads specified by nreads will be extracted using the yield() function as below in steps until there are no more reads available
```
fq1 <- yield(f1))
fq2 <- yield(f2))
```
•	For each set of extracted reads, a for-loop will run over sets of patient barcodes
•	In the for-loop, the grepl() string searching function will search for the select barcode on the second mate given the patient barcode is on the reverse pair. If a barcode is present, the reads will be written to a ‘fastq’ file for selected patient barcode.  Also, the corresponding reads from the other mate will be written to the other ‘fastq’ file. 
```
selected_reads = grepl(pattern, sread(fq2))
fqsub1 <- fq1[selected_reads]
fqsub2 <- fq2[selected_reads]
writeFastq(fqsub1,{fq_patient_barcodeZ_plate_barcodeY_illumina_barcodeX_pair1_filepath},mode="a",compress=FALSE)
writeFastq(fqsub2,{fq_patient_barcodeZ_plate_barcode_Y_illumina_barcodeX_pair2_filepath},mode="a",compress=FALSE)
```
•	Final ‘fastq’ files from patient demultiplexing will be used for alignments

6. Rename resulting ‘fastq’ files. This step can be done at the end also. However it may be good practice to change the names of the fastq files to patient IDs using the patient mapping CSV file early on. 

## Mapping Reads and Determining Gene Counts

After the files have been demultiplexed, we can then proceed with the process of mapping reads to genes using RSubread Bioconductor R package using the following steps. 

1. Constructing ‘fasta’ reference file

Prior to mapping, we will need to assemble the  nCoV, Human RPP30, and SARS CoV sequences into one ‘fasta’ file. This can be done using cat in command line:
```
cat ncov.fa hs_rpp30.fa sars_cov.fa > cov_ref.fa
```
2. Build the RSubread index for alignments

•	Use created ‘fasta’ reference file to build the mapping index as below
```
buildindex(basename="./reference_index",reference=”cov_ref.fa”)
```
3. Create annotation file for gene counts
•	A Simplified Annotation Format (SAF) file should be created for each of the genes in the format below.
```
GeneID  Chr Start End Strand
```

4.  Map ‘fastq’ reads and determine gene counts for each patient to nCoV, Human RPP30, and SARS CoV
•	Using built index and selected ‘fastq’ files for a patient, align reads to genes as below. 
```
align.stat <- align(index = "./reference_index", readfile1 = reads1, readfile2 = reads2, output_file = "./{Rsubread_alignment}.bam")
```
•	Using alignment bam files as well as the SAF file, obtain gene counts using the featureCounts() function as below.  
```
featureCounts(files=“./{Rsubread_alignment}.bam”, annot.ext=saf_file)
```
•	The gene counts can then be assembled all patients into matrix that can be used further analyses.

## Running the Routines

Of note, users can choose different ways of running these routines. The straightforward way of running the routines is demultiplex successively by each set of barcodes and then perform the alignments as a loop. However this can result in a large number of files waiting to be alignment. One solution that we will incorporate is the parLappy() function from the snow R package to parallelize the alignment steps. Additionally, the user may start by splitting files based on Illumina barcodes and then demultiplex and align files with one type of Illumina barcode before proceeding to the next Illumina barcodes. The different routines can be run based on the system’s needs and specifications.  




