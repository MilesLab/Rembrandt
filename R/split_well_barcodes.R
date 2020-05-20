#' Split by well barcodes
#'
#' Split fastq files by the well barcodes provided
#'
#' @param pair1_filepath file path of the first paired end file
#' @param pair2_filepath file path of the second paired end file
#' @param nreads number of reads to iterate by
#' @param output_prefix The output prefix for the output files
#' @param reverse_barcode A character vector with the reverse barcodes
#'
#' @export
#'

split_well_barcodes <- function(pair1_filepath, pair2_filepath,
                                 reverse_barcode, nreads = 100,
                                 output_prefix = "outputfile"){

reverse_barcoded_files1 = c()
reverse_barcoded_files2 = c()
for(i in 1:length(reverse_barcode)){
    reverse_barcoded_files1[i] = paste(output_prefix, "_split_", reverse_barcode[i], "_1.fq", sep = "")
    reverse_barcoded_files2[i] = paste(output_prefix, "_split_", reverse_barcode[i], "_2.fq", sep = "")
}


f1 <- ShortRead::FastqStreamer(pair1_filepath, nreads)
f2 <- ShortRead::FastqStreamer(pair2_filepath, nreads)
fq1 <- ShortRead::yield(f1)
fq2 <- ShortRead::yield(f2)
on.exit(close(fq1))
on.exit(close(fq2))
while(length(fq1)){


  for(i in 1:length(reverse_barcode)){

    selected_reads = grepl(pattern = reverse_barcode[i], ShortRead::sread(fq2))
    fqsub1 <- fq1[selected_reads]
    fqsub2 <- fq2[selected_reads]

    file.name1 = reverse_barcoded_files1[i]
    file.name2 = reverse_barcoded_files2[i]
    ShortRead::writeFastq(fqsub1,file.name1,mode="a",compress=FALSE)
    ShortRead::writeFastq(fqsub2,file.name2,mode="a",compress=FALSE)
  }
  fq1 <- ShortRead::yield(f1)
  fq2 <- ShortRead::yield(f2)


}

split_well_barcodes = data.frame(reverse_barcode,
                                  reverse_barcoded_files1,
                                  reverse_barcoded_files2,
                                  stringsAsFactors = F)

return(split_well_barcodes)
}



