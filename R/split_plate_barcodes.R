#' Split by plate barcodes
#'
#' Split fastq files by the plate barcodes provided
#'
#' @param pair1_filepath file path of the first paired end file
#' @param pair2_filepath file path of the second paired end file
#' @param nreads number of reads to iterate by
#' @param output_prefix The output prefix for the output files
#' @param forward_barcode A character vector with the foward barcodes
#'
#' @export
#'

split_plate_barcodes <- function(pair1_filepath, pair2_filepath,
                                 forward_barcode, nreads = 100,
                                 output_prefix = "outputfile"){

forward_barcoded_files1 = c()
forward_barcoded_files2 = c()
for(i in 1:length(forward_barcode)){
    forward_barcoded_files1[i] = paste(output_prefix, "_split_", forward_barcode[i], "_1.fq", sep = "")
    forward_barcoded_files2[i] = paste(output_prefix, "_split_", forward_barcode[i], "_2.fq", sep = "")
}


f1 <- ShortRead::FastqStreamer(pair1_filepath, nreads)
f2 <- ShortRead::FastqStreamer(pair2_filepath, nreads)
fq1 <- ShortRead::yield(f1)
fq2 <- ShortRead::yield(f2)

while(length(fq1)){


  for(i in 1:length(forward_barcode)){

    selected_reads = grepl(pattern = forward_barcode[i], ShortRead::sread(fq1))
    fqsub1 <- fq1[selected_reads]
    fqsub2 <- fq2[selected_reads]

    file.name1 = forward_barcoded_files1[i]
    file.name2 = forward_barcoded_files2[i]
    ShortRead::writeFastq(fqsub1,file.name1,mode="a",compress=FALSE)
    ShortRead::writeFastq(fqsub2,file.name2,mode="a",compress=FALSE)
  }
  fq1 <- ShortRead::yield(f1)
  fq2 <- ShortRead::yield(f2)


}

split_plate_barcodes = data.frame(forward_barcode,
                                  forward_barcoded_files1,
                                  forward_barcoded_files2,
                                  stringsAsFactors = F)

return(split_plate_barcodes)
}



