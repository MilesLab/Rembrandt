% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_well_barcodes.R
\name{split_well_barcodes}
\alias{split_well_barcodes}
\title{Split by well barcodes}
\usage{
split_well_barcodes(
  pair1_filepath,
  pair2_filepath,
  reverse_barcode,
  nreads = 100,
  output_prefix = "outputfile"
)
}
\arguments{
\item{pair1_filepath}{file path of the first paired end file}

\item{pair2_filepath}{file path of the second paired end file}

\item{reverse_barcode}{A character vector with the reverse barcodes}

\item{nreads}{number of reads to iterate by}

\item{output_prefix}{The output prefix for the output files}
}
\description{
Split fastq files by the well barcodes provided
}
