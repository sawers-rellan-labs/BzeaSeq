% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/introgression_plots.R
\name{convert_to_matrix}
\alias{convert_to_matrix}
\title{Convert Run-Length Encoded Data to Matrix Format}
\usage{
convert_to_matrix(
  processed_data,
  bin_size = 1e+06,
  chromosomes = paste0("chr", 1:10),
  fill_missing = TRUE
)
}
\arguments{
\item{processed_data}{Data frame with processed introgression segments}

\item{bin_size}{Integer. Genomic bin size in base pairs (default: 1000000)}

\item{chromosomes}{Character vector. Chromosomes to include (default: chr1-chr10)}

\item{fill_missing}{Logical. Whether to fill missing bins with "REF" (default: TRUE)}
}
\value{
List containing:
\item{matrix}{Matrix with samples as columns, genomic bins as rows}
\item{bin_metadata}{Data frame with bin position information}
\item{chromosome_lengths}{Data frame with chromosome length information}
}
\description{
Converts segment-based introgression data to a matrix format suitable
for ComplexHeatmap visualization and other matrix-based analyses.
}
\details{
Converts variable-length segments to fixed-size genomic bins for
matrix-based visualization. Each bin is assigned the genotype of
overlapping segments. Missing bins are filled with "REF" by default.
}
\examples{
\dontrun{
# Convert to matrix format
matrix_data <- convert_to_matrix(processed_data, bin_size = 1000000)

# Access the genotype matrix
genotype_matrix <- matrix_data$matrix

# Get bin metadata for annotations
bin_info <- matrix_data$bin_metadata
}

}
