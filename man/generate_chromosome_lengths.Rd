% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/introgression_plots.R
\name{generate_chromosome_lengths}
\alias{generate_chromosome_lengths}
\title{Generate Chromosome Length Data}
\usage{
generate_chromosome_lengths(assembly = "B73v5")
}
\arguments{
\item{assembly}{Character. Genome assembly version (currently only "B73v5" supported)}
}
\value{
Data frame with columns:
\item{chrom}{Chromosome name (chr1-chr10)}
\item{length}{Chromosome length in base pairs}
}
\description{
Returns exact chromosome lengths for the B73 reference genome assembly.
}
\details{
Chromosome lengths are from Zm-B73-REFERENCE-NAM-5.0.fa.fai.
Values represent the exact nuclear chromosome lengths used for
bin calculations and genome-wide visualizations.
}
\examples{
chrom_lengths <- generate_chromosome_lengths()
print(chrom_lengths)

}
