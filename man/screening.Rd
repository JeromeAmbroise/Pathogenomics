% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screening.R
\name{screening}
\alias{screening}
\title{Screen for a set of genes of interst (e.g. virulence genes) from a (draft of complete) genome}
\usage{
screening(genomePath, genesPath, lengthconf, identconf, outputdir)
}
\arguments{
\item{genomePath}{A string of character corresponding to the path of the FASTA file of the genome to be screened.}

\item{genesPath}{A vector of string of character corresponding to the paths of the FASTA files of the genes of interst (e.g. virulence genes).}

\item{lengthconf}{The minimum accepted length of the hit (in percentage of the gene of interest).}

\item{identconf}{The minimum accepted identity of the hit (in percentage of identity).}

\item{outputdir}{The directory where output files are stored.}
}
\value{
result is a vector containing 0 and 1. This vector has the length of the number of genes of interest.
The name of the elements of this vector corresponds to the name of the gene of interest.
}
\description{
screening returns a vector indicating if the genome matches each gene of interst.
}
\examples{
genomePath <- list.files(system.file("extdata/genome", package = "Pathogenomics"),full.names=T)
Ngenomes <- length(genomePath)
genesPath <- list.files(system.file("extdata/virulence/Escherichia-coli/Virulencefinder", package = "Pathogenomics"),full.names=T)
Ngenes <- length(genesPath)
screening(genomePath= genomePath[1],genesPath=genesPath,lengthconf = 60,identconf = 90,outputdir = '7-virulence-Irenge')
}
