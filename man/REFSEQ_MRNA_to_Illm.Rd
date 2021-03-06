% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BGSC.R
\name{REFSEQ_MRNA_to_Illm}
\alias{REFSEQ_MRNA_to_Illm}
\title{Convert REFSEQ_MRNA string from DAVID to Illumina ids}
\usage{
REFSEQ_MRNA_to_Illm(string, asString = FALSE)
}
\arguments{
\item{string}{with REFSEQ_MRNA ids example:  "NM_032169, NM_000903, NM_015913"}

\item{asString}{if TRUE return a comma seperagted string with Illumina Ids}
}
\value{
a \code{vector} with Illumina Ids
}
\description{
Convert REFSEQ_MRNA string from DAVID to Illumina ids
}
\author{
Claus Weinholdt
}
