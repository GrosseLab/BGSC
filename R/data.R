
#' Expression data set of Glioblastoma cell line SF767 
#'
#' A dataset containing the EListRaw objekt from limma::read.ilmn(files="Einzelanalyse_nonorm_nobkgd_SF767-1-799-6.txt", ctrlfiles="ControlProbeProfile_SF767-1-799-6.txt", other.columns=c("Detection","BEAD_STDERR","Avg_NBEADS") )
#' @format lists : "source" ,"E"   ,    "genes" ,  "other" ,  "targets"
#' @format "E" is a data frame with 48209 rows and 6 variables:
#' \describe{
#'   \item{SF 767 1}{expression}
#'   \item{SF 767 2}{expression}
#'   ...
#'   \item{SF 767 6}{expression}
#' }
#' 
"ExpData"


#' Expression data set of all qPCRs
#'
#' A dataset containing the
#' 
#' @format A data frame with 30 rows and 14 variables:
#' \describe{
#'   \item{Probe}{Name of celline}
#'   \item{NR}{Sample number}
#'   \item{CKAP2L_ctGen.ctRef}{deltaCt of CKAP2L}
#'   \item{CKAP2L_relative.Werte}{relativ expression of CKAP2L}
#'   ...
#'   \item{GALNSctGen.ctRef}{deltaCt of GALNS}
#'   \item{GALNSrelative.Werte}{relativ expression of GALNS}
#' }
#' 
"RawQPCRsf"

#' Expression data set of all qPCRs
#' GEOD.adf <-  fread('Illumina HumanHT-12_V4.0_A-GEOD-13475.adf.txt')
#' setnames(GEOD.adf,'Reporter Database Entry [genbank]','genbank')
#' setnames(GEOD.adf,'Reporter Name','illumina_humanht_12_v4')
#' REFSEQ_MRNA <- sapply(GEOD.adf[['genbank']], function(x) stringr::str_split(x,'\\.') )
#' GEOD.adf$REFSEQ_MRNA <- sapply(REFSEQ_MRNA, "[[", 1 ) 
#' setkey(GEOD.adf,'illumina_humanht_12_v4')
#' Illumina_to_REFSEQ_MRNA <- GEOD.adf
#' devtools::use_data(Illumina_to_REFSEQ_MRNA)
#' A dataset containing the
#' 
#' @format A data frame with 47323 rows and 5 variables:
#'
"Illumina_to_REFSEQ_MRNA" 
