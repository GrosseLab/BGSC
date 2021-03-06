
#' Expression data set of Glioblastoma cell line SF767 
#'
#' A dataset containing the EListRaw objekt from limma::read.ilmn(files="Einzelanalyse_nonorm_nobkgd_SF767-1.txt", ctrlfiles="ControlProbeProfile_SF767-1.txt", other.columns=c("Detection","BEAD_STDERR","Avg_NBEADS") )
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


#' Delta CT values of EGF treatment for SF767 and LNZ308
#'
#' A dataset containing the
#' 
#' 
#' @format A data frame with 12 rows and 9 variables:
#' \describe{
#'   \item{ID}{qPCR ID}
#'   \item{CellLine}{Cell line ID}
#'   \item{Treatment}{Treatment}
#'   \item{GAPDH}{deltaCt of GAPDH}
#'   \item{CKAP2L}{deltaCt of CKAP2L}
#'   \item{ROCK1}{deltaCt of ROCK1}
#'   \item{TPR}{deltaCt of TPR}
#'   \item{ALDH4A1}{deltaCt of ALDH4A1}
#'   \item{GALNS}{deltaCt of GALNS}
#' }
#' 
"qPCR_SF767_LNZ308"

#' Delta CT values of RNAi treatment for SF767 
#'
#' A dataset containing the
#' 
#' 
#' @format A data frame with 18 rows and 10 variables:
#' \describe{
#'   \item{CellLine}{Cell line ID}
#'   \item{Treatment}{Treatment}
#'   \item{GAPDH}{deltaCt of GAPDH}
#'   \item{CKAP2L}{deltaCt of CKAP2L}
#'   \item{ROCK1}{deltaCt of ROCK1}
#'   \item{TPR}{deltaCt of TPR}
#'   \item{ALDH4A1}{deltaCt of ALDH4A1}
#'   \item{GALNS}{deltaCt of GALNS}
#'   \item{CLCA2}{deltaCt of CLCA2}
#' }
#' 
"qPCR_SF767"

#' Illumina to REFSEQ_MRNA tabel
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
