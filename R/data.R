#' Overlapping classcodes
#'
#' classcodes used to classify assembled transfrags present in the tracking
#' files from the GffCompare utility were classified in overlapping
#' a reference transcript (in same strand) and non-overlapping
#'
#' @format A character vector of length 8
"overlapping_class_codes"

#' Non-overlapping classcodes
#'
#' classcodes used to classify assembled transfrags present in the tracking
#' files from the GffCompare utility were classified in overlapping
#' a reference transcript (in same strand) and non-overlapping
#'
#' @format A character vector of length 8
"non_overlapping_class_codes"

#' GTF
#'
#' Toy example of a mm39 mouse GTF file as imported with rtracklayer::readGFF\cr
#' NOTE that, when reading another GTF, not all variables
#' need to be equivalent to the ones in this example.\cr
#' Important columns are: "seqid", "start", "end", "strand",
#' "transcript_id", "gene_id"...\cr
#' This particular GTF example is a subset with 100 genes and 224 transcripts
#' form a merged transcriptome generated with gffCompare
#' from StringTie assemblies on 3 LSK total-RNA replicates
#'
#' @format A data.frame with 19 variables and 1901 observations
"gtf_df"
