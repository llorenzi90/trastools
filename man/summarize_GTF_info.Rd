% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_to_extract_info_from_gtf.R
\name{summarize_GTF_info}
\alias{summarize_GTF_info}
\title{Summarize GTF Information}
\usage{
summarize_GTF_info(gtf_df, retrieve_exons = T, retrieve_introns = T)
}
\arguments{
\item{gtf_df}{A data frame containing GTF data. It should include columns such as \code{type}, \code{seqid}, \code{start}, \code{end}, \code{strand}, \code{transcript_id}, and \code{gene_id}.}

\item{retrieve_exons}{A logical value indicating whether to include detailed exon information in the output. Default is \code{TRUE}.}

\item{retrieve_introns}{A logical value indicating whether to include detailed intron information in the output. Default is \code{TRUE}.}
}
\value{
A list containing:
\item{summary}{A named vector with summary statistics including
number of transcripts, genes, exons, and introns,
as well as average and median exon and intron lengths, and other metrics.}
\item{EPT}{A data frame with the number of exons per transcript.}
\item{EPG}{A data frame with information on genes, including the maximum number of exons per gene and the type of gene (monoexonic or multiexonic).}
\item{exons_info}{A data frame with detailed information about exons (if \code{retrieve_exons} is \code{TRUE}).}
\item{introns_info}{A data frame with detailed information about introns (if \code{retrieve_introns} is \code{TRUE}).}
}
\description{
Summarizes information from a GTF (Gene Transfer Format) data frame,
including statistics about transcripts, genes, exons, and introns.
It calculates various metrics such as the number of transcripts and genes,
exon and intron lengths, and the ratio of exons to introns.
The function can also provide detailed information about exons and introns if requested.
}
\examples{
# Summarize GTF information
result <- summarize_GTF_info(gtf_df)
print(result$summary)
print(result$EPT)
print(result$EPG)
print(result$exons_info)
print(result$introns_info)
}
