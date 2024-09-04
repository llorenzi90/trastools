#' @import utils
utils::globalVariables(c("V2", "V4", "Nsamples"))

#' Get number of samples in which each transcript is assembled
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#'
#' @return A vector with the number of samples for each transcript
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_n_samples_per_transcript_from_tracking(tracking)
get_n_samples_per_transcript_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  apply(tracking[,cols],1,
        function(x)sum(x!="-"))

}

#' Get number of exons for each assembled transcript
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#'
#' @return A vector with the number of exons for each transcript
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_n_exons_per_transcript(tracking)
get_n_exons_per_transcript <- function(tracking,cols=c(5:ncol(tracking))){

  single_tr_info=apply(tracking[,cols],1,function(x)x[x!="-"][1])
  n_exons=sapply(strsplit(single_tr_info,split = "\\|"),function(x)x[3])
  return(n_exons)
}


#' Split samples information columns into gene_id and transcrip_id
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#' @param qnames Optional character vector with sample names
#'
#' @return Modified tracking table with extra columns from separated gene and transcript for each sample
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' split_samples_info(tracking)
split_samples_info <- function(tracking,cols=c(5:ncol(tracking)),qnames=NULL){
  if(is.null(qnames)){
    qnames=paste0("q",1:length(cols))
  }
  names(cols)=qnames
  c1=cols[1]
  cc=c1
  for (qn in names(cols)) {

    tracking=tidyr::separate(tracking,col = cc,into = c(paste0(qn,"_gene"),paste0(qn,"_transcript")),
                      sep = "\\|",extra = "drop",fill = "right")
    tracking[,paste0(qn,"_gene")]=gsub("q[1-9*]:","",tracking[,paste0(qn,"_gene")])
    cc=cc+2
  }
  return(tracking)
}

# gene level functions ----

#' Get sample occurrance per gene
#'
#' For each 'XLOC...' loci returns for each sample TRUE/FALSE if it is
#' present or not in that sample
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#'
#' @return A table with a logical column for each sample and gene_id
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_sample_occurrance_per_gene_from_tracking(tracking)
get_sample_occurrance_per_gene_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  sample_occurrance_per_gene=lapply(cols, function(sa){
    stats::aggregate(tracking[,sa],by=list(gene_id=tracking$V2),function(x)any(x!="-"))
  })
  gene_id=sample_occurrance_per_gene[[1]]$gene_id

  sample_occurrance_per_gene <- as.data.frame(sapply(sample_occurrance_per_gene,function(x)x$x))
  sample_occurrance_per_gene$gene_id=gene_id

  return(sample_occurrance_per_gene)
}


#' Get number of samples per gene from tracking
#'
#' For each 'XLOC...' loci returns the number of samples in which it appears
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#'
#' @return A named integer vector with the number of samples in which each gene_id appears
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_n_samples_per_gene_from_tracking(tracking)
get_n_samples_per_gene_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  sopg=get_sample_occurrance_per_gene_from_tracking(tracking,cols)
  Nsamps=rowSums(sopg[,1:length(cols)])
  names(Nsamps)=sopg$gene_id
  return(Nsamps)
}

#' Get maximum number of samples in which individual transcripts from each gene are assembled
#'
#' Retrieves the maximum number of samples in which any transcript from each
#' 'XLOC...' gene is assembled. This gives an idea of how reproducible each
#' merged loci is between samples.
#'
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#'
#' @return A table with gene_id and the maximum number of samples in which individual transcripts from each gene are assembled
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_Max_n_samples_per_gene_from_tracking(tracking)
get_Max_n_samples_per_gene_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  tracking$Nsamples=get_n_samples_per_transcript_from_tracking(tracking,cols)
  Nsamps=tracking%>%dplyr::group_by(V2) %>%dplyr::summarise(maxNsamps=max(Nsamples))
  colnames(Nsamps)[1]="gene_id"
  return(Nsamps)
}

#' Check and report overlap with reference transcripts at gene level
#'
#' Overlap is defined first at transcript level. If the classcode is one of \code{overlapping_class_codes}
#' If any transcript in an 'XLOC...' loci is an overlap, then the gene is an overlap.
#'
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param overlapping_classcodes Character vector, classcodes to consider as overlapping, default: "=" "c" "e" "o" "j" "k" "m" "n" "p"
#' @return A table with gene_id and logical for overlap with reference transcripts
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_overlapRef_gene_level(tracking)
get_overlapRef_gene_level <- function(tracking,overlapping_classcodes = overlapping_class_codes){
  oref=tracking$V4%in%overlapping_classcodes
  orefGL=stats::aggregate(oref,by=list(gene_id=tracking$V2),function(x)any(x))
  return(orefGL)
}

# extract expression values ----

#' Get expression values in TPM or FPKM from gffCompare tracking
#'
#' Retrieves a table with the desired expression values (TPM default, or FPKM) from each sample's information column in the provided
#' tracking file
#'
#' @param tracking A gffCompare tracking file-like data.frame that matches transcripts up between samples
#' @param measure Character: "tpm" (default) or "fpkm"
#' @param cols Column numbers corresponding to samples in the tracking file, defaults to 5:ncol(tracking)
#'
#' @return Table of expression values ('NA' if transcript not present in sample)
#' with one column per sample and one column with the mean values across samples ('NA's removed)
#' @export
#'
#' @examples
#' tracking <- structure(list(
#' V1 = c("TCONS_00000001", "TCONS_00000002", "TCONS_00000003",
#' "TCONS_00000004", "TCONS_00000005", "TCONS_00000006"),
#' V2 = c("XLOC_000001",
#'        "XLOC_000003",
#'        "XLOC_000004",
#'        "XLOC_000005",
#'        "XLOC_000009",
#'        "XLOC_000009"
#' ),
#' V3 = c("XLOC_002036|TCONS_00004598", "-", "-", "XLOC_002044|TCONS_00004614",
#'           "XLOC_000074|TCONS_00000075", "XLOC_000074|TCONS_00000080"),
#' V4 = c("i", "u", "r", "x", "=", "="),
#' V5 = c("q1:STRG.1|STRG.1.1|1|0.091235|0.267011|6.848781|410",
#'        "q1:STRG.5|STRG.5.1|1|0.066681|0.195150|5.005556|360",
#'        "q1:STRG.10|STRG.10.1|1|0.084226|0.246497|6.322581|248",
#'        "q1:STRG.13|STRG.13.1|1|0.079929|0.233920|6.000000|251",
#'        "q1:STRG.32|STRG.32.1|9|8.963122|26.231630|672.835327|2474",
#'        "q1:STRG.32|STRG.32.2|8|0.459790|1.345631|34.515125|2372"
#'        ),
#' V6 = c("-", "-", "-", "-",
#' "q2:STRG.11|STRG.11.1|9|8.893389|29.440718|746.519531|2430",
#' "q2:STRG.11|STRG.11.2|8|0.305824|1.012402|25.671185|2328"
#' ),
#' V7 = c("-", "-", "-", "-",
#' "q3:STRG.23|STRG.23.2|9|8.601234|23.619555|632.080566|2429",
#' "q3:STRG.23|STRG.23.1|8|0.289132|0.793974|21.247471|2594"
#' )),
#' row.names = c(NA, 6L), class = "data.frame")
#' get_expression_values(tracking)
get_expression_values <- function(tracking, measure=c("tpm","fpkm"),cols=5:ncol(tracking)){
  measure=measure[1]
  print(measure)
  targetcol=ifelse(measure=="tpm",5,ifelse(measure=="fpkm",4,
                                           stop("expression measure must be 'tpm' or 'fpkm' ")))
  exp_mat=apply(tracking[,cols],2,function(s)sapply(strsplit(s,split = "\\|"),function(x)x[targetcol]))
  rownames(exp_mat)=tracking$V1
  mean_expr=apply(exp_mat,1, function(x)mean(as.numeric(x[!is.na(x)])))
  exp_mat=as.data.frame(exp_mat)
  colnames(exp_mat) <- paste0("samp_",colnames(exp_mat))
  exp_mat$mean_expr=mean_expr
  colnames(exp_mat)[ncol(exp_mat)]=paste0("mean_",measure)
  return(exp_mat)
}
