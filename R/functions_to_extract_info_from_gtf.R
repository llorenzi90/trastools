# Compute non-redundant exonic length of genes

#' Get per gene merged exons
#'
#' For each gene, overlapping regions
#' from its different exons are collapsed into
#' unique genomic regions.
#'
#' @param gtf_df A data frame containing the GTF (Gene Transfer Format) data.
#' This should include columns such as `seqid`, `start`, `end`, `strand`. \cr
#' GTF files can be loaded
#' with function `rtracklayer::readGFF`
#'
#' @return A data.frame with all non-redundant exonic regions covered by each gene
#' @export
#'
#' @examples
#' head(gtf_df)
#' get_per_gene_merged_exons(gtf_df)
get_per_gene_merged_exons <- function(gtf_df){

  gtf_df=gtf_df %>% dplyr::filter(type=="exon")

  txdb <-GenomicRanges::makeGRangesFromDataFrame(as.data.frame(gtf_df),
                                                 keep.extra.columns = T)

  per_gene_merged_exons <-
    as.data.frame(GenomicRanges::reduce(GenomicRanges::split(txdb,txdb$gene_id)))

  return(per_gene_merged_exons)
}

#' Compute non-redundant exonic length of genes
#'
#' For each gene, overlapping regions
#' from its different exons are first collapsed into
#' unique genomic regions. The 'non-redundant exonic
#' length' is then computed as the sum of the lengths
#' of the genomic regions covered by a gene's exonic
#' sequences.
#'
#'
#' @inheritParams get_per_gene_merged_exons
#'
#' @return A data.frame with gene_id and exonic_length columns
#' @export
#'
#' @examples
#'  head(gtf_df)
#'  compute_gene_non_redundant_exonic_length(gtf_df)
#'
compute_gene_non_redundant_exonic_length <- function(gtf_df){

  per_gene_merged_exons <- get_per_gene_merged_exons(gtf_df)

  gene_exonic_length <- stats::aggregate(per_gene_merged_exons$width,
                                  by=list(gene_id=per_gene_merged_exons$group_name),
                                  function(x)sum(x))
  colnames(gene_exonic_length)[2]="exonic_length"

  return(gene_exonic_length)
}



#' Calculates exonic length of each transcript
#'
#' @inheritParams get_per_gene_merged_exons
#'
#' @return a data.frame with transcript_id and length columns
#' @export
#'
#' @examples
#' compute_per_transcript_exonic_length(gtf_df)
compute_per_transcript_exonic_length <- function(gtf_df){
  gtf_df=gtf_df %>% dplyr::filter(type=="exon")
  gtf_df$lens=gtf_df$end - gtf_df$start+1
  tr_len=gtf_df%>%dplyr::group_by(transcript_id)%>%dplyr::summarise(len=sum(lens))
  return(tr_len)
}


# function to compute overlap with repeats including total gene fraction
# covered by repeats:

# get_overlap_with_repeats_from_GTF=function(gtf, repeats){
#   require(bedtoolsr)
#   require(tidyverse)
#   pgme=get_per_gene_merged_exons(gtf)
#   pgme=pgme[,c(3:5,2,6,7)]
#   pgme=pgme[order(pgme$seqnames,pgme$start),]
#   exonic_length=pgme%>%dplyr::group_by(group_name)%>%
#     dplyr::summarise(exonic_length=sum(width))
#   ol=bt.intersect(a = pgme,b=repeats,wo = T)
#   ol=ol%>% rowwise()%>% mutate(st=max(c(V2,V8)),en=min(c(V3,V9)))
#   olsm=ol%>%dplyr::group_by(V4)%>%
#     dplyr::summarise(repeats=paste(unique(V10[order(-V13)]),
#                             collapse = ","),
#               repeats_classes=paste(unique(V11[order(-V13)]),
#                                     collapse = ","),
#               t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
#                 start = st,end = en)))))
#   olsm=dplyr::left_join(olsm,exonic_length,by=c(V4="group_name"))
#   olsm$repeat.fraction=olsm$t.overlap_length/olsm$exonic_length
#   return(olsm)
# }


#' Get the number of exons per transcript or gene
#'
#' The functions `exons_per_transcript` and `exons_per_gene` retrieve
#' the exon count for each transcript or gene, respectively. At gene level,
#' the number of unique exons is retrieved, i.e. identical exons shared by different
#' transcripts are only counted once.
#' Additionally, an annotation column ('biotype' by default)
#' may be retrieved for each transcript or gene if present in then input gtf
#'
#' @inheritParams get_per_gene_merged_exons
#' @param anno_col An annotation column ('biotype' by default) you may want to retrieve together with the exon number
#'
#' @return A data.frame with transcript_id/gene_id, N_exons and 'anno_col' (if present in gtf)
#' @export
#'
#' @examples
#' exons_per_transcript(gtf_df)
exons_per_transcript <- function(gtf_df, anno_col="biotype"){
  ept <-  gtf_df%>%dplyr::filter(type=="exon")%>%
    dplyr::group_by(transcript_id)%>%
    dplyr::summarise(N_exons = dplyr::n())

  if(anno_col%in%colnames(gtf_df)){
    gtf_df=as.data.frame(gtf_df)
    ept$biotype=gtf_df[,anno_col][match(ept$transcript_id,
                                     gtf_df$transcript_id)]

    colnames(ept)[ncol(ept)]=anno_col
  }
  return(ept)

}

#' @rdname exons_per_transcript
#' @param gene_col Name of the column to use to summarise at gene level, default: "gene_id"
#'
exons_per_gene=function(gtf_df, gene_col="gene_id", anno_col="biotype"){
  gtf_df=as.data.frame(gtf_df)
  gtf_df$exon_id=paste0(gtf_df$seqid,":",
                        gtf_df$start,"-",
                        gtf_df$end,":",
                        gtf_df$strand)

  colnames(gtf_df)[colnames(gtf_df)==gene_col]="gene_id"
  epg <-  gtf_df%>%dplyr::filter(type=="exon")%>%dplyr::group_by(gene_id)%>%
    dplyr::summarise(N_exons=length(unique(exon_id)))

  if(anno_col%in%colnames(gtf_df)){
    gtf_df=as.data.frame(gtf_df)
    epg$biotype=gtf_df[,anno_col][match(epg$gene_id,
                                     gtf_df[,gene_col])]

    colnames(epg)[ncol(epg)]=anno_col
  }
  colnames(epg)[colnames(epg)=="gene_id"]=gene_col

  return(epg)

}


#' Get number of transcripts per gene
#'
#' @inheritParams exons_per_transcript
#' @param anno_col An annotation column ('biotype' by default)
#' you may want to retrieve together with the transcript number
#'
#'
#' @return A data.frame with gene_id, N_transcripts and
#' 'anno_col' (if present in gtf)
#' @export
#'
#' @examples
#' transcripts_per_gene(gtf_df)
transcripts_per_gene=function(gtf_df, gene_col="gene_id", anno_col="biotype"){
  gtf_df=as.data.frame(gtf_df)
  colnames(gtf_df)[colnames(gtf_df)==gene_col]="gene_id"
  tpg <-  gtf_df%>%dplyr::group_by(gene_id)%>%
    dplyr::summarise(N_transcripts=length(unique(transcript_id)))

  if(anno_col%in%colnames(gtf_df)){
    gtf_df=as.data.frame(gtf_df)
    tpg$biotype=gtf_df[,anno_col][match(tpg$gene_id,
                                     gtf_df[,gene_col])]

    colnames(tpg)[ncol(tpg)]=anno_col
  }
  colnames(tpg)[colnames(tpg)=="gene_id"]=gene_col

  return(tpg)

}


#' Get genomic range by gene
#'
#' Retrieves the genomic coordinates that each gene
#' spans from the start of the first (most 5') exon to the end of the last
#' (most 3') exon. The function returns a data frame in BED format with
#' the genomic ranges for each gene.
#'
#' @inheritParams get_per_gene_merged_exons
#' @param by A character string specifying the column name to group by.
#' Default is `"gene_id"`. This column should contain the gene identifiers.
#'
#' @return A data frame with the genomic ranges of each gene in BED format,
#'  containing columns: `seqid`, `start`, `end`, `gene_id`, `strand`, `width`.
#' @importFrom dplyr group_by mutate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#'
#' @examples
#' # Example GTF data frame
#' # Get genomic range by gene
#' result <- get_genomic_range_by_gene(gtf_df)
#' print(result)
get_genomic_range_by_gene <- function(gtf_df,by="gene_id"){
  df=as.data.frame(gtf_df)
  df$gene_id=df[,by]
  df=df[!is.na(df$gene_id),]
  df <- df %>% dplyr::group_by(gene_id) %>%
    dplyr::mutate(start=min(start),
                  end=max(end))
  df=df[!duplicated(df$gene_id),]

  genecoords=GenomicRanges::makeGRangesFromDataFrame(df[,c("seqid" , "start","end","strand","gene_id")],
                                                     keep.extra.columns = T)

  genecoords=as.data.frame(genecoords)
  colnames(genecoords)[6]=by
  #change format to bed file
  genecoords=genecoords[,c(1:3,6,4,5)]
  return(genecoords)
}

#' Summarize GTF Information
#'
#' Summarizes information from a GTF (Gene Transfer Format) data frame,
#' including statistics about transcripts, genes, exons, and introns.
#' It calculates various metrics such as the number of transcripts and genes,
#' exon and intron lengths, and the ratio of exons to introns.
#' The function can also provide detailed information about exons and introns if requested.
#'
#' @param gtf_df A data frame containing GTF data. It should include columns such as `type`, `seqid`, `start`, `end`, `strand`, `transcript_id`, and `gene_id`.
#' @param retrieve_exons A logical value indicating whether to include detailed exon information in the output. Default is `TRUE`.
#' @param retrieve_introns A logical value indicating whether to include detailed intron information in the output. Default is `TRUE`.
#'
#' @return A list containing:
#' \item{summary}{A named vector with summary statistics including
#' number of transcripts, genes, exons, and introns,
#' as well as average and median exon and intron lengths, and other metrics.}
#' \item{EPT}{A data frame with the number of exons per transcript.}
#' \item{EPG}{A data frame with information on genes, including the maximum number of exons per gene and the type of gene (monoexonic or multiexonic).}
#' \item{exons_info}{A data frame with detailed information about exons (if `retrieve_exons` is `TRUE`).}
#' \item{introns_info}{A data frame with detailed information about introns (if `retrieve_introns` is `TRUE`).}
#'
#' @importFrom dplyr filter group_by arrange reframe left_join
#' @importFrom stats median
#' @export
#'
#' @examples
#' # Summarize GTF information
#' result <- summarize_GTF_info(gtf_df)
#' print(result$summary)
#' print(result$EPT)
#' print(result$EPG)
#' print(result$exons_info)
#' print(result$introns_info)
summarize_GTF_info <- function(gtf_df, retrieve_exons=T,retrieve_introns=T){
  gtf_df_exons=gtf_df%>%dplyr::filter(type=="exon")
  N_tr=length(unique(gtf_df$transcript_id))
  N_genes=length(unique(gtf_df$gene_id))
  mean_tr_per_gene=N_tr/N_genes
  gtf_df_exons$exon_id=paste0(gtf_df_exons$seqid,":",
                           gtf_df_exons$start,"-",
                           gtf_df_exons$end,":",
                           gtf_df_exons$strand)
  N_exons=length(unique(gtf_df_exons$exon_id))

  # add exon number
  gtf_df_exonNum_plus=gtf_df_exons%>%dplyr::filter(strand=="+")%>%dplyr::group_by(transcript_id) %>%
    dplyr::reframe(exon_id,exon_number=1:dplyr::n())
  gtf_df_exonNum_minus=gtf_df_exons%>%dplyr::filter(strand=="-")%>%dplyr::group_by(transcript_id) %>%
    dplyr::reframe(exon_id,exon_number=dplyr::n():1)
  gtf_df_exonNum=rbind(gtf_df_exonNum_minus,gtf_df_exonNum_plus)
  gtf_df_exons$exon_number=gtf_df_exonNum$exon_number[match(gtf_df_exons$exon_id,
                                                      gtf_df_exonNum$exon_id)]

  # extract introns
  ttstart=gtf_df%>%dplyr::filter(type=="exon") %>% dplyr::group_by(transcript_id)%>%
    dplyr::arrange(start)%>%
    dplyr::reframe(start=end[-length(end)] + 1)
  ttend=gtf_df%>%dplyr::filter(type=="exon") %>% dplyr::group_by(transcript_id)%>%
    dplyr::arrange(start)%>%
    dplyr::reframe(end=start[-1] - 1)

  introns=ttstart
  introns$end=ttend$end

  intron_info=gtf_df_exons%>%dplyr::select(seqid,
                                        strand,
                                        transcript_id)
  intron_info=intron_info%>%dplyr::filter(!duplicated(intron_info))
  gtf_df_introns=dplyr::left_join(introns,intron_info)
  gtf_df_introns$intron_id=paste0(gtf_df_introns$seqid,":",gtf_df_introns$start,"-",
                               gtf_df_introns$end,":",gtf_df_introns$strand)
  gtf_df_intronNum_plus=gtf_df_introns%>%dplyr::filter(strand=="+")%>%dplyr::group_by(transcript_id) %>%
    dplyr::reframe(intron_id,intron_number=1:dplyr::n())
  gtf_df_intronNum_minus=gtf_df_introns%>%dplyr::filter(strand=="-")%>%dplyr::group_by(transcript_id) %>%
    dplyr::reframe(intron_id,intron_number=dplyr::n():1)
  gtf_df_intronNum=rbind(gtf_df_intronNum_minus,gtf_df_intronNum_plus)

  gtf_df_introns=dplyr::left_join(gtf_df_introns,gtf_df_intronNum)

  N_introns=length(unique(gtf_df_introns$intron_id))

  exon_intron_Ratio=N_exons/N_introns

  # exon and intron length
  gtf_df_exons$length=gtf_df_exons$end -gtf_df_exons$start +1
  summary(gtf_df_exons$length[!duplicated(gtf_df_exons$exon_id)])
  median_exon_length=stats::median(gtf_df_exons$length[!duplicated(gtf_df_exons$exon_id)])
  mean_exon_length=mean(gtf_df_exons$length[!duplicated(gtf_df_exons$exon_id)])
  gtf_df_introns$length=gtf_df_introns$end -gtf_df_introns$start +1
  summary(gtf_df_introns$length[!duplicated(gtf_df_introns$intron_id)])
  median_intron_length=stats::median(gtf_df_introns$length[!duplicated(gtf_df_introns$intron_id)])
  mean_intron_length=mean(gtf_df_introns$length[!duplicated(gtf_df_introns$intron_id)])

  # number of unassigned strand
  no_strand=sum(gtf_df$strand[!duplicated(gtf_df$transcript_id)]=="*")

  # mean exons per transcript
  EPT=exons_per_transcript(gtf = gtf_df)
  mean_EPT=mean(EPT$N_exons)

  # fraction monoexonic transcripts
  frac_monoexonic_tr=sum(EPT$N_exons==1)/nrow(EPT)

  # fraction monoexonic genes
  EPT$gene_id=gtf_df$gene_id[match(EPT$transcript_id,
                                gtf_df$transcript_id)]
  EPG=EPT%>%dplyr::group_by(gene_id)%>%dplyr::summarise(max_exons=max(N_exons),
                                          exonic_type=ifelse(any(N_exons!=1),"multiexonic","monoexonic"),
                                          any_mono=any(N_exons==1),
                                          N_tr=length(unique(transcript_id)))

  frac_monoexonic_gn=sum(EPG$exonic_type=="monoexonic")/nrow(EPG)

  median_tr_per_gene=stats::median(EPG$N_tr)
  # fraction of genes with monoexonic transcripts
  frac_genes_with_monoexonic_isoform=sum(EPG$any_mono)/nrow(EPG)
  # histogram N_max_exons vs N_genes
  #hist(EPG$max_exons)
  rlist = list(N_tr,N_genes,
             N_exons,
             mean_exon_length,
             median_exon_length,
             N_introns,
             mean_intron_length,
             median_intron_length,
             exon_intron_Ratio,
             mean_tr_per_gene,
             median_tr_per_gene,
             no_strand,
             mean_EPT,
             frac_monoexonic_tr,
             frac_monoexonic_gn,
             frac_genes_with_monoexonic_isoform

  )
  names(rlist)=c("N_transcripts",
                 "N_genes",
                 "N_exons",
                 "mean_exon_length",
                 "median_exon_length",
                 "N_introns",
                 "mean_intron_length",
                 "median_intron_length",
                 "exon_intron_Ratio",
                 "mean_tr_per_gene",
                 "median_tr_per_gene",
                 "transcript_with_no_strand",
                 "mean_exons_per_transcript",
                 "fraction_monoexonic_transcripts",
                 "fraction_monoexonic_genes",
                 "fraction_genes_with_monoexonic_isoform")
  rlist=list(rlist,EPT,EPG)
  names(rlist)=c("summary","EPT","EPG")

  if(retrieve_exons&retrieve_introns){
    rlist=c(rlist,list(gtf_df_exons,gtf_df_introns))
    names(rlist)=c("summary","EPT","EPG","exons_info","introns_info")

  }else{
      if(retrieve_exons){
        rlist=c(rlist,list(gtf_df_exons))
        names(rlist)=c("summary","EPT","EPG","exons_info")
      }else{
        if(retrieve_introns){
          rlist=c(rlist,list(gtf_df_introns))
          names(rlist)=c("summary","EPT","EPG","introns_info")
        }
      }
    }

  return(rlist)

}

#tlist=get_GTF_info(gtf)

# function idea ----
### annotate GTF ----
# the idea is to use the tracking file and the GTF and retrieve
# an annotated GTF data frame with all the analyses I want to include
# starting by
