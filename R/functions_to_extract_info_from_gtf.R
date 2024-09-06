# Compute non-redundant exonic length of genes

#' Get per gene merged exons
#'
#' For each gene, overlapping regions
#' from its different exons are collapsed into
#' unique genomic regions.
#'
#' @param gtf_df
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
#' @param gtf_df A data.frame resulting from reading a gtf file, e.g by
#'   rtracklayer::readGFF
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



compute_per_transcript_exonic_length <- function(gtf_df){
  gtf_df=gtf_df %>% dplyr::filter(type=="exon")
  gtf_df$lens=gtf_df$end - gtf_df$start+1
  tr_len=gtf_df%>%dplyr::group_by(transcript_id)%>%dplyr::summarise(len=sum(lens))
  return(tr_len)
}
get_per_gene_merged_transcripts <- function(gtf_df,gene_id_col="gene_id"){
  require(tidyverse)
  require(GenomicFeatures)
  gtf_df_exons=gtf_df %>% dplyr::filter(type=="transcript")

  txdb <-makeGRangesFromDataFrame(as.data.frame(gtf_df_exons),keep.extra.columns = T)

  per_gene_merged_transcripts <- as.data.frame(IRanges::reduce(split(txdb,txdb$gene_id)))

  return(per_gene_merged_transcripts)
}
# function to compute overlap with repeats including total gene fraction
# covered by repeats:

get_overlap_with_repeats_from_GTF=function(gtf){
  require(bedtoolsr)
  require(tidyverse)
  pgme=get_per_gene_merged_exons(gtf)
  pgme=pgme[,c(3:5,2,6,7)]
  pgme=pgme[order(pgme$seqnames,pgme$start),]
  exonic_length=pgme%>%dplyr::group_by(group_name)%>%
    dplyr::summarise(exonic_length=sum(width))
  ol=bt.intersect(a = pgme,b=repeats,wo = T)
  ol=ol%>% rowwise()%>% mutate(st=max(c(V2,V8)),en=min(c(V3,V9)))
  olsm=ol%>%dplyr::group_by(V4)%>%
    dplyr::summarise(repeats=paste(unique(V10[order(-V13)]),
                            collapse = ","),
              repeats_classes=paste(unique(V11[order(-V13)]),
                                    collapse = ","),
              t.overlap_length=sum(IRanges::width(IRanges::reduce(IRanges(
                start = st,end = en)))))
  olsm=left_join(olsm,exonic_length,by=c(V4="group_name"))
  olsm$repeat.fraction=olsm$t.overlap_length/olsm$exonic_length
  return(olsm)
}


exons_per_transcript=function(gtf, anno_col="biotype"){
  ept <-  gtf%>%dplyr::filter(type=="exon")%>%dplyr::group_by(transcript_id)%>%
    dplyr::summarise(N_exons=n())

  if(anno_col%in%colnames(gtf)){
    gtf=as.data.frame(gtf)
    ept$biotype=gtf[,anno_col][match(ept$transcript_id,
                                     gtf$transcript_id)]

    colnames(ept)[ncol(ept)]=anno_col
  }
  return(ept)

}


exons_per_gene=function(gtf, gene_col="gene_id", anno_col="biotype"){
  gtf$exon_id=paste0(gtf$seqid,":",gtf$start,"-",gtf$end,":",gtf$strand)
  gtf=as.data.frame(gtf)
  colnames(gtf)[colnames(gtf)==gene_col]="gene_id"
  epg <-  gtf%>%dplyr::filter(type=="exon")%>%dplyr::group_by(gene_id)%>%
    dplyr::summarise(N_exons=length(unique(exon_id)))

  if(anno_col%in%colnames(gtf)){
    gtf=as.data.frame(gtf)
    epg$biotype=gtf[,anno_col][match(epg$gene_id,
                                     gtf[,gene_col])]

    colnames(epg)[ncol(epg)]=anno_col
  }
  colnames(epg)[colnames(epg)=="gene_id"]=gene_col

  return(epg)

}


transcripts_per_gene=function(gtf, gene_col="gene_id", anno_col="biotype"){
  gtf=as.data.frame(gtf)
  colnames(gtf)[colnames(gtf)==gene_col]="gene_id"
  tpg <-  gtf%>%dplyr::group_by(gene_id)%>%
    dplyr::summarise(N_transcripts=length(unique(transcript_id)))

  if(anno_col%in%colnames(gtf)){
    gtf=as.data.frame(gtf)
    tpg$biotype=gtf[,anno_col][match(tpg$gene_id,
                                     gtf[,gene_col])]

    colnames(tpg)[ncol(tpg)]=anno_col
  }
  colnames(tpg)[colnames(tpg)=="gene_id"]=gene_col

  return(tpg)

}


require(data.table)
require(GenomicRanges)
get_genomic_regions=function(gffdf,by="gene_id"){
  df=as.data.frame(gffdf)
  df$gene_id=df[,by]
  df=df[!is.na(df$gene_id),]
  df=as.data.table(df)
  df[,':='(start=min(start),end=max(end)),by=gene_id]
  df=df[!duplicated(df$gene_id),]
  genecoords=makeGRangesFromDataFrame(df[,c("seqid" , "start","end","strand","gene_id")],keep.extra.columns = T)
  genecoords=as.data.frame(genecoords)
  colnames(genecoords)[6]=by
  #change format to bed file
  genecoords=genecoords[,c(1:3,6,4,5)]
  return(genecoords)
}

get_GTF_info <- function(GTF, retrieve_exons=T,retrieve_introns=T){
  require(tidyverse)
  GTF_exons=GTF%>%dplyr::filter(type=="exon")
  N_tr=length(unique(GTF$transcript_id))
  N_genes=length(unique(GTF$gene_id))
  mean_tr_per_gene=N_tr/N_genes
  GTF_exons$exon_id=paste0(GTF_exons$seqid,":",
                           GTF_exons$start,"-",
                           GTF_exons$end,":",
                           GTF_exons$strand)
  N_exons=length(unique(GTF_exons$exon_id))

  # add exon number
  GTF_exonNum_plus=GTF_exons%>%dplyr::filter(strand=="+")%>%dplyr::group_by(transcript_id) %>%
    reframe(exon_id,exon_number=1:n())
  GTF_exonNum_minus=GTF_exons%>%dplyr::filter(strand=="-")%>%dplyr::group_by(transcript_id) %>%
    reframe(exon_id,exon_number=n():1)
  GTF_exonNum=rbind(GTF_exonNum_minus,GTF_exonNum_plus)
  GTF_exons$exon_number=GTF_exonNum$exon_number[match(GTF_exons$exon_id,
                                                      GTF_exonNum$exon_id)]

  # extract introns
  ttstart=GTF%>%dplyr::filter(type=="exon") %>% dplyr::group_by(transcript_id)%>%
    arrange(start)%>%
    reframe(start=end[-length(end)] + 1)
  ttend=GTF%>%dplyr::filter(type=="exon") %>% dplyr::group_by(transcript_id)%>%
    arrange(start)%>%
    reframe(end=start[-1] - 1)

  introns=ttstart
  introns$end=ttend$end

  intron_info=GTF_exons%>%dplyr::select(seqid,
                                        strand,
                                        transcript_id)
  intron_info=intron_info%>%dplyr::filter(!duplicated(intron_info))
  GTF_introns=left_join(introns,intron_info)
  GTF_introns$intron_id=paste0(GTF_introns$seqid,":",GTF_introns$start,"-",
                               GTF_introns$end,":",GTF_introns$strand)
  GTF_intronNum_plus=GTF_introns%>%dplyr::filter(strand=="+")%>%dplyr::group_by(transcript_id) %>%
    reframe(intron_id,intron_number=1:n())
  GTF_intronNum_minus=GTF_introns%>%dplyr::filter(strand=="-")%>%dplyr::group_by(transcript_id) %>%
    reframe(intron_id,intron_number=n():1)
  GTF_intronNum=rbind(GTF_intronNum_minus,GTF_intronNum_plus)

  GTF_introns=left_join(GTF_introns,GTF_intronNum)

  N_introns=length(unique(GTF_introns$intron_id))

  exon_intron_Ratio=N_exons/N_introns

  # exon and intron length
  GTF_exons$length=GTF_exons$end -GTF_exons$start +1
  summary(GTF_exons$length[!duplicated(GTF_exons$exon_id)])
  median_exon_length=median(GTF_exons$length[!duplicated(GTF_exons$exon_id)])
  mean_exon_length=mean(GTF_exons$length[!duplicated(GTF_exons$exon_id)])
  GTF_introns$length=GTF_introns$end -GTF_introns$start +1
  summary(GTF_introns$length[!duplicated(GTF_introns$intron_id)])
  median_intron_length=median(GTF_introns$length[!duplicated(GTF_introns$intron_id)])
  mean_intron_length=mean(GTF_introns$length[!duplicated(GTF_introns$intron_id)])

  # number of unassigned strand
  no_strand=sum(GTF$strand[!duplicated(GTF$transcript_id)]=="*")

  # mean exons per transcript
  EPT=exons_per_transcript(gtf = GTF)
  mean_EPT=mean(EPT$N_exons)

  # fraction monoexonic transcripts
  frac_monoexonic_tr=sum(EPT$N_exons==1)/nrow(EPT)

  # fraction monoexonic genes
  EPT$gene_id=GTF$gene_id[match(EPT$transcript_id,
                                GTF$transcript_id)]
  EPG=EPT%>%dplyr::group_by(gene_id)%>%dplyr::summarise(max_exons=max(N_exons),
                                          exonic_type=ifelse(any(N_exons!=1),"multiexonic","monoexonic"),
                                          any_mono=any(N_exons==1),
                                          N_tr=length(unique(transcript_id)))

  frac_monoexonic_gn=sum(EPG$exonic_type=="monoexonic")/nrow(EPG)

  median_tr_per_gene=median(EPG$N_tr)
  # fraction of genes with monoexonic transcripts
  frac_genes_with_monoexonic_isoform=sum(EPG$any_mono)/nrow(EPG)
  # histogram N_max_exons vs N_genes
  hist(EPG$max_exons)
  rlist=list(N_tr,N_genes,
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

  # if(retrieve_exons&retrieve_introns){
  #   rlist=list(rlist,EPT,EPG,GTF_exons,GTF_introns)
  # }else{
  #     if(retrieve_exons){
  #       rlist=list(rlist,EPT,EPG,GTF_exons)
  #
  #     }else{
  #       if(retrieve_introns){
  #         rlist=list(rlist,EPT,EPG,GTF_introns)
  #
  #       }
  #     }
  #   }

  rlist=list(rlist,EPT,EPG,GTF_exons,GTF_introns)
  names(rlist)=c("summary","EPT","EPG","exons_info","introns_info")
  return(rlist)

}

#tlist=get_GTF_info(gtf)

