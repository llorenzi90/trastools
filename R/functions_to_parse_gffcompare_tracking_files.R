#' @import utils
utils::globalVariables(c("V2", "V4", "Nsamples"))
get_n_samples_per_transcript_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  apply(tracking[,cols],1,
        function(x)sum(x!="-"))

}

get_n_exons_per_transcript <- function(tracking,cols=c(5:ncol(tracking))){

  single_tr_info=apply(tracking[,cols],1,function(x)x[x!="-"][1])
  n_exons=sapply(strsplit(single_tr_info,split = "\\|"),function(x)x[3])
  return(n_exons)
}


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
get_sample_occurrance_per_gene_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  sample_occurrance_per_gene=lapply(cols, function(sa){
    stats::aggregate(tracking[,sa],by=list(gene_id=tracking$V2),function(x)any(x!="-"))
  })
  gene_id=sample_occurrance_per_gene[[1]]$gene_id

  sample_occurrance_per_gene <- as.data.frame(sapply(sample_occurrance_per_gene,function(x)x$x))
  sample_occurrance_per_gene$gene_id=gene_id

  return(sample_occurrance_per_gene)
}


get_n_samples_per_gene_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  sopg=get_sample_occurrance_per_gene_from_tracking(tracking,cols)
  Nsamps=rowSums(sopg[,1:length(cols)])
  names(Nsamps)=sopg$gene_id
  return(Nsamps)
}

get_Max_n_samples_per_gene_from_tracking <- function(tracking,cols=c(5:ncol(tracking))){
  tracking$Nsamples=get_n_samples_per_transcript_from_tracking(tracking,cols)
  Nsamps=tracking%>%dplyr::group_by(V2) %>%dplyr::summarise(maxNsamps=max(Nsamples))
  colnames(Nsamps)[1]="gene_id"
  return(Nsamps)
}
#test1=split_samples_info(tracking = merged_refs_track[1:100,])
#single_tr_info=apply(tracking_SL[,6:8],1,function(x)x[x!="-"][1])
#n_exons=sapply(strsplit(single_tr_info,split = "\\|"),function(x)x[3])
get_overlapRef_gene_level <- function(tracking,overlapping_classcodes = overlapping_class_codes){
  oref=tracking$V4%in%overlapping_classcodes
  orefGL=stats::aggregate(oref,by=list(gene_id=tracking$V2),function(x)any(x))
  return(orefGL)
}

# extract expression values ----

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
