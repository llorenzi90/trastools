## code to prepare `DATASET` dataset goes here
#overlapping_class_codes=c("=","c","e","o","j","k","m","n","p") # old one
#non_overlapping_class_codes=c("i","s","r","u","x","y",".") # old one

# It is better to not consider p as an overlap, because tecnically it is not

overlapping_class_codes=c("=","j","k","m","n","c","e","o")
non_overlapping_class_codes=c("p","i","s","u","r","x","y",".")
#
usethis::use_data(overlapping_class_codes,overwrite = T)
usethis::use_data(non_overlapping_class_codes,overwrite = T)
usethis::use_data(overlapping_class_codes, non_overlapping_class_codes, internal = TRUE,overwrite = T)

# GTF dataset
gtf_df=readRDS("inst/extdata/gtf.rds")
usethis::use_data(gtf_df,overwrite = T)

