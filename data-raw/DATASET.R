## code to prepare `DATASET` dataset goes here
overlapping_class_codes=c("=","c","e","o","j","k","m","n","p")
non_overlapping_class_codes=c("i","s","r","u","x","y",".")

usethis::use_data(overlapping_class_codes,overwrite = T)
usethis::use_data(non_overlapping_class_codes,overwrite = T)
