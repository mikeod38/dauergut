SNPs <- read.delim(file.path(pathname,"data/WS254.CHII.NIL59.variants.gff2"), header = FALSE) %>% dplyr::select(c(1:5,9))
colnames(SNPs) <- c("Chromosome", "feature", "type", "start", "end", "desc")

### for max interval
SNPs %>% dplyr::filter(grepl("CB4856",desc)) %>% nrow ### count of SNPs
SNPs %>% dplyr::filter(grepl("CB4856",desc)) %>% dplyr::filter(grepl("AAChange",desc)) %>% nrow ### missense

# for min interval 
SNPs %>% dplyr::filter(grepl("CB4856",desc),start <= 13087000, start >= 12991000) %>% nrow ### count of SNPs
SNPs %>% dplyr::filter(grepl("CB4856",desc),start <= 13087000, start >= 12991000) %>% dplyr::filter(grepl("AAChange",desc)) %>% nrow ### missense
