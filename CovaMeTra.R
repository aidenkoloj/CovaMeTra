library(tidyverse)
library(readr)
library(dplyr)
library(readxl)
library(lsa)
library(reshape2)

### Load in data -------------------------------------------------------------------------


metabolites = melt(read.csv("NIDDK_soup_neg_11_02.csv"))
rna_seq = read_csv("NIDDK_RNA_seq_with_replicate_ID.csv")

### Combine and tidy data -------------------------------------------------------------------------

metabolites = metabolites %>% rename(feature = variable, abundance = value, replicate = Replicate) %>% select(replicate, feature, abundance)
metabolites$feature <- paste(metabolites$feature, sep="_","metabolite")

rna_seq = rna_seq %>% rename(feature = ext_gene, abundance = normalized_counts) %>% select(replicate, feature, abundance)

metra_data <- rbind(metabolites, rna_seq)


# CovaMeTra -------------------------------------------------------------------------
#
## Takes 'x' a gene name or feature name and returns highest correlated genes/features. 
## Metabolite features must have "_metabolite" on the end 
#
#-------------------------------------------------------------------------

covametra = function(x){
  
  input_data <- 
    metra_data %>% 
    #Select columns of interest
    filter(feature == x) 
  
  df <- data.frame(a = c(input_data$feature), b = c(input_data$abundance)) 
  df <- bind_rows(replicate((nrow(metra_data)/24), df, simplify = FALSE))
  combined <- cbind(metra_data, df) 
  combined_r2 <- combined %>% rename(input_names = a, input_abundance = b) %>% group_by(feature) %>% 
    summarise(r_2 = cor(abundance,input_abundance), method = "pearson") 
  combined_r2$r_2 = combined_r2$r_2*combined_r2$r_2
  combined_r2 <- combined_r2 %>% 
    arrange(desc(r_2)) 
  return(combined_r2)
  
  
}
#-------------------------------------------------------------------------
# To search for covariance of genes, type a gene name, to search for metabolites
# include "_metabolite" after the feature. 

results = covametra("daf-22")

