# Functions for microbial ecology

library("phyloseq")
library("microbiome")
library("vegan")
library("tidyverse")
library("ggplot2")
library("plotrix")
library("iNEXT")
library("purrr")
library("EcolUtils")
library("geosphere")
library("nlme")
library("RColorBrewer")
library('reshape2')
library('stringr')
library("DESeq2")
library("pairwiseAdonis")
library("ggvenn")
library("lme4")
library('emmeans')
library('decontam')
library('ggordiplots')
library('ANCOMBC')
library("ggdendro")
library("car")
library("qiime2R")
library("ggpattern")
library("cowplot")
library("metagMisc")
library("betareg")

nonbact_filter <- function(phylo_data) {
bdf <- subset_taxa(phylo_data, (Order!="Chloroplast") | is.na(Order)) # 48 taxa
bdf2 <- subset_taxa(bdf, (Family!="Mitochondria") | is.na(Family)) # 41 taxa
bdf3 <- subset_taxa(bdf2, (Kingdom != "Archaea") | is.na(Kingdom)) # 23 taxa
bdf4 <- subset_taxa(bdf3, !is.na(Phylum)) # 34 taxa
}

prevalence_table <- function(phylo_data) {
  prevalence_table_sample = apply(X = otu_table(phylo_data),
                                MARGIN = 1,
                                FUN = function(x){sum(x > 0)})
  
  prevalence_table_sample2 = data.frame(Prevalence = prevalence_table_sample,
                                     TotalAbundance = taxa_sums(phylo_data),
                                     tax_table(phylo_data))
  
  prevalence_table_sample2_sort <- prevalence_table_sample2 %>% arrange(desc(Prevalence))
  
  return(prevalence_table_sample2_sort)
}



low_read_filter <- function(phylo_data, prev_table, read_depth){
  
  filter_otus <- prev_table%>% 
    filter(TotalAbundance > read_depth) %>%
    dplyr::select(OTU)
  filter_otus <- as.vector(filter_otus[,1])
  
  data2 <- prune_taxa(filter_otus, phylo_data)
  
  return(data2)
}

pcoa_plotter <- function(dist_matrix, samp_dat, cent_variable, join_cent, shape_var){

  beta_dis <- betadisper(d = dist_matrix, group = cent_variable)
  b2 <- permutest(beta_dis)
  s1 <- scores(beta_dis)
  
 #extract scores and join with sample_data
  sample.wcr <- tibble::rownames_to_column(samp_dat, "sample")
  scores_pcoa_wcr <- as.data.frame(s1$sites)
  scores_pcoa_wcr2 <- tibble::rownames_to_column(scores_pcoa_wcr, "sample")
  scores_pcoa_wcr3 <- inner_join(scores_pcoa_wcr2, sample.wcr, by = "sample")
  
  centroids_wcr <- scores(beta_dis, display = "centroids")
  centroids_wcr2 <- as.data.frame(centroids_wcr)
  centroids_wcr3 <- tibble::rownames_to_column(centroids_wcr2, join_cent)
  
  plot_data_wcr <- inner_join(scores_pcoa_wcr3, centroids_wcr3, by = join_cent)
  
  rare_wcr_figure <- plot_data_wcr %>% 
    ggplot(aes(x = PCoA1.x, y = PCoA2.x, colour = .data[[join_cent]], shape = .data[[shape_var]])) +
    geom_point() +
    geom_segment(aes(x = PCoA1.x, y = PCoA2.x, xend = PCoA1.y, yend = PCoA2.y)) +
    theme_classic()

  return(list(rare_wcr_figure, b2))

}


prev_tables <- function(phylobj){
  prevelance_table = apply(X = otu_table(phylobj),
                                    MARGIN = 1,
                                    FUN = function(x){sum(x > 0)})
  
  prevelance_table = data.frame(Prevalence = prevelance_table,
                                         TotalAbundance = taxa_sums(phylobj),
                                         tax_table(phylobj))
  
  # mean abundance of phylum
  phylum_abundance <- plyr::ddply(prevelance_table, "Phylum", function(df1){
    data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
  
  percent_phylum_abundance <- phylum_abundance %>% 
    mutate(percent_total = (total_abundance / sum(colSums(otu_table(phylobj)))) * 100) %>% 
    dplyr::arrange(desc(percent_total))
  
  #mean abundance of class
  class_abundance <- plyr::ddply(prevelance_table, "Class", function(df1){
    data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
  
  percent_class_abundance <- class_abundance %>% 
    mutate(percent_total = (total_abundance / sum(colSums(otu_table(phylobj)))) * 100) %>% 
    dplyr::arrange(desc(percent_total))
  
  # mean abundance of order
  order_abundance <- plyr::ddply(prevelance_table, "Order", function(df1){
    data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
  
  percent_order_abundance <- order_abundance %>% 
    mutate(percent_total = (total_abundance / sum(colSums(otu_table(phylobj)))) * 100) %>% 
    dplyr::arrange(desc(percent_total))
  
  # mean abundance of family
  family_abundance <- plyr::ddply(prevelance_table, "Family", function(df1){
    data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
  
  percent_family_abundance <- family_abundance %>% 
    mutate(percent_total = (total_abundance / sum(colSums(otu_table(phylobj)))) * 100) %>% 
    dplyr::arrange(desc(percent_total))
  
  # mean abundance of genus
  genus_abundance  <- plyr::ddply(prevelance_table, "Genus", function(df1){
    data.frame(mean_prevalence=mean(df1$Prevalence), total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
  
  percent_genus_abundance <- genus_abundance %>% 
    mutate(percent_total = (total_abundance / sum(colSums(otu_table(phylobj)))) * 100) %>% 
    dplyr::arrange(desc(percent_total))
  
  percent_abundance <- list(percent_phylum_abundance, percent_class_abundance, percent_order_abundance, percent_family_abundance, percent_genus_abundance)
  
  return(percent_abundance)
  
}

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
