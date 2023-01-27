# Loading libraries

library(CNAqc)
library(ggplot2)
library(dplyr)
library(ggrepel)


# Default functions for trimming and plotting data

# Trimming FACET output
trimming_cna = function(cna) {
  cna_trimmed = cna %>% select(chrom, start, end, tcn.em, lcn.em, Purity)
  colnames(cna_trimmed) = c("chr", "from", "to", "Major", "minor", "Purity")
  cna_trimmed$chrom = paste0("chr", cna_trimmed$chrom)
  cna_trimmed$CCF = 1
  cna_trimmed$Purity = NULL
  return(cna_trimmed)
}

# Getting purity information from FACET output
getting_purity = function(cna) {
  purity = 1
  if (is.na(cna$Purity[1]) != TRUE) {
    purity = cna$Purity[1]
  }  
  return(purity)
}

# Trimming maf file
trimming_mutations = function(mutations) {
  onco_mut_list = list("Missense_Mutation", "Frame_Shift_Del", "In_Frame_Del", "Nonsense_Mutation", "Splice_Site", "In_Frame_Ins")
  mutations = mutations[3:nrow(mutations),]
  mutations_trimmed = mutations %>% select(V1, V5, V6, V7, V9, V11, V13, V40, V42)
  colnames(mutations_trimmed) = c("gene", "chr", "from", "to", "mut_type", "ref", "alt", "DP", "NV")
  mutations_trimmed$VAF = as.numeric(mutations_trimmed$NV)/as.numeric(mutations_trimmed$DP)
  mutations_trimmed$chr = paste0("chr", mutations_trimmed$chr)
  mutations_trimmed$gene_name = paste0(mutations_trimmed$gene, ":", mutations_trimmed$mut_type)
  mutations_trimmed = mutations_trimmed %>% filter(mut_type %in% onco_mut_list)
  return(mutations_trimmed)
}

adding_mut_id = function(x) {
  x$mutations$mut_id = paste0(x$mutations$chr, ":", x$mutations$from, ":", x$mutations$to, ":", x$mutations$ref, ":", x$mutations$alt)
  return(x)
}

# Default plotting parameters

N = 5000
chromosomes = paste0('chr', c(1:22, 'X', 'Y'))
bl_plot = CNAqc:::blank_genome(ref = "hg38", chromosomes = chromosomes)
cols = c("Unique" = "red", "Shared" = "blue", "All_Shared" = "magenta", "PT_LN1_Shared"="orange", "PT_LN2_Shared" = "green", "LN1_LN2_Shared" = "blue")
shapes = c("Unique" = 21, "Shared" = 15, "All_Shared" = 15, "PT_LN1_Shared"= 16, "PT_LN2_Shared" = 17, "LN1_LN2_Shared" = 18)
top_shapes = c("high_in_both" = 15, "high_in_res" = 16, "high_in_non_res" = 17)

# Following function is modified from "plot_gw_vaf" from CNAqc package
plot_nonsyno_mut_VAF = function(x, name){
  cols = c("Unique" = "grey", "Shared" = "red", "All_Shared" = "red", "PT_LN1_Shared"="orange", "PT_LN2_Shared" = "green", "LN1_LN2_Shared" = "blue")
  shapes = c("Unique" = 20, "Shared" = 15, "All_Shared" = 15, "PT_LN1_Shared"= 16, "PT_LN2_Shared" = 17, "LN1_LN2_Shared" = 18)
  
  mutations = x$mutations %>% dplyr::filter(chr %in% chromosomes)
  mutations = CNAqc:::relative_to_absolute_coordinates(x,
                                                       mutations)
  # VAF stats
  quant = quantile(mutations$DP, probs = c(.1, .99))
  Q1_VAF = quantile(mutations$VAF, probs=0.25)
  Q3_VAF = quantile(mutations$VAF, probs=0.75)
  
  mutations = mutations %>%
    filter(DP > quant[1], DP < quant[2])
  
  if (nrow(mutations) > N)
    mutations = mutations %>% sample_n(N)
  
  # Drop out non-shared mutations
  mutations = mutations %>% dplyr::filter(is.na(shared) == FALSE)
  # Count each shared groups to mark on legend
  group_count = mutations %>% group_by(shared) %>% summarize(count = n())
  labels = paste0(group_count$shared, " (n=", group_count$count, ")")
  
  plot = bl_plot + geom_point(data = mutations,
                                       aes(x = from, y = VAF, shape = shared, color = shared),
                                       size = 3, alpha = 0.7) +
    geom_hline(yintercept = Q1_VAF, size = .4, # Q1 line
                        linetype = 'dashed',
                        color = 'darkred') +
    geom_hline(yintercept = Q3_VAF, size = .4, # Q3 line
                        linetype = 'dashed',
                        color = 'darkred') + 
    theme(axis.ticks.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text=element_text(size=18, face="bold"),
                   axis.title=element_text(size=18,face="bold"),
                   plot.title = element_text(size=26),
                   legend.text=element_text(size = 16, face = "bold")) +
    ylim(0,quantile(mutations$VAF, probs=0.98)) +
    labs(title = name, x= "Chromosome", y = "VAF") +
    guides(fill = "none") + 
    scale_color_manual(name = "", values = cols, labels = labels) + 
    scale_shape_manual(name = "", values = shapes, labels = labels)
  return(plot)
}
