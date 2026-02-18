####################
# title:"Watershed training weight plotting"
# author:"Sherry Yang"
# date:"2026-02-16"
####################

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(data.table)
library(optparse)
library(ggh4x)
library(forcats)

arguments <- parse_args(OptionParser(
  usage = "%prog [options]", description = "Watershed command line args",
  option_list = list(
    make_option(c("-m", "--model_object_rds"), default = NULL, help = "the rds model object"),
    make_option(c("-i", "--combined_annotation_input"), default = NULL, help = "csv file of the combined annoation from eval_prep script"),
    make_option(c("-o", "--output_file_png"), default = NULL, help = "output weight plotting graph name"),
    make_option(c("-t", "--omic_type"), default = NULL, help = "omic type input")
  )
))

model_object <- arguments$model_object_rds
combined_annotation <- arguments$combined_annotation_input
output <- arguments$output_file_png
omic <- arguments$omic_type

ws_obj <- readRDS(model_object)
inputs <- fread(combined_annotation)

colnames(inputs) <- tolower(colnames(inputs))

feats <- colnames(inputs)[3:(dim(inputs)[2] - 2)]
feat_df <- data.frame(feature = feats, theta = as.numeric(ws_obj$model_params$theta))
t1 <- feat_df %>%
  filter(is.finite(theta)) %>%
  arrange(desc(theta))
t1$feature <- gsub("unknown_", "", t1$feature, ignore.case = TRUE)
t1$feature <- trimws(t1$feature, which = "both")

# Lowercase all annotation type vectors for consistent matching
coding <- tolower(c("SV_5prime_exon_truncation", "coding_fraction_affected", "SV_3prime_exon_truncation", "exon_variant"))
conservation <- tolower(c("max_cpgpct", "top10_LINSIGHT", "mean_GC_content", "top10_phastCON", "top10_CADD"))
noncoding <- tolower(c(
  "nc_dist2gene", "tss_flank_impacted", "tes_flank_impacted",
  "downstream_noncoding_variant", "upstream_noncoding_variant", "is_ABC_SV"
))
regulatory <- tolower(c("num_enhancers_cell_types", "num_TADs", "remap_crm_score"))
SV <- tolower(c(
  "af", "SVTYPE_DEL_CNV", "Allele", "SVTYPE_INV", "SVTYPE_INS",
  "SVTYPE_DUP_CNV", "dCN", "SVTYPE_DUP", "SVTYPE_DEL", "length"
))
chromHMM <- tolower(c(
  "TxEnhW", "EnhA2", "PromBiv", "EnhW1", "TxEnh3", "Tx5", "TxEnh5", "PromD2",
  "EnhW2", "EnhA1", "ReprPC", "Quies", "EnhAc", "PromD1", "ZNF_Rpts", "DNase",
  "PromP", "TxReg", "Tx", "EnhAF", "PromU", "TxWk", "TssA", "Het", "Tx3"
))
vep <- tolower(c(
  "non_coding_transcript_variant", "intron_variant", "stop_lost", "frameshift_variant",
  "regulatory_region_ablation", "TFBS_ablation", "MODERATE", "NMD_transcript_variant",
  "feature_truncation", "feature_elongation", "start_retained_variant", "start_lost_HIGH",
  "splice_polypyrimidine_tract_variant", "LOW", "MODIFIER", "regulatory_region_variant",
  "TF_binding_site_variant", "non_coding_transcript_exon_variant", "3_prime_UTR_variant",
  "TFBS_amplification", "regulatory_region_amplification", "5_prime_UTR_variant"
))

matches_any <- function(feature, patterns) {
  sapply(feature, function(f) any(sapply(patterns, function(pat) grepl(pat, f, ignore.case = TRUE))))
}

t1$feature <- tolower(t1$feature)

t1$annotation_type <- case_when(
  matches_any(t1$feature, coding) ~ "coding",
  matches_any(t1$feature, conservation) ~ "conservation track",
  matches_any(t1$feature, chromHMM) ~ "encode chromhmm 25 states",
  matches_any(t1$feature, noncoding) ~ "noncoding",
  matches_any(t1$feature, regulatory) ~ "regulatory",
  matches_any(t1$feature, SV) ~ "sv",
  matches_any(t1$feature, vep) ~ "vep",
  TRUE ~ "other" # unmatched features are assigned to 'other'
)

t1$subcat <- case_when(
  grepl("tss", t1$feature, ignore.case = TRUE) ~ "tss",
  grepl("tes", t1$feature, ignore.case = TRUE) ~ "tes",
  grepl("gene_body", t1$feature, ignore.case = TRUE) ~ "gene_body",
  TRUE ~ "nonregion_specific"
)
t1$facet_group <- paste0(t1$annotation_type, "_", t1$subcat)


t1 <- t1 %>%
  group_by(annotation_type, subcat) %>%
  arrange(desc(theta), .by_group = TRUE) %>%
  mutate(feature = fct_inorder(feature)) %>%
  ungroup()

filename <- paste0(output, "_model_weights.tsv")
fwrite(t1, filename)
# Define colors (lowercase keys)
colors <- c(
  "coding" = "#e9602c",
  "conservation track" = "#eb346f",
  "encode chromhmm 25 states" = "#7845e3",
  "noncoding" = "#4a89f7",
  "regulatory" = "#298c8c",
  "sv" = "#9e1e63",
  "vep" = "#f2a125",
  "other" = "#3dd08178"
)

p <- ggplot(t1, aes(x = theta, y = feature, fill = annotation_type)) +
  geom_col() +
  facet_grid(rows = vars(annotation_type), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0, size = 6),
    axis.text.y = element_text(size = 5),
    legend.position = "none",
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(x = "Effect Size", y = "Annotations", title = paste0("Model Weights ", omic))


print(p)
ggsave(output, plot = p, width = 8, heigh = 6, dpi = 300)
