####################
# title:"multi-ome analysis: AD/PD risk gene colocalization with Watershed variants"
# author:"Sherry Yang"
# date:"2026-02-16"
####################

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(biomaRt)
library(tidyverse)
library(patchwork)
library(ggplot2)


## Figure 3d/e
sample_swap <- fread("~/syang/ADRC/meta_data/ADRC_SAMS.metadata.master_table.sample_swaps_annotated.tsv")
sample_swap <- sample_swap[which(sample_swap$DROP_SAMPLE == FALSE)]
rare_variant_info_expression <- fread("~/syang/Watershed-SV/run_watershed_combinations/rna_105/combined_annotation_pre_merge_noimpute.medZ.csv")
rare_variant_info_expression$gene_id <- gsub(":.*", "", rare_variant_info_expression$GeneName)
rare_variant_info_expression$sniffle <- sub(".*:", "", rare_variant_info_expression$GeneName)


rare_variant_info_prot <- fread("~/syang/Watershed-SV/run_watershed_combinations/prot_25/unknown/combined_annotation_pre_merge_noimpute.medZ.csv")
rare_variant_info_prot$gene_id <- gsub(":.*", "", rare_variant_info_prot$GeneName)
rare_variant_info_prot$sniffle <- sub(".*:", "", rare_variant_info_prot$GeneName)
rare_variant_info_prot

rare_variant_info_csf_prot <- fread("~/syang/Watershed-SV/run_watershed_combinations/csf_pro_25/combined_annotation_pre_merge_noimpute.medZ.csv")
rare_variant_info_csf_prot$gene_id <- gsub(":.*", "", rare_variant_info_csf_prot$GeneName)
rare_variant_info_csf_prot$sniffle <- sub(".*:", "", rare_variant_info_csf_prot$GeneName)


### loading files
coloc_expression <- fread("~/tannerj/ADRComics/QTL_mapping/ADRC_coloc/brain_gwas/output/coloc/2025-01-28_11-57-33.462186_GTEx_coloc/GTEx_coloc.AD_PD_gene_lists.tsv")
coloc_expression$feature <- sub("\\..*", "", coloc_expression$feature)
coloc_expression_0.5_conf <- coloc_expression[which(coloc_expression$AD_max_pp4 >= 0.5 | coloc_expression$PD_max_pp4 >= 0.5), ]

### RNA
expression_prioritized_variants <- fread("~/syang/Watershed-SV/Watershed_strobe/model_RIVER_number_of_dimensions_pc_105_1_posterior_probability.txt")
expression_with_zscore <- fread("~/syang/ADRC/copy_unprocessed_data/output_filtered_rna_pc_105.csv")
expression_prioritized_variants$ID <- sub(".*:", "", expression_prioritized_variants$sample_names)
expression_prioritized_variants$PIDN <- as.integer(gsub(":.*", "", expression_prioritized_variants$sample_names))
expression_prioritized_0.6_conf <- expression_prioritized_variants[which(expression_prioritized_variants$Watershed_posterior_outlier_signal_1 >= 0.6), ]


coloc_expression_watershed_prioritized <- merge(coloc_expression_0.5_conf, expression_prioritized_0.6_conf, by.x = "feature", by.y = "ID")
expression_merged_by_patient <- coloc_expression_watershed_prioritized %>%
  left_join(sample_swap, by = "PIDN") %>%
  dplyr::select(feature, PIDN, gene_name, AD_max_pp4, PD_max_pp4, best_coloc_tissue, Watershed_posterior_outlier_signal_1, Cohort, Age, Dx, APOE)
expression_merged_by_rare_variant <- expression_merged_by_patient %>%
  left_join(rare_variant_info_expression, by = c("PIDN" = "SubjectID", "feature" = "gene_id")) %>%
  dplyr::select(feature, PIDN, gene_name, AD_max_pp4, PD_max_pp4, best_coloc_tissue, Watershed_posterior_outlier_signal_1, Cohort, Age, Dx, APOE, sniffle)
expression_merged_by_rare_variant_over_under_outliers_annotated <- expression_merged_by_rare_variant %>%
  left_join(expression_with_zscore, by = c("PIDN" = "PIDN", "feature" = "gene_id")) %>%
  dplyr::select(feature, PIDN, gene_name.x, AD_max_pp4, PD_max_pp4, best_coloc_tissue, Watershed_posterior_outlier_signal_1, Age, Dx, APOE, sniffle, zscore, has_expression_outlier)

write.csv(expression_merged_by_rare_variant_over_under_outliers_annotated, "merged_by_rare_variant_over_under_expression_outliers_coloc_0.5_watershed_0.6.csv")


### plasma protein
reference_prot <- fread("~/tannerj/ADRComics/plasma_proteomics/maggie_dan_normalized/Proteome_Annotations.csv")
filtered_reference_prot <- reference_prot %>%
  group_by(Symbol, Definition) %>%
  filter(n() == 1) %>%
  ungroup()

reference_file <- fread("~/syang/ADRC/pre_processing_scripts/MANE.GRCh38.v1.3.summary.txt", sep = "\t")
reference_file$Ensembl_Gene <- str_replace(reference_file$Ensembl_Gene, "\\..*", "")
reference_with_gene <- filtered_reference_prot %>%
  left_join(reference_file, by = c("Symbol" = "symbol")) %>%
  dplyr::select(Protein, Symbol, Definition, Ensembl_Gene)

coloc_prot <- fread("~/tannerj/ADRComics/QTL_mapping/ADRC_coloc/brain_gwas/output/coloc/2025-01-31_21-55-27.234243_deCODE_UKBB_PQTLs/protein_pqtls_colocalizations.AD_PD_gene_list.txt")
coloc_prot$feature <- sub("\\..*", "", coloc_prot$feature)
coloc_0.5_conf_prot <- coloc_prot[which(coloc_prot$AD_max_pp4 >= 0.5 | coloc_prot$PD_max_pp4 >= 0.5), ]


coloc_0.5_conf_prot$matched_symbol <- NA
coloc_0.5_conf_prot$matched_ensembl <- NA

# Loop through each Symbol in reference_with_gene
for (i in seq_len(nrow(reference_with_gene))) {
  symbol <- reference_with_gene$Symbol[i]
  ensembl_gene <- reference_with_gene$Ensembl_Gene[i]

  matched_indices <- grepl(symbol, coloc_0.5_conf_prot$feature, fixed = TRUE)

  coloc_0.5_conf_prot$matched_symbol[matched_indices] <- ifelse(
    is.na(coloc_0.5_conf_prot$matched_symbol[matched_indices]),
    symbol,
    paste(coloc_0.5_conf_prot$matched_symbol[matched_indices], symbol, sep = "; ")
  )

  coloc_0.5_conf_prot$matched_ensembl[matched_indices] <- ifelse(
    is.na(coloc_0.5_conf_prot$matched_ensembl[matched_indices]),
    ensembl_gene,
    paste(coloc_0.5_conf_prot$matched_ensembl[matched_indices], ensembl_gene, sep = "; ")
  )
}

# Filter rows where a match was found
convert_coloc_0.5_conf_prot <- coloc_0.5_conf_prot[!is.na(coloc_0.5_conf_prot$matched_symbol), ]
convert_coloc_0.5_conf_prot %>% filter(!is.na(matched_ensembl))
convert_coloc_0.5_conf_prot <- convert_coloc_0.5_conf_prot %>%
  separate_rows(matched_symbol, matched_ensembl, sep = ";") %>%
  as.data.frame()


### CSF protein
csf_prot_conversion <- fread("~/syang/csf_soma7k_pQTL_aptamer_info.csv")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"), filters = "entrezgene_id", values = csf_prot_conversion$EntrezGeneID, mart = ensembl)
conversion$entrezgene_id <- as.character(conversion$entrezgene_id)
csf_prot_conversion <- merge(csf_prot_conversion, conversion, by.x = "EntrezGeneID", by.y = "entrezgene_id", all.x = TRUE)
csf_prot_conversion <- csf_prot_conversion[!is.na(csf_prot_conversion$ensembl_gene_id), ]

csf_coloc_AD <- fread("~/tannerj/ADRComics/QTL_mapping/ADRC_coloc/brain_gwas/output/coloc/2025-02-15_12-57-43.610250_CSF_pqtl_neurogwas_coloc/Bellenguez_AD_GWAS_sumamry_stats_formatted_tsv_gz_coloc_status.txt")
csf_coloc_AD <- csf_coloc_AD[csf_coloc_AD$eqtl_file == "CSF_Cruchaga_cis_QTLs_txt_gz", ]
csf_coloc_AD$gene_id <- csf_prot_conversion$ensembl_gene_id[match(csf_coloc_AD$feature, csf_prot_conversion$Analytes)]
csf_coloc_AD$prot_symbol <- csf_prot_conversion$EntrezGeneSymbol[match(csf_coloc_AD$feature, csf_prot_conversion$Analytes)]
csf_coloc_AD$feature <- ifelse(is.na(csf_coloc_AD$gene_id), gsub("0+$", "", csf_coloc_AD$feature), csf_coloc_AD$feature)
csf_coloc_AD$gene_id <- csf_prot_conversion$ensembl_gene_id[match(csf_coloc_AD$feature, csf_prot_conversion$Analytes)]
csf_coloc_AD$coloc_type <- "AD"

csf_coloc_PD1 <- fread("~/tannerj/ADRComics/QTL_mapping/ADRC_coloc/brain_gwas/output/coloc/2025-02-15_12-57-43.610250_CSF_pqtl_neurogwas_coloc/Kim_JJ_PD_GWAS_summary_stats_formatted_tsv_gz_coloc_status.txt")
csf_coloc_PD2 <- fread("~/tannerj/ADRComics/QTL_mapping/ADRC_coloc/brain_gwas/output/coloc/2025-02-15_12-57-43.610250_CSF_pqtl_neurogwas_coloc/Nalls_PD_GWAS_summary_stats_formatted_tsv_gz_coloc_status.txt")
csf_coloc_PD_set <- rbind(csf_coloc_PD1, csf_coloc_PD2)
colnames(csf_coloc_PD_set) <- colnames(csf_coloc_AD)[1:11]
csf_coloc_PD_set <- csf_coloc_PD_set[csf_coloc_PD_set$eqtl_file == "CSF_Cruchaga_cis_QTLs_txt_gz", ]
csf_coloc_PD_set$gene_id <- csf_prot_conversion$ensembl_gene_id[match(csf_coloc_PD_set$feature, csf_prot_conversion$Analytes)]
csf_coloc_PD_set$prot_symbol <- csf_prot_conversion$EntrezGeneSymbol[match(csf_coloc_PD_set$feature, csf_prot_conversion$Analytes)]
csf_coloc_PD_set$feature <- ifelse(is.na(csf_coloc_PD_set$gene_id), gsub("0+$", "", csf_coloc_PD_set$feature), csf_coloc_PD_set$feature)
csf_coloc_PD_set$gene_id <- csf_prot_conversion$ensembl_gene_id[match(csf_coloc_PD_set$feature, csf_prot_conversion$Analytes)]
csf_coloc_PD_set <- csf_coloc_PD_set[!is.na(csf_coloc_PD_set$gene_id), ]
csf_coloc_PD_set$coloc_type <- "PD"

csf_coloc_AD_PD <- rbind(csf_coloc_AD, csf_coloc_PD_set)

csf_prot_prioritized_variants <- fread("~/syang/Watershed-SV/run_watershed_combinations/csf_pro_25/model_RIVER_number_of_dimensions_1_on_csf_prot_pc_25_posterior_probability.txt")
csf_prot_prioritized_variants$ID <- sub(".*:", "", csf_prot_prioritized_variants$sample_names)
csf_prot_prioritized_variants$PIDN <- as.integer(gsub(":.*", "", csf_prot_prioritized_variants$sample_names))
csf_prot_prioritized_0.6_conf <- csf_prot_prioritized_variants[which(csf_prot_prioritized_variants$Watershed_posterior_outlier_signal_1 >= 0.6), ]

csf_prot_with_zscore <- fread("~/syang/ADRC/copy_unprocessed_data/output_filtered_csf_prot_pc_25.csv")

csf_coloc_watershed_prioritized <- merge(csf_coloc_AD_PD, csf_prot_prioritized_0.6_conf, by.x = "gene_id", by.y = "ID")
csf_prot_merged_by_patient <- csf_coloc_watershed_prioritized %>%
  left_join(sample_swap, by = "PIDN") %>%
  dplyr::select(gene_id, PIDN, eqtl_file, clpp_h4, coloc_type, Watershed_posterior_outlier_signal_1, Age, Dx, APOE, prot_symbol) %>%
  filter(clpp_h4 > 0.5)

csf_prot_with_zscore$PIDN <- as.character(csf_prot_with_zscore$PIDN)
csf_prot_merged_by_rare_variant <- csf_prot_merged_by_patient %>% left_join(rare_variant_info_csf_prot, by = c("PIDN" = "SubjectID", "gene_id" = "gene_id"), relationship = "many-to-many")
csf_prot_merged_by_rare_variant <- csf_prot_merged_by_rare_variant %>% dplyr::select(gene_id, PIDN, prot_symbol, clpp_h4, coloc_type, Watershed_posterior_outlier_signal_1, Age, Dx, APOE, sniffle)
csf_prot_merged_by_rare_variant$PIDN <- as.character(csf_prot_merged_by_rare_variant$PIDN)
csf_prot_merged_by_rare_variant_over_under_outliers_annotated <- csf_prot_merged_by_rare_variant %>%
  left_join(csf_prot_with_zscore, by = c("PIDN" = "PIDN", "gene_id" = "gene_id"), relationship = "many-to-many") %>%
  dplyr::select(gene_id, PIDN, prot_symbol, clpp_h4, coloc_type, Watershed_posterior_outlier_signal_1, Age, Dx, APOE, sniffle, zscore, outlier_status)
write.csv(csf_prot_merged_by_rare_variant_over_under_outliers_annotated, "csf_prot_merged_by_rare_variant_over_under_outliers_coloc_0.5_watershed_0.6.csv")

### plasma protein data
protein_prioritized_variants <- fread("~/syang/Watershed-SV/Watershed_strobe/model_RIVER_number_of_dimensions_prot_pc_25_1_posterior_probability.txt")
protein_with_zscore <- fread("~/syang/ADRC/copy_unprocessed_data/output_filtered_prot_pc_25.csv")

protein_prioritized_variants$ID <- sub(".*:", "", protein_prioritized_variants$sample_names)
protein_prioritized_variants$PIDN <- as.integer(gsub(":.*", "", protein_prioritized_variants$sample_names))
protein_prioritized_0.6_conf <- protein_prioritized_variants[which(protein_prioritized_variants$Watershed_posterior_outlier_signal_1 >= 0.6), ]


coloc_prot_watershed_prioritized <- merge(convert_coloc_0.5_conf_prot, protein_prioritized_0.6_conf, by.x = "matched_ensembl", by.y = "ID")
prot_merged_by_patient <- coloc_prot_watershed_prioritized %>%
  left_join(sample_swap, by = "PIDN") %>%
  dplyr::select(feature, PIDN, AD_max_pp4, PD_max_pp4, Watershed_posterior_outlier_signal_1, matched_ensembl, matched_symbol, Cohort, Age, Dx, APOE)
prot_merged_by_rare_variant <- prot_merged_by_patient %>%
  left_join(rare_variant_info_prot, by = c("PIDN" = "SubjectID", "matched_ensembl" = "gene_id"), relationship = "many-to-many") %>%
  dplyr::select(feature, PIDN, AD_max_pp4, PD_max_pp4, Watershed_posterior_outlier_signal_1, matched_ensembl, matched_symbol, Cohort, Age, Dx, APOE, sniffle)
prot_merged_by_rare_variant$PIDN <- as.character(prot_merged_by_rare_variant$PIDN)
prot_merged_by_rare_variant_over_under_outliers_annotated <- prot_merged_by_rare_variant %>%
  left_join(protein_with_zscore, by = c("PIDN" = "PIDN", "matched_ensembl" = "gene_id")) %>%
  dplyr::select(feature, PIDN, AD_max_pp4, PD_max_pp4, Watershed_posterior_outlier_signal_1, matched_symbol, Age, Dx, APOE, sniffle, zscore, has_prot_outlier)
write.csv(prot_merged_by_rare_variant_over_under_outliers_annotated, "prot_merged_by_rare_variant_over_under_expression_outliers_coloc_0.5_watershed_0.6.csv")

rare_sv_nearby_expressed_gene_counts <- rare_variant_info_expression %>% distinct(SubjectID, gene_id)
rare_sv_nearby_expressed_gene <- rare_sv_nearby_expressed_gene_counts %>% count(SubjectID)
rare_sv_nearby_expressed_gene_by_patient <- rare_sv_nearby_expressed_gene %>% nrow()

rare_sv_nearby_colocalized_gene_counts <- inner_join(coloc_expression_0.5_conf, rare_variant_info_expression, by = c("feature" = "gene_id")) %>% distinct(SubjectID, feature)
rare_sv_nearby_colocalized_gene <- rare_sv_nearby_colocalized_gene_counts %>% count(SubjectID)
rare_sv_nearby_colocalized_gene_by_patient <- rare_sv_nearby_colocalized_gene %>% nrow()

rare_sv_prioritized_by_watershed_counts <- inner_join(rare_variant_info_expression, expression_prioritized_0.6_conf, by = c("gene_id" = "ID", "SubjectID" = "PIDN")) %>% distinct(SubjectID, gene_id)
rare_sv_prioritized_by_watershed_number <- rare_sv_prioritized_by_watershed_counts %>% count(SubjectID)
rare_sv_prioritized_by_watershed_number_by_patient <- rare_sv_prioritized_by_watershed_number %>% nrow()






expression_outliers_z_2 <- fread("~/syang/Watershed-SV/run_watershed_combinations/rna_105/rna_formatted_for_watershed_pc_105.tsv")
expression_outliers_z_2 <- expression_outliers_z_2 %>% filter(abs(zscore) > 2)
rare_sv_nearby_outlier_z_2_counts <- inner_join(rare_variant_info_expression, expression_outliers_z_2, by = c("gene_id" = "gene", "SubjectID" = "ID")) %>% distinct(SubjectID, gene_id)
rare_sv_nearby_outlier_z_2_number <- rare_sv_nearby_outlier_z_2_counts %>% count(SubjectID)
rare_sv_nearby_outlier_z_2_by_patient <- rare_sv_nearby_outlier_z_2_number %>% nrow()

expression_near_colocalized_gene_and_prioritized_by_watershed_counts <- inner_join(rare_sv_nearby_colocalized_gene_counts, expression_prioritized_0.6_conf, by = c("SubjectID" = "PIDN", "feature" = "ID")) %>% distinct(SubjectID, feature)
expression_near_colocalized_gene_and_prioritized_by_watershed_number <- expression_near_colocalized_gene_and_prioritized_by_watershed_counts %>% count(SubjectID)
expression_near_colocalized_gene_and_prioritized_by_watershed_number_by_patient <- expression_near_colocalized_gene_and_prioritized_by_watershed_number %>% nrow()


convert_coloc_0.5_conf_prot <- convert_coloc_0.5_conf_prot[convert_coloc_0.5_conf_prot$AD_max_pp4 > 0.5 | convert_coloc_0.5_conf_prot$PD_max_pp4 > 0.5, ]
rare_variant_info_prot$SubjectID <- as.character(rare_variant_info_prot$SubjectID)
protein_prioritized_0.6_conf$PIDN <- as.character(protein_prioritized_0.6_conf$PIDN)


rare_sv_nearby_prot_gene_counts <- rare_variant_info_prot %>% distinct(SubjectID, gene_id)
rare_sv_nearby_prot_gene <- setDT(rare_sv_nearby_prot_gene_counts %>% count(SubjectID))
rare_sv_nearby_prot_gene_by_patient <- rare_sv_nearby_prot_gene %>% nrow()

rare_sv_nearby_colocalized_prot_counts <- inner_join(convert_coloc_0.5_conf_prot, rare_variant_info_prot, by = c("matched_ensembl" = "gene_id")) %>% distinct(SubjectID, matched_ensembl)
rare_sv_nearby_colocalized_prot <- setDT(rare_sv_nearby_colocalized_prot_counts %>% count(SubjectID))
rare_sv_nearby_colocalized_gene_prot_by_patient <- rare_sv_nearby_colocalized_prot %>% nrow()

rare_sv_prioritized_by_watershed_prot_counts <- inner_join(rare_variant_info_prot, protein_prioritized_0.6_conf, by = c("gene_id" = "ID", "SubjectID" = "PIDN")) %>% distinct(SubjectID, gene_id)
rare_sv_prioritized_by_watershed_prot_number <- rare_sv_prioritized_by_watershed_prot_counts %>% count(SubjectID)
rare_sv_prioritized_by_watershed_prot_by_patient <- rare_sv_prioritized_by_watershed_prot_number %>% nrow()






protein_outliers_z_2 <- fread("~/syang/Watershed-SV/formatted_data/prot_formatted_for_watershed_pc_25.tsv")
rare_sv_nearby_outlier_z_2_prot_counts <- inner_join(rare_variant_info_prot, protein_outliers_z_2, by = c("gene_id" = "gene", "SubjectID" = "ID")) %>%
  filter(abs(zscore) > 2) %>%
  distinct(SubjectID, gene_id)
rare_sv_nearby_outlier_z_2_number_prot <- rare_sv_nearby_outlier_z_2_prot_counts %>% count(SubjectID)
rare_sv_nearby_outlier_z_2_prot_by_patient <- rare_sv_nearby_outlier_z_2_number_prot %>% nrow()

prot_near_colocalized_gene_and_prioritized_by_watershed_counts <- inner_join(rare_sv_nearby_colocalized_prot_counts, protein_prioritized_0.6_conf, by = c("SubjectID" = "PIDN", "matched_ensembl" = "ID")) %>% distinct(SubjectID, matched_ensembl)
prot_near_colocalized_gene_and_prioritized_by_watershed_number <- setDT(prot_near_colocalized_gene_and_prioritized_by_watershed_counts %>% count(SubjectID))
prot_near_colocalized_gene_and_prioritized_by_watershed_by_patient <- prot_near_colocalized_gene_and_prioritized_by_watershed_number %>% nrow()



csf_coloc_conf_0.5_AD_PD <- csf_coloc_AD_PD %>% filter(clpp_h4 > 0.5)
rare_variant_info_csf_prot$SubjectID <- as.character(rare_variant_info_csf_prot$SubjectID)
csf_prot_prioritized_0.6_conf$PIDN <- as.character(csf_prot_prioritized_0.6_conf$PIDN)

rare_sv_nearby_csf_prot_gene_counts <- rare_variant_info_csf_prot %>% distinct(SubjectID, gene_id)
rare_sv_nearby_csf_prot_gene <- rare_sv_nearby_csf_prot_gene_counts %>% count(SubjectID)
rare_sv_nearby_csf_prot_gene_by_patient <- rare_sv_nearby_csf_prot_gene %>% nrow()

rare_sv_nearby_colocalized_csf_prot_counts <- inner_join(csf_coloc_conf_0.5_AD_PD, rare_variant_info_csf_prot, by = c("gene_id" = "gene_id")) %>% distinct(SubjectID, gene_id)
rare_sv_nearby_colocalized_csf_prot <- rare_sv_nearby_colocalized_csf_prot_counts %>% count(SubjectID)
rare_sv_nearby_colocalized_csf_prot_by_patient <- rare_sv_nearby_colocalized_csf_prot %>% nrow()


rare_sv_prioritized_by_watershed_csf_prot_counts <- inner_join(rare_variant_info_csf_prot, csf_prot_prioritized_0.6_conf, by = c("gene_id" = "ID", "SubjectID" = "PIDN")) %>% distinct(SubjectID, gene_id)
rare_sv_prioritized_by_watershed_csf_prot_number <- rare_sv_prioritized_by_watershed_csf_prot_counts %>% count(SubjectID)
rare_sv_prioritized_by_watershed_csf_prot_by_patient <- rare_sv_prioritized_by_watershed_csf_prot_number %>% nrow()







csf_protein_outliers_z_2 <- fread("~/syang/Watershed-SV/formatted_data/csf_prot_formatted_for_watershed_pc_25.tsv")
csf_protein_outliers_z_2$ID <- as.character(csf_protein_outliers_z_2$ID)
rare_sv_nearby_outlier_z_2_csf_prot_counts <- inner_join(rare_variant_info_csf_prot, csf_protein_outliers_z_2, by = c("gene_id" = "gene", "SubjectID" = "ID")) %>%
  filter(abs(zscore) > 2) %>%
  distinct(SubjectID, gene_id)
rare_sv_nearby_outlier_z_2_csf_prot <- rare_sv_nearby_outlier_z_2_csf_prot_counts %>% count(SubjectID)
rare_sv_nearby_outlier_z_2_csf_prot_by_patient <- rare_sv_nearby_outlier_z_2_csf_prot %>% nrow()


csf_prot_near_colocalized_gene_and_prioritized_by_watershed_counts <- inner_join(rare_sv_nearby_colocalized_csf_prot_counts, csf_prot_prioritized_0.6_conf, by = c("SubjectID" = "PIDN", "gene_id" = "ID")) %>% distinct(SubjectID, gene_id)
csf_prot_near_colocalized_gene_and_prioritized_by_watershed <- csf_prot_near_colocalized_gene_and_prioritized_by_watershed_counts %>% count(SubjectID)
csf_prot_near_colocalized_gene_and_prioritized_by_watershed_by_patient <- csf_prot_near_colocalized_gene_and_prioritized_by_watershed %>% nrow()



### assemble dataframe for plotting
x_values <- c(
  rare_sv_nearby_expressed_gene[, 2],
  rare_sv_nearby_colocalized_gene[, 2],
  rare_sv_nearby_outlier_z_2_number[, 2],
  rare_sv_prioritized_by_watershed_number[, 2],
  expression_near_colocalized_gene_and_prioritized_by_watershed_number[, 2]
)

x_values_patients <- c(
  rare_sv_nearby_expressed_gene_by_patient,
  rare_sv_nearby_colocalized_gene_by_patient,
  rare_sv_nearby_outlier_z_2_by_patient,
  rare_sv_prioritized_by_watershed_number_by_patient,
  expression_near_colocalized_gene_and_prioritized_by_watershed_number_by_patient
)

x_values_prot <- c(
  rare_sv_nearby_prot_gene[, 2],
  rare_sv_nearby_colocalized_prot[, 2],
  rare_sv_nearby_outlier_z_2_number_prot[, 2],
  rare_sv_prioritized_by_watershed_prot_number[, 2],
  prot_near_colocalized_gene_and_prioritized_by_watershed_number[, 2]
)

x_values_patients_prot <- c(
  rare_sv_nearby_prot_gene_by_patient,
  rare_sv_nearby_colocalized_gene_prot_by_patient,
  rare_sv_nearby_outlier_z_2_prot_by_patient,
  rare_sv_prioritized_by_watershed_prot_by_patient,
  prot_near_colocalized_gene_and_prioritized_by_watershed_by_patient
)

x_values_csf_prot <- c(
  rare_sv_nearby_csf_prot_gene[, 2],
  rare_sv_nearby_colocalized_csf_prot[, 2],
  rare_sv_nearby_outlier_z_2_csf_prot[, 2],
  rare_sv_prioritized_by_watershed_csf_prot_number[, 2],
  csf_prot_near_colocalized_gene_and_prioritized_by_watershed[, 2]
)

x_values_patients_csf_prot <- c(
  rare_sv_nearby_csf_prot_gene_by_patient,
  rare_sv_nearby_colocalized_csf_prot_by_patient,
  rare_sv_nearby_outlier_z_2_csf_prot_by_patient,
  rare_sv_prioritized_by_watershed_csf_prot_by_patient,
  csf_prot_near_colocalized_gene_and_prioritized_by_watershed_by_patient
)





### plotting figure 3e
set.seed(42)
categories <- c("expression", "plasma_protein", "csf_protein")
filter_levels <- c("rareSVs nearby expressed gene", "rareSVs nearby colocalized gene, p>0.5", "rareSVs nearby outlier z>2", "rareSVs prioritized by Watershed p>0.6", "rareSVs nearby colocalized gene and prioritized by Watershed")

df_expression <- data.frame(
  `Disease variant filters` = rep(filter_levels, times = sapply(x_values, length)),
  `Number SV-Gene pairs prioritized` = unlist(x_values),
  `Filter Category` = "expression"
)

df_protein <- data.frame(
  `Disease variant filters` = rep(filter_levels, times = sapply(x_values_prot, length)),
  `Number SV-Gene pairs prioritized` = unlist(x_values_prot),
  `Filter Category` = "plasma_protein"
)

df_csf <- data.frame(
  `Disease variant filters` = rep(filter_levels, times = sapply(x_values_csf_prot, length)),
  `Number SV-Gene pairs prioritized` = unlist(x_values_csf_prot),
  `Filter Category` = "csf_protein"
)

df_csf$Disease.variant.filters <- factor(df_csf$Disease.variant.filters, levels = filter_levels)
df <- bind_rows(df_expression, df_protein)
df$Disease.variant.filters <- factor(df$Disease.variant.filters, levels = filter_levels)

warm_terracotta <- "#E07A5F"

# Define color palette
category_colors <- c(
  "expression" = "#eb346f",
  "plasma_protein" = "#7845e3",
  "csf_protein" = "#4a89f7"
)
common_theme <- theme(
  axis.text.x = element_text(size = 14, color = "black"),
  axis.title.x = element_text(size = 14, color = "black"),
  axis.line = element_line(size = 1.5)
)


# Make a function to build a data frame from values
make_df <- function(values, category_label) {
  data.frame(
    `Disease variant filters` = factor(rep(filter_levels, times = sapply(values, length)), levels = filter_levels),
    `Number SV-Gene pairs prioritized` = unlist(values),
    `Filter Category` = category_label
  )
}

# Assemble data frames
df_expression <- make_df(x_values, "expression")
df_protein <- make_df(x_values_prot, "plasma_protein")
df_csf_protein <- make_df(x_values_csf_prot, "csf_protein")

# Combine all into one
df <- bind_rows(df_expression, df_protein, df_csf_protein)
df$Filter.Category <- factor(df$Filter.Category, levels = c("expression", "plasma_protein", "csf_protein"))

# Create a composite y-axis label
df$Filter_Row <- paste(df$Filter.Category, df$Disease.variant.filters, sep = " - ")

# Ensure the y-axis shows Filter_Row in correct order (expression on top)
df <- df %>%
  arrange(factor(df$Filter.Category, levels = c("expression", "plasma_protein", "csf_protein")), Disease.variant.filters)

# Then set the levels accordingly (reversed to plot top-down)
df$Filter_Row <- factor(df$Filter_Row, levels = rev(unique(df$Filter_Row)))


pstrip <- ggplot(df, aes(
  x = Number.SV.Gene.pairs.prioritized,
  y = Filter_Row,
  color = Filter.Category
)) +
  geom_jitter(height = 0.2, width = 0, size = 3, alpha = 0.6) +
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    color = "black",
    alpha = 0.75,
  ) +
  scale_x_log10() +
  scale_color_manual(values = category_colors) +
  theme_minimal() +
  labs(
    x = "Number of SV-Gene pairs prioritized",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, color = "black"),
    axis.line = element_line(size = 1)
  )

print(pstrip)
ggsave("~/tannerj/pstrip.pdf")


# Patient count bar plot
df_expression_patients <- data.frame(
  `Disease variant filters` = c("rareSVs nearby expressed gene", "rareSVs nearby colocalized gene, p>0.5", "rareSVs nearby outlier z>2", "rareSVs prioritized by Watershed p>0.6", "rareSVs nearby colocalized gene and prioritized by Watershed"),
  `Patient Count` = x_values_patients,
  `Filter Category` = "expression"
)

df_expression_patients$Disease.variant.filters <- factor(df_expression_patients$Disease.variant.filters, levels = filter_levels)

# Data frame for "protein" patient counts
df_protein_patients <- data.frame(
  `Disease variant filters` = c("rareSVs nearby expressed gene", "rareSVs nearby colocalized gene, p>0.5", "rareSVs nearby outlier z>2", "rareSVs prioritized by Watershed p>0.6", "rareSVs nearby colocalized gene and prioritized by Watershed"),
  `Patient Count` = x_values_patients_prot,
  `Filter Category` = "plasma_protein"
)

df_protein_patients$Disease.variant.filters <- factor(df_protein_patients$Disease.variant.filters, levels = filter_levels)

df_csf_protein_patients <- data.frame(
  `Disease variant filters` = c("rareSVs nearby expressed gene", "rareSVs nearby colocalized gene, p>0.5", "rareSVs nearby outlier z>2", "rareSVs prioritized by Watershed p>0.6", "rareSVs nearby colocalized gene and prioritized by Watershed"),
  `Patient Count` = x_values_patients_csf_prot,
  `Filter Category` = "csf_protein"
)

df_csf_protein_patients$Disease.variant.filters <- factor(df_csf_protein_patients$Disease.variant.filters, levels = filter_levels)

# Combine patient data
df_patients <- bind_rows(df_expression_patients, df_protein_patients, df_csf_protein_patients)

# Extract the order of Filter_Row from the main df (used in pstrip)
desired_order <- unique(df$Filter_Row)

# Create matching Filter_Row column in df_patients to align with df
df_patients$Filter_Row <- paste(df_patients$Filter.Category, df_patients$Disease.variant.filters, sep = " - ")

# Set Filter_Row as factor with the desired order
df_patients$Filter_Row <- factor(df_patients$Filter_Row, levels = desired_order)

# Now set Disease.variant.filters based on that ordered Filter_Row
df_patients <- df_patients %>%
  arrange(Filter_Row) %>%
  mutate(Disease.variant.filters = factor(Disease.variant.filters, levels = unique(Disease.variant.filters)))
df_patients$Disease.variant.filters


library(magrittr)
df_patients$Disease.variant.filters %<>% factor(rev(levels(df$Disease.variant.filters)))
df_patients$Filter.Category <- factor(df_patients$Filter.Category, levels = c("expression", "plasma_protein", "csf_protein"))
p2 <- ggplot(df_patients, aes(x = Patient.Count, y = Disease.variant.filters, fill = Filter.Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  labs(x = "Patient count", y = NULL) +
  theme(legend.position = "none", axis.title.y = element_blank(), strip.text = element_blank(), axis.text.y = element_blank()) +
  facet_wrap(~Filter.Category, ncol = 1)

p2 <- p2 + common_theme

pstrip + p2 + plot_layout(width = c(3, 1))
combined_plots <- pstrip + p2 + plot_layout(width = c(3, 1))
plot_grid(plotlist = list(pstrip + theme(legend.position = "none"), p2), rel_widths = c(.75, .25), nrow = 1)
ggplot2::ggsave("~/tannerj/ADRComics/summary_plots_for_colocalation_and_rar_variants_by_counts_and_patients.pdf", width = 14, height = 9)



filters <- c(
  "rareSVs nearby colocalized gene and prioritized by Watershed",
  "rareSVs prioritized by Watershed p>0.6",
  "rareSVs nearby outlier z>2",
  "rareSVs nearby colocalized gene, p>0.5",
  "rareSVs nearby expressed gene"
)


# Create patient count data frames
df_expression_patients <- data.frame(
  Disease.variant.filters = filters,
  Patient.Count = x_values_patients,
  Filter.Category = "expression"
)

df_protein_patients <- data.frame(
  Disease.variant.filters = filters,
  Patient.Count = x_values_patients_prot,
  Filter.Category = "plasma_protein"
)


df_csf_protein_patients <- data.frame(
  Disease.variant.filters = filters,
  Patient.Count = x_values_patients_csf_prot,
  Filter.Category = "csf_protein"
)

# Combine all patient data
df_patients <- bind_rows(df_expression_patients, df_protein_patients, df_csf_protein_patients)

df_patients$Disease.variant.filters <- factor(df_patients$Disease.variant.filters, levels = rev(filters)) # rev for top-to-bottom order
df_patients$Filter.Category <- factor(df_patients$Filter.Category, levels = c("expression", "plasma_protein", "csf_protein"))

# Create patient count bar plot
p2 <- ggplot(df_patients, aes(x = Patient.Count, y = Disease.variant.filters, fill = Filter.Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = category_colors) +
  facet_wrap(~Filter.Category, ncol = 1) +
  labs(x = "Patient count", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    strip.text = element_blank(),
    axis.text.y = element_blank()
  ) +
  common_theme

print(p2)
# Combine plots
combined_plots_bar <- pstrip + p2 + plot_layout(widths = c(3, 1))

# Display and save
print(combined_plots_bar)
ggsave("summary_plots_for_colocalation_and_rar_variants_by_counts_and_patients.png",
  plot = combined_plots_bar, width = 14, height = 9, dpi = 300
)


### plotting figure 3d
set.seed(123)

total_colocs <- c(csf_coloc_genes, expression_colocs) %>% unique()
csf_coloc_genes <- csf_coloc_AD_PD %>%
  filter(clpp_h4 > .5) %>%
  pull(gene_id) %>%
  unique()
prot_coloc_genes <- coloc_prot %>%
  filter(clpp_h4 > .5) %>%
  pull(gene_id) %>%
  unique()
csf_prot_prioritized_variants <- csf_prot_prioritized_variants %>% rename(gene_id = ID)

csf_prot_prioritized_variants$coloc <- ifelse(csf_prot_prioritized_variants$gene_id %in% total_colocs, yes = "coloc", no = "background")
csf_prot_prioritized_variants %>% filter(coloc == "coloc", Watershed_posterior_outlier_signal_1 > .6)
csf.plot <- ggplot(csf_prot_prioritized_variants, aes(Watershed_posterior_outlier_signal_1, fill = coloc)) +
  geom_histogram(bins = 200) +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = c(1, 10, 1000, 10000)) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line(linewidth = .25)) +
  scale_fill_manual(values = c("#808184", "#8C5AA3")) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red") +
  xlab("")

prot_prioritized_variants <- protein_prioritized_variants %>% rename(gene_id = ID)
prot_prioritized_variants$coloc <- ifelse(prot_prioritized_variants$gene_id %in% total_colocs, yes = "coloc", no = "background")
prot.plot <- ggplot(prot_prioritized_variants, aes(Watershed_posterior_outlier_signal_1, fill = coloc)) +
  geom_histogram(bins = 200) +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = c(1, 10, 1000, 10000)) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line(linewidth = .25)) +
  scale_fill_manual(values = c("#808184", "#8C5AA3")) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red") +
  xlab("")

expr.plot <- ggplot(expression_prioritized_variants, aes(Watershed_posterior_outlier_signal_1, fill = coloc)) +
  geom_histogram(bins = 200) +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = c(1, 10, 1000, 10000)) +
  theme(legend.position = "none") +
  xlab("") +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line(linewidth = .25)) +
  scale_fill_manual(values = c("#808184", "#8C5AA3")) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red")
expr.plot

watershed_posteriors <- rbind(
  expression_prioritized_variants %>% mutate(ome = "RNA"),
  prot_prioritized_variants %>% mutate(ome = "Plasma Protein"),
  csf_prot_prioritized_variants %>% mutate(ome = "CSF Protein")
)
watershed_posteriors$ome %<>% factor(levels = c("RNA", "Plasma Protein", "CSF Protein"))
ggplot(watershed_posteriors, aes(Watershed_posterior_outlier_signal_1, fill = coloc)) +
  geom_histogram(bins = 200) +
  theme_classic() +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = c(1, 10, 1000, 10000)) +
  theme(legend.position = "none") +
  xlab("") +
  facet_wrap(~ome, ncol = 1) +
  theme(panel.grid.major = element_line(), panel.grid.minor = element_line(linewidth = .25), panel.border = element_rect(fill = NA)) +
  scale_fill_manual(values = c("#808184", "#8C5AA3")) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red")
ggsave("~/tannerj/ADRComics/watershed.distribution.prioritization.pdf", height = 6, width = 5)
expression_prioritized_variants <- expression_prioritized_variants %>% rename(gene_id = ID)
expression_prioritized_variants$coloc <- ifelse(expression_prioritized_variants$gene_id %in% total_colocs, yes = "coloc", no = "background")
