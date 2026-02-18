####################
# title:"rare structural variants and multi-ome reads processing"
# author:"Sherry Yang (adapted from Tanner Jensen)"
# date:"2026-02-16"
####################

# This script aims to process the relevant omes through pc corrections for largest rare variant enrichment.
# Steps include covariates remove/correction and rare variant enrichment for ranges of pc used

library(data.table)
library(corrplot)
library(tidyverse)
library(magrittr)
library(bestNormalize)
library(PCAtools)
library(limma)
library(Matrix)
library(rtracklayer)
library(data.table)
library(tidyverse)
library(plyranges)
library(magrittr)
library(dplyr)
library(purrr)
library(RNOmni)
library(data.table)
library(tidyverse)
library(magrittr)
library(PCAtools)
library(broom)
library(biomaRt)
library(limma)
library(dplyr)
library(ggplot2)

########## for rareSV processing and analysis

args <- commandArgs(trailingOnly = TRUE)
af_thresholds <- as.numeric(args[1])
flanking_windows <- as.numeric(args[2])
input_file <- args[3]
output_file <- as.character(args[4])


svs_file <- input_file
current_dir <- getwd()
parent_parent_dir <- dirname(dirname(dirname(current_dir)))
gff_file_path <- file.path(parent_parent_dir, "tannerj/ADRComics/ADRC_PBMC_scRNA/adrc/gencode.v32.annotation.gff3.gz")
gff <- readGFFAsGRanges(gff_file_path)
genes <- gff[gff$type == "gene", ]
svs <- fread(svs_file)
swapfile_path <- file.path(parent_parent_dir, "syang/ADRC/ADRC_SAMS.metadata.master_table.sample_swaps_annotated.tsv")
swapfile <- fread(swapfile_path)

# formatting and qc filter
svs <- svs[svs$gt > 0, ]
svs$length <- svs$end - svs$start
## any gene size over 10 million will simply be unrealistic, so excluding that for quality control
svs <- svs[abs(svs$length) < 10000000, ]

svs_gr <- makeGRangesFromDataFrame(svs, keep.extra.columns = TRUE)


for (flanking_window in flanking_windows) {
    for (af_threshold in af_thresholds) {
        start(genes) <- start(genes) - flanking_window
        end(genes) <- end(genes) + flanking_window

        genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)
        genes.subset <- genes %>% select(type, gene_id, gene_type, gene_name)

        svs_gr %<>% join_overlap_left(genes.subset)
        svs_gr_df <- as.data.frame(svs_gr) %>%
            select(variant_id, gene_id, gene_type, gene_name, ID, gt, AF, PIDN, length) %>%
            filter(gene_type == "protein_coding")
        svs_gr_df$gene_id <- sub("\\..*$", "", svs_gr_df$gene_id)
        svs_gr_df$PIDN <- as.character(svs_gr_df$PIDN)

        sv_genes <- na.omit(unique(svs_gr_df$gene_id))
        sv_samples <- na.omit(unique(svs_gr_df$ID))
        sv_sample_gene_mat <- matrix(nrow = length(sv_genes), ncol = length(sv_samples), 0)
        colnames(sv_sample_gene_mat) <- sv_samples
        rownames(sv_sample_gene_mat) <- sv_genes


        sv_outlier_gene_df <- data.frame(gene_id = sv_genes, sv_sample_gene_mat) %>%
            pivot_longer(cols = c(-gene_id), names_to = "ID", values_to = "has_sv_outlier") %>%
            left_join(as.data.frame(svs_gr_df))
        sv_outlier_gene_df <- sv_outlier_gene_df[which(!is.na(sv_outlier_gene_df$AF)), ]
        sv_outlier_gene_df_unique <- sv_outlier_gene_df %>%
            group_by(gene_id, ID) %>%
            select(gene_id, ID, has_sv_outlier, AF, length) %>%
            summarize(AF = map_dbl(AF, ~ .x[which.min(.x)]))
        sv_outlier_gene_df_unique$has_sv_outlier <- ifelse(sv_outlier_gene_df_unique$AF <= af_threshold, 1, 0)
        sv_outlier_gene_df_unique_pidn <- sv_outlier_gene_df_unique %>%
            left_join(swapfile %>% select(LRS_ID, PIDN, DROP_SAMPLE), by = c("ID" = "LRS_ID")) %>%
            filter(DROP_SAMPLE == FALSE)

        sv_outlier_gene_df_unique_pidn_subset <- sv_outlier_gene_df_unique_pidn %>% select(gene_id, PIDN, AF, has_sv_outlier)

        fwrite(sv_outlier_gene_df_unique_pidn_subset, output_file, row.names = FALSE)
    }
}






########## for plasma protein processing
raw_prot <- fread("~/tannerj/ADRComics/plasma_proteomics/maggie_dan_normalized/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv")
prot_metadat <- fread("~/tannerj/ADRComics/plasma_proteomics/maggie_dan_normalized/Chip_Phen_Clock_03_02_at_15_48.csv")
prot_metadat_select <- prot_metadat %>%
    select(SampleId, PIDN, ADRC_ID, Visit, Gender, Diagnosis_group, Age_at_Prot_Draw, Age_at_Draw_Difference, PlateId) %>%
    mutate(
        Age = Age_at_Prot_Draw + Age_at_Draw_Difference, Sex_M = as.numeric(Gender == "M"), dx_AD = as.numeric(Diagnosis_group %in% c("AD", "MCI", "MCI-AD")),
        dx_PD = as.numeric(Diagnosis_group %in% c("LBD", "MCI-PD", "PD", "PDD"))
    ) %>%
    distinct()
prot_metadat_select %<>% group_by(SampleId) %>% slice(1)

raw_prot <- raw_prot[raw_prot$SampleType == "Sample", ] # remove control / calibration samples
raw_prot <- raw_prot[raw_prot$SampleId %in% prot_metadat_select$SampleId, ]

table(raw_prot$SampleId)[which(table(raw_prot$SampleId) > 1)]

raw_prot.mat <- raw_prot[, 4:ncol(raw_prot)] %>%
    as.matrix() %>%
    t()
cormat <- cor(raw_prot.mat)
mean_sample_cors <- rowMeans(cormat)
ggplot(NULL, aes(x = mean_sample_cors)) +
    geom_histogram(bins = 100) +
    xlim(.75, 1) +
    geom_vline(xintercept = .81, color = "red", linetype = "dashed")
ggsave("portein_correlation_matrix.pdf")

exclude_samples <- which(mean_sample_cors < quantile(mean_sample_cors, .005))
cormat[cormat < .7] <- .7

raw_prot <- raw_prot[-exclude_samples, ]
raw_prot.mat_filt <- raw_prot.mat[, -exclude_samples]

# find variable proteins
library(Matrix)
row_means <- rowMeans(raw_prot.mat)
row_var <- apply(raw_prot.mat, 1, var)
CV <- apply(raw_prot.mat, 1, sd) / rowMeans(raw_prot.mat) * 100
table(CV > 50)
var_probes <- which(CV > 50)

raw_prot.mat_normalized <- apply(raw_prot.mat_filt, 1, function(x) boxcox(x)$x.t) %>% t()
pc.result <- pca(raw_prot.mat_normalized)
pc.result$metadata <- data.frame(SampleId = raw_prot$SampleId)
pc.result$metadata %<>% left_join(prot_metadat_select)
pc.result$metadata$decade <- round(pc.result$metadata$Age / 10) * 10

biplot(pc.result, colby = "PlateId", lab = NULL)
biplot(pc.result, colby = "Gender", lab = NULL)
biplot(pc.result, colby = "Diagnosis_group", lab = NULL, legendPosition = "bottom")
biplot(pc.result, colby = "decade", lab = NULL, legendPosition = "bottom")
ggsave("./normalized_protein_biplot.pdf")

pc_ranges <- seq(5, 200, 5)

for (number_pcs_correct in pc_ranges) {
    covariates <- data.frame(pc.result$metadata %>% select(Age, Sex_M, dx_AD, dx_PD), pc.result$rotated[, 1:number_pcs_correct])
    raw_prot.mat_normalized.corrected <- removeBatchEffect(raw_prot.mat_normalized, covariates = covariates)
    raw_prot.mat_normalized.corrected.scaled <- raw_prot.mat_normalized.corrected %>%
        t() %>%
        scale() %>%
        t()

    zscores <- raw_prot.mat_normalized.corrected.scaled %>% melt()
    colnames(zscores) <- c("probeID", "SampleId", "zscore")
    zscores$SampleId <- pc.result$metadata$SampleId[zscores$SampleId]

    sample_info <- prot_metadat_select %>% select(SampleId, PIDN, ADRC_ID)
    zscores %<>% left_join(sample_info) %>% select(-SampleId)

    zscore.individual_collapse <- zscores %>%
        group_by(PIDN, probeID, ADRC_ID) %>%
        summarize(zscore = median(zscore))
    probe_meta <- data.frame(probeID = unique(zscores$probeID)) %>% mutate(gene = gsub(probeID, pattern = "(\\w+)\\..*", replacement = "\\1"))
    zscore.individual_collapse %<>% left_join(probe_meta)


    output_filename <- paste0(
        "~/syang/ADRC/unprocessed_data/ADRC.plasma_protein.boxcox_normalized.corrected.",
        number_pcs_correct, "PCs.individual_collapsed_zscores.tsv"
    )

    fwrite(zscore.individual_collapse, file = output_filename, sep = "\t", row.names = F, col.names = T, quote = F)
}

number_outliers <- do.call(rbind, lapply(seq(2, 5, .5), function(n) {
    zscore.individual_collapse %>%
        group_by(PIDN) %>%
        summarize(num_outliers = sum(abs(zscore) > n)) %>%
        mutate(z_threshold = n)
}))
number_over <- do.call(rbind, lapply(seq(2, 5, .5), function(n) {
    zscore.individual_collapse %>%
        group_by(PIDN) %>%
        summarize(num_outliers = sum(zscore > n)) %>%
        mutate(z_threshold = n, direction = "over")
}))
number_under <- do.call(rbind, lapply(seq(2, 5, .5), function(n) {
    zscore.individual_collapse %>%
        group_by(PIDN) %>%
        summarize(num_outliers = sum(zscore < (-1 * n))) %>%
        mutate(z_threshold = n, direction = "under")
}))

number_direction <- rbind(number_over, number_under)
ggplot(number_outliers, aes(factor(z_threshold), num_outliers, fill = z_threshold)) +
    geom_boxplot() +
    scale_y_log10()
ggsave("number_outliers.pdf")

ggplot(number_direction, aes(factor(z_threshold), num_outliers, fill = direction, alpha = z_threshold)) +
    geom_boxplot() +
    scale_y_log10()
ggsave("number_outliers.pdf")


global_outliers <- number_outliers %>%
    filter(z_threshold == 3) %>%
    mutate(global_outlier_z = as.numeric(scale(num_outliers))) %>%
    filter(global_outlier_z > 3) %>%
    pull(PIDN)
pc.result$metadata$global_outlier <- pc.result$metadata$PIDN %in% global_outliers

biplot(pc.result, colby = "global_outlier", legendPosition = "bottom", lab = NULL)
ggsave("biplot.global_outliers.pdf")

write.table(zscore.individual_collapse, file = "./ADRC.plasma_protein.boxcox_normalized.corrected.individual_collapsed_zscores.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
global_outliers

zscore.individual_collapse %>% filter(gene == "HSPA1A", zscore < -2)


### CSF Protein processing
csf.meta.full <- fread("~/tannerj/ADRComics/somascan_csf_Proteomics/Metadata_CSF_demographics_May2021_v2_additionalQC_1.csv") %>% filter(SampleType == "Sample")
csf.meta <- csf.meta.full %>% filter(Study %in% c("ADRC", "Wagner"))
csf.meta$PIDN <- as.integer(csf.meta$PIDN)
csf.meta %<>% filter(!is.na(PIDN))

prot.meta <- fread("~/tannerj/ADRComics/somascan_csf_Proteomics/ProteinMetadata_with_LOD1.csv")
csf_samples <- unique(csf.meta$PIDN)

LR.mat <- fread("~/tannerj/ADRComics/QTL_mapping/genotypes/ADRC.cohort_combined.SV_genotype.matrix.common_MAF_0.05.txt")
LR_samples <- as.integer(colnames(LR.mat)[6:ncol(LR.mat)])
SR_LRS_matched <- fread("~/tannerj/ADRComics/QTL_mapping/susie_finemapping/methylation_finemapping/SRS_matched_snps.fam")$V1

table(csf_samples %in% LR_samples)
table(csf_samples %in% SR_LRS_matched)

raw_data <- fread("~/tannerj/ADRComics/somascan_csf_Proteomics/CSFProts.log10.csv")
raw_data[, 1:10]

csf.mat <- as.matrix(raw_data[raw_data$SampleType == "Sample", 3:ncol(raw_data)])
rownames(csf.mat) <- raw_data$SampleId[raw_data$SampleType == "Sample"]
head(csf.mat[, 1:10])

# remove duplicate samples, get baseline only
sample_ids_map <- csf.meta.full %>%
    filter(!PIDN %in% c("", "#N/A")) %>%
    group_by(PIDN) %>%
    summarize(best_id = SampleId[which.max(ConnectivityZscore)])

sample_ids_keep <- sample_ids_map$best_id
csf.meta <- csf.meta.full %>% filter(SampleId %in% sample_ids_keep)
csf.mat <- csf.mat[sample_ids_keep, ]

# rank norm for QTL mapping
vst_scaled <- apply(csf.mat, 1, RankNorm)
pca_scaled <- pca(vst_scaled)
var_scaled <- as.vector(vst_scaled) %>% var()
vst_npc <- chooseGavishDonoho(vst_scaled, var.explained = pca_scaled$sdev^2, noise = var_scaled)

# # unnormalized PC for
table(csf.meta$PlateId)
colnames(csf.meta)
table(csf.meta$Gender)
summary(csf.meta$Storage_time)
table(csf.meta$PlateId, csf.meta$SlideId, csf.meta$ScannerID)
covariates <- c("Age", "Gender", "PlateId", "Storage_time")


csf_pca <- pca(csf.mat %>% t())
var_explained <- data.frame(var_explained = (csf_pca$sdev^2) / sum(csf_pca$sdev^2), PC_number = 1:length(csf_pca$sdev))

pca_covariates <- csf_pca$rotated[, 1:20]
csf.meta %<>% as.data.frame
rownames(csf.meta) <- csf.meta$SampleId
csf.meta <- csf.meta[rownames(csf.mat), ]

csf.meta$Storage_time[is.na(csf.meta$Storage_time)] <- median(csf.meta$Storage_time, na.rm = T)
csf.meta$Gender[is.na(csf.meta$Gender)] <- "F"
csf.meta$Age[is.na(csf.meta$Age)] <- median(csf.meta$Age, na.rm = T)

covariates.mat <- as.matrix(cbind(pca_covariates, csf.meta$Gender == "M", csf.meta$Age, csf.meta$Storage_time))
covariates.mat
csf_residualized <- removeBatchEffect(csf.mat %>% t(), csf.meta$PlateId, covariates = covariates.mat)
csf_residualized.scaled <- csf_residualized %>%
    t() %>%
    scale() %>%
    t()
global_outliers <- colSums(abs(csf_residualized.scaled) > 2) > 500


pca_data <- data.frame(
    Sample = rownames(pca_covariates),
    PC1 = pca_covariates[, 1],
    PC2 = pca_covariates[, 2],
    Batch = csf.meta$PlateId,
    Group = interaction(csf.meta$Gender, csf.meta$Age, csf.meta$Storage_time)
)

ggplot(pca_data, aes(x = PC1, y = PC2, color = csf.meta$Gender)) +
    geom_point(size = 3) +
    labs(
        title = "PCA Plot After Batch Effect Correction by Gender",
        x = "Principal Component 1",
        y = "Principal Component 2"
    ) +
    theme_minimal()
ggsave("pca_corrected_by_gender.pdf")

# PCA plot colored by Storage_time
ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(csf.meta$Storage_time))) +
    geom_point(size = 3) +
    labs(
        title = "PCA Plot After Batch Effect Correction by Storage Time",
        x = "Principal Component 1",
        y = "Principal Component 2"
    ) +
    theme_minimal()
ggsave("pca_corrected_by_storage_time.pdf")


ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(csf.meta$Age))) +
    geom_point(size = 3) +
    labs(
        title = "PCA plot after batch effect correction by Age",
        x = "PC1",
        y = "PC2"
    ) +
    theme_minimal()
ggsave("pca_corrected_by_age.pdf")

# remove global outliers and rerun
csf.mat.no_global_outliers <- csf.mat[!global_outliers, ]
csf_pca <- pca(csf.mat.no_global_outliers %>% t())

pc_ranges <- seq(5, 200, 5)
for (num_pc_correct in pc_ranges) {
    pca_covariates <- csf_pca$rotated[, num_pc_correct]
    csf.meta.no_global_outliers <- csf.meta[rownames(csf.mat.no_global_outliers), ]
    covariates.mat <- as.matrix(cbind(pca_covariates, csf.meta.no_global_outliers$Gender == "M", csf.meta.no_global_outliers$Age, csf.meta.no_global_outliers$Storage_time))
    csf.mat.no_global_outliers.residualized <- removeBatchEffect(csf.mat.no_global_outliers %>% t(), csf.meta.no_global_outliers$PlateId, covariates = covariates.mat)
    csf.mat.no_global_outliers.residualized.scaled <- csf.mat.no_global_outliers.residualized %>%
        t() %>%
        scale() %>%
        t()

    csf_residualized.scaled.LR_samples <- csf.mat.no_global_outliers.residualized.scaled[, csf.meta.no_global_outliers$PIDN %in% LR_samples]
    colnames(csf_residualized.scaled.LR_samples) <- csf.meta.no_global_outliers$PIDN[csf.meta.no_global_outliers$PIDN %in% LR_samples]
    colnames(csf_residualized.scaled.LR_samples)

    csf_zscores <- csf_residualized.scaled.LR_samples %>% melt()
    colnames(csf_zscores) <- c("SomaID", "PIDN", "zscore")

    ### biomaRt conversion
    protein.meta <- fread("~/tannerj/ADRComics/somascan_csf_Proteomics/ProteinMetadata_with_LOD1.csv")
    protein.meta$EntrezGeneSymbol
    csf_zscores %<>% left_join(protein.meta %>% dplyr::select(Key_2, EntrezGeneSymbol) %>% rename(SomaID = Key_2))
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_mapping <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "hgnc_symbol", values = csf_zscores$EntrezGeneSymbol, mart = mart)
    csf_zscores_withENSG <- csf_zscores %<>% left_join(gene_mapping, by = c("EntrezGeneSymbol" = "hgnc_symbol")) %<>% distinct()
    csf_zscores_withENSG <- csf_zscores_withENSG %>% drop_na()
    colnames(csf_zscores_withENSG) <- c("SomaID", "PIDN", "zscore", "gene_symbol", "gene_id")
    csf_zscores_withENSG <- csf_zscores_withENSG %>%
        mutate(outlier_status = case_when(
            zscore > 2 ~ 1,
            zscore < -2 ~ -1,
            TRUE ~ 0
        ))
    csf_zscores_withENSG$has_csf_protein_outlier <- abs(csf_zscores_withENSG$outlier_status)
    output_filename <- paste0("~/syang/ADRC/copy_unprocessed_data/output_filtered_csf_prot_pc_", num_pc_correct, ".csv")
    fwrite(csf_zscores_withENSG, file = output_filename, col.names = T, sep = "\t", quote = F)
}


### RNA expression processing: RNA reads for different pcs has been created via create_pseudobulk_profiles.R script


### now assign binary labels to omic outliers

library(rtracklayer)
library(data.table)
library(tidyverse)
library(plyranges)
library(magrittr)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
z_score_threshold <- as.numeric(args[1])
input_file <- args[2]
output_file <- args[3]
rna <- fread(input_file)

rna$has_expression_outlier <- ifelse(abs(rna$zscore) > z_score_threshold, 1, 0)
rna_subset <- rna %>% select(gene_id, PIDN, zscore, has_expression_outlier)

fwrite(rna_subset, file = output_file, row.name = FALSE)






### Lastly outlier pair analysis for every ome pair: RNA-protein plasma, RNA-CSF plasma, run all combinations to find optimal pcs with maximum number of outliers
library(dplyr)
library(data.table)
library(epitools)

args <- commandArgs(trailingOnly = TRUE)
path1 <- args[1]
label1 <- args[2]
path2 <- args[3]
label2 <- args[4]
output_file <- args[5]



rename_columns <- function(data, type) {
    if (type == "ome1") {
        colnames(data) <- gsub("^has_.*_outlier$", "has_ome1_outlier", colnames(data))
    } else if (type == "ome2") {
        colnames(data) <- gsub("^has_.*_outlier$", "has_ome2_outlier", colnames(data))
    }

    return(data)
}


pairwise_comparisons <- function(dataset1, dataset2, combined_table, label1, label2) {
    common_genes <- intersect(dataset1[["gene_id"]], dataset2[["gene_id"]])

    copy_result <- combined_table[combined_table$gene_id %in% common_genes, ]

    # Create a count table for the selected outlier columns
    copy_result_count_table <- copy_result[, c("gene_id", "PIDN", "has_ome1_outlier", "has_ome2_outlier")]
    copy_result_count_table[is.na(copy_result_count_table)] <- 0

    # Double positive cases
    copy_result_count_table_double_positive <- copy_result_count_table[copy_result_count_table[["has_ome1_outlier"]] == 1 & copy_result_count_table[["has_ome2_outlier"]] == 1, ]

    # Summarize counts
    summary_counts <- table(abs(copy_result[["has_ome1_outlier"]]) == 1, abs(copy_result[["has_ome2_outlier"]]) == 1)
    print(summary_counts)

    # Collapsing by gene_id and ID, determining whether there's any outlier for each modality
    gene_collapsed <- copy_result[, .(outlier1 = any(get("has_ome1_outlier")), outlier2 = any(get("has_ome2_outlier"))), by = .(gene_id, PIDN)]
    print(gene_collapsed)

    # Filter for genes with outliers in both modalities
    genes_to_test <- gene_collapsed[outlier1 == TRUE & outlier2 == TRUE, unique(gene_id)]
    gene_collapsed_filtered <- gene_collapsed[gene_id %in% genes_to_test]

    summary_counts_double_positive <- table(gene_collapsed_filtered$outlier1 == TRUE, gene_collapsed_filtered$outlier2 == TRUE)

    # Double positive summary counts
    if (nrow(gene_collapsed_filtered) == 0 || all(is.na(gene_collapsed_filtered$outlier1)) || all(is.na(gene_collapsed_filtered$outlier2))) {
        return(list(
            label1 = label1,
            label2 = label2,
            summary_counts = NA,
            summary_counts_double_positive = NA,
            glm_result = NA,
            glm_result_reversed = NA,
            fisher_result = NA,
            fisher_result_double_positive = NA
        ))
    }

    # Logistic regression models
    glm_result <- tryCatch(
        glm(outlier1 ~ outlier2, family = "binomial", data = gene_collapsed_filtered),
        error = function(e) {
            return(NA)
        }
    )

    glm_result_reversed <- tryCatch(
        glm(outlier2 ~ outlier1, family = "binomial", data = gene_collapsed_filtered),
        error = function(e) {
            return(NA)
        }
    )

    # Fisher's exact test
    fisher_result <- tryCatch(
        fisher.test(summary_counts),
        error = function(e) {
            return(NA)
        }
    )

    fisher_result_double_positive <- tryCatch(
        fisher.test(summary_counts_double_positive),
        error = function(e) {
            return(NA)
        }
    )

    return(list(
        label1 = label1,
        label2 = label2,
        summary_counts = summary_counts,
        summary_counts_double_positive = summary_counts_double_positive,
        glm_result = glm_result,
        glm_result_reversed = glm_result_reversed,
        fisher_result = fisher_result,
        fisher_result_double_positive = fisher_result_double_positive
    ))
}



print("reading files")
ome1_data <- fread(path1, sep = ",")
label_ome1 <- label1
ome2_data <- fread(path2, sep = ",")
label_ome2 <- label2

print("finished reading files")


# Rename columns for ome1 and ome2 dynamically
ome1_data <- rename_columns(ome1_data, "ome1") # Always rename current ome1 as ome1
ome2_data <- rename_columns(ome2_data, "ome2") # Always rename current ome2 as ome2

print("generalize to ome1 and ome2")


ome1_data <- unique(ome1_data, by = c("gene_id", "PIDN"))
ome2_data <- unique(ome2_data, by = c("gene_id", "PIDN"))

ome2_data$PIDN <- as.integer(ome2_data$PIDN)
ome1_data$PIDN <- as.integer(ome1_data$PIDN)
# Merge the datasets
result <- merge(ome1_data, ome2_data, by = c("gene_id", "PIDN"), all = TRUE)
result_copy <- as.data.table(result)

# Replace NA values with 0 for outlier columns
result_copy$has_ome1_outlier[is.na(result_copy$has_ome1_outlier)] <- 0
result_copy$has_ome2_outlier[is.na(result_copy$has_ome2_outlier)] <- 0


result_copy[, any(has_ome1_outlier > 0 | has_ome2_outlier > 0), by = gene_id]
# Perform pairwise comparison
comparison_result <- pairwise_comparisons(ome1_data, ome2_data, result_copy, label1, label2)

saveRDS(comparison_result, output_file)
