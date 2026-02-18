####################
# title:"Watershed vs. GAM performance comparison"
# author:"Sherry Yang"
# date:"2026-02-16"
####################

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(optparse)


arguments <- parse_args(OptionParser(
  usage = "%prog [options]", description = "Watershed command line args",
  option_list = list(
    make_option(c("-e", "--eval_object"), default = NULL, help = "the evaluate rds model object"),
    make_option(c("-o", "--output_file_png"), default = NULL, help = "output model performance graph name")
  )
))

results <- arguments$eval_object
out_png <- arguments$output_file_png

aucs <- c()
gam_aucs <- c()

ws_obj <- readRDS(results)
expr_aucs <- c(aucs, ws_obj$auc[[1]]$evaROC$watershed_pr_auc)
expr_gam_aucs <- c(gam_aucs, ws_obj$auc[[1]]$evaROC$GAM_pr_auc)


df_ws_rna <- data.frame(AUC = expr_aucs, Outlier = "Expression", Model = "Watershed")
df_gam_rna <- data.frame(AUC = expr_gam_aucs, Outlier = "Expression", Model = "GAM")
df <- rbind(df_ws_rna, df_gam_rna) # df_ws_pro, df_gam_pro

pr_df_rna <- data.frame(
  recall = ws_obj$auc[[1]]$evaROC$watershed_recall,
  precision = ws_obj$auc[[1]]$evaROC$watershed_precision,
  outlier = "RNA", model = sprintf("WS (AUC=%.3f)", ws_obj$auc[[1]]$evaROC$watershed_pr_auc)
)

pr_df_rna_gam <- data.frame(
  recall = ws_obj$auc[[1]]$evaROC$GAM_recall,
  precision = ws_obj$auc[[1]]$evaROC$GAM_precision,
  outlier = "RNA", model = sprintf("GAM (AUC=%.3f)", AUC = ws_obj$auc[[1]]$evaROC$GAM_pr_auc)
)

pr_df <- rbind(pr_df_rna_gam, pr_df_rna)
pr_df <- pr_df %>% mutate(name = sprintf("%s %s", outlier, model))

p <- ggplot(pr_df, aes(x = recall, y = precision, color = name)) +
  geom_point(size = 1) +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.9), legend.background = element_rect(fill = "grey97"), legend.title = element_blank()) +
  theme(aspect.ratio = 2 / 3, text = element_text(size = 10)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  ggtitle("Precision/recall curve")

ggsave(out_png, plot = p, width = 6, height = 4)
