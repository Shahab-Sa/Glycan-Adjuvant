#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT:          BCR Diversity & Gene Usage Engine
# PURPOSE:         Accepts command line arguments to process pairs of BCR datasets, 
#                  calculating Alpha Diversity, Abundance, and V-Gene Fold-Change.
#
# AUTHOR:          Shahab Saghaei (Research Fellow & Data Manager)
# AFFILIATIONS:    Harvard Medical School | Broad Institute | Ragon Institute
# CONTACT:         https://www.linkedin.com/in/shahabsa/
#
# PUBLICATION:     "A glycan-based adjuvant expands the breadth and duration of 
#                  protection of mRNA-based vaccines." (Nature Immunology)
# CONTAINER IMAGE: Immcantation Framework (docker.io/immcantation/suite:4.5.0)
# LICENSE:         BSD 3-Clause
# ==============================================================================

suppressPackageStartupMessages({
  library(airr)
  library(alakazam)
  library(dplyr)
  library(ggplot2)
  library(tools)
  library(tidyr)
})

# Prevent R from generating a default Rplots.pdf file during non-interactive execution
pdf(NULL)

# --- 1. PARSE COMMAND LINE ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Requires exactly 6 arguments: file_1 file_2 name_1 name_2 main_out_dir dataset_prefix")
}

file_1         <- args[1]
file_2         <- args[2]
name_1         <- args[3]
name_2         <- args[4]
main_out_dir   <- args[5]
dataset_prefix <- toupper(args[6]) # Ensure cohort names conform to uppercase standard

# --- 2. CREATE SUBFOLDERS (Organized by Analysis Type) ---
base_dir  <- file.path(main_out_dir, "bcr_summary", dataset_prefix)
vgene_dir <- file.path(main_out_dir, "v_gene_usage", dataset_prefix)
div_dir   <- file.path(main_out_dir, "diversity_analysis", dataset_prefix)

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vgene_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(div_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# EXPLICIT COLOR PALETTE
# We assign specific hex codes directly to the conditions passed via CLI 
# to guarantee that the name_2 is always red (#F8766D) and  
# the name_1 is always teal (#00BFC4) across all diversity plots, 
# -------------------------------------------------------------------------
my_colors <- setNames(c("#00BFC4", "#F8766D"), c(name_1, name_2))

# Function: enhance_plot
# Purpose: Adjusts the line thickness and ribbon transparency for Alakazam 
#          diversity plots to improve aesthetic readability in final publications.
enhance_plot <- function(plot_obj, line_size = 1.5, ribbon_alpha = 0.1) {
  for (i in seq_along(plot_obj$layers)) {
    if (inherits(plot_obj$layers[[i]]$geom, "GeomRibbon")) {
      plot_obj$layers[[i]]$aes_params$alpha <- ribbon_alpha
    }
    if (inherits(plot_obj$layers[[i]]$geom, "GeomLine")) {
      plot_obj$layers[[i]]$aes_params$size <- line_size
    }
  }
  return(plot_obj)
}

# --- 3. DATA LOADING & PRE-PROCESSING ---
cat("\nLoading and filtering data...\n")

db_1 <- airr::read_rearrangement(file_1) %>% 
  mutate(condition = name_1) %>%
  mutate(cell_id = paste(condition, cell_id, sep = "_")) %>%
  filter(locus == "IGH")

db_2 <- airr::read_rearrangement(file_2) %>% 
  mutate(condition = name_2) %>%
  mutate(cell_id = paste(condition, cell_id, sep = "_")) %>%
  filter(locus == "IGH")

bcr_data <- bind_rows(db_1, db_2)
cat(sprintf("Total IGH rows before productive filter: %d\n", nrow(bcr_data)))

bcr_data <- bcr_data %>% filter(productive)
cat(sprintf("Total rows after productive filter: %d\n", nrow(bcr_data)))

# Filter out cells expressing multiple heavy chains to preserve strict single-cell confidence
multi_heavy <- table(filter(bcr_data, locus == "IGH")$cell_id)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]
bcr_data <- filter(bcr_data, !cell_id %in% multi_heavy_cells)
cat(sprintf("Total rows after filtering multiple heavy chains: %d\n", nrow(bcr_data)))

total_seqs <- bcr_data %>% group_by(condition) %>% summarise(total_sequences = n())
print(total_seqs)

# --- 4. CLONE SIZE FREQUENCIES ---
cat("\nComputing clone sizes...\n")
size_freq_wide <- bcr_data %>%
  group_by(condition, clone_id) %>%
  summarise(size = n(), .groups = "drop") %>%
  group_by(condition, size) %>%
  summarise(num_clones = n(), .groups = "drop") %>%
  arrange(size) %>%
  pivot_wider(names_from = condition, values_from = num_clones, values_fill = 0) %>%
  rename_with(~ paste0("num_clones_", .x), -size) %>%
  select(size, starts_with("num_clones_"))

write.csv(size_freq_wide, file.path(base_dir, "clone_size_frequencies_wide.csv"), row.names = FALSE)

# --- 5. ABUNDANCE ANALYSIS ---
cat(sprintf("\nRunning Abundance Estimation...\n"))
abund <- estimateAbundance(bcr_data, group = "condition", nboot = 10000)

write.table(abund@abundance, file = file.path(div_dir, "abund_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

abund_plot <- plot(abund, colors = my_colors, silent = TRUE) %>%
  enhance_plot(line_size = 1.5, ribbon_alpha = 0.6) +
  theme(
    axis.text.x = element_text(size = 32), axis.text.y = element_text(size = 32),
    axis.title  = element_text(size = 38), legend.title = element_text(size = 32),
    legend.text  = element_text(size = 32), plot.title = element_text(size = 40, hjust = 0.5)
  )
ggsave(file.path(div_dir, "abund_plot.jpg"), abund_plot, width = 10, height = 8, dpi = 300)

# --- 6. V-GENE USAGE & HEATMAP ---
cat("\nAnalyzing V-Gene Usage...\n")
gene_usage_v <- countGenes(bcr_data, gene="v_call", groups="condition", mode="gene") %>% ungroup()

# -------------------------------------------------------------------------
# DENOMINATOR ALIGNMENT
# Downstream heatmap functions automatically coerce group names to factors, 
# forcing an alphabetical visual hierarchy. To ensure our numerical Fold-Change 
# matrices perfectly sync with the visual orientation of those plots, 
# we dynamically extract the reference condition using identical factor coercion.
# -------------------------------------------------------------------------
pipeline_factors <- levels(as.factor(gene_usage_v$condition))
ref_condition <- pipeline_factors[1]   
comp_condition <- pipeline_factors[2]  

gene_usage_wide <- gene_usage_v %>%
  select(gene, condition, seq_freq) %>%
  pivot_wider(names_from = condition, values_from = seq_freq, values_fill = 0)

gene_fc <- gene_usage_wide %>%
  mutate(log2fc = log2((!!sym(comp_condition) + 1e-6)/(!!sym(ref_condition) + 1e-6)))

top_pos <- gene_fc %>% filter(log2fc > 0) %>% arrange(desc(log2fc)) %>% slice_head(n = 25)
top_neg <- gene_fc %>% filter(log2fc < 0) %>% arrange(desc(log2fc)) %>% slice_tail(n = 25)
top_genes_fc <- c(top_pos$gene, top_neg$gene)

if(any(duplicated(top_genes_fc))) {
    cat("Warning: Duplicated genes found in Top Fold Change selection.\n")
}

gene_usage_top <- gene_usage_v %>% filter(gene %in% top_genes_fc)
gene_usage_top$gene <- factor(gene_usage_top$gene, levels = rev(top_genes_fc))

fc_labels <- gene_fc %>%
  filter(gene %in% top_genes_fc) %>%
  mutate(dummy_x = max(as.numeric(as.factor(gene_usage_top$condition))) + 1,
         gene = factor(gene, levels = rev(top_genes_fc)))

base_font_pt <- 10
geom_text_size <- base_font_pt / ggplot2::.pt

p_heatmap_fc_label <- ggplot() +
  geom_tile(data = gene_usage_top, aes(x = condition, y = gene, fill = seq_freq),
            color = "white", linewidth = 0.1) +
  geom_text(data = fc_labels,
            aes(x = dummy_x, y = gene, label = sprintf("%.3f", log2fc)),
            hjust = 0, size = geom_text_size, family = "sans", fontface = "plain") +
  scale_fill_gradient(low = "white", high = "red",
                      labels = scales::percent_format(accuracy = 0.1)) +
  scale_x_discrete(expand = expansion(add = c(0.5, 3))) +
  theme_minimal(base_size = base_font_pt, base_family = "sans") +
  theme(
    axis.text.y = element_text(size = base_font_pt, family = "sans", face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = base_font_pt, family = "sans", face = "plain"),
    axis.title = element_blank(), panel.grid = element_blank(),
    legend.title = element_text(size = base_font_pt, family = "sans", face = "plain"),
    legend.text = element_text(size = base_font_pt, family = "sans", face = "plain"),
    plot.margin = margin(t = 20, r = 80, b = 10, l = 10)
  ) +
  labs(title = paste0("Heatmap of V-Gene Usage (", comp_condition, " vs ", ref_condition, ")"),
       fill = "Frequency") +
  coord_fixed(ratio = 0.5)

ggsave(file.path(vgene_dir, "heatmap_top25pos_top25neg_fc_final.jpg"), plot = p_heatmap_fc_label, width = 12, height = 14, dpi = 300)
write.table(gene_usage_v, file = file.path(vgene_dir, "gene_usage_freq.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_fc, file = file.path(vgene_dir, "gene_usage_rank_fc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- 7. ALPHA DIVERSITY ---
cat("\nRunning Alpha Diversity Comparisons...\n")

# Non-uniform (no rarefaction)
div_curve <- alphaDiversity(bcr_data, group = "condition", clone = "clone_id", 
                            min_q = 0, max_q = 4, step_q = 0.1, 
                            ci = 0.95, nboot = 10000, uniform = FALSE)

write.csv(div_curve@diversity, file.path(div_dir, "alpha_diversity_hill_curve_non_uniform.csv"), row.names = FALSE)

p_div_curve <- plotDiversityCurve(div_curve, colors = my_colors, legend_title = "Condition") %>%
  enhance_plot(line_size = 1.5, ribbon_alpha = 0.6) +
  theme(
    axis.text.x = element_text(size = 32), axis.text.y = element_text(size = 32),
    axis.title  = element_text(size = 38), legend.title = element_text(size = 32),
    legend.text  = element_text(size = 32), plot.title = element_text(size = 40, hjust = 0.5)
  )
ggsave(file.path(div_dir, "alpha_diversity_hill_curve_non_uniform.jpg"), plot = p_div_curve, width = 10, height = 8, dpi = 300)

# Uniform (rarefaction)
rare_curve <- alphaDiversity(bcr_data, group = "condition", clone = "clone_id", 
                             min_q = 0, max_q = 4, step_q = 0.1, 
                             ci = 0.95, nboot = 10000, uniform = TRUE)

write.csv(rare_curve@diversity, file.path(div_dir, "alpha_diversity_hill_curve_uniform.csv"), row.names = FALSE)

p_rare_curve <- plotDiversityCurve(rare_curve, colors = my_colors, legend_title = "Condition") %>%
  enhance_plot(line_size = 1.5, ribbon_alpha = 0.6) +
  theme(
    axis.text.x = element_text(size = 32), axis.text.y = element_text(size = 32),
    axis.title  = element_text(size = 38), legend.title = element_text(size = 32),
    legend.text  = element_text(size = 32), plot.title = element_text(size = 40, hjust = 0.5)
  )
ggsave(file.path(div_dir, "alpha_diversity_hill_curve_uniform.jpg"), plot = p_rare_curve, width = 10, height = 8, dpi = 300)

# --- 8. EXPORT STATISTICS ---
cat("\n--- Diversity Statistics (Rarefied) ---\n")
options(tibble.print_max = Inf, tibble.print_min = Inf)

cat("Full Hill Diversity (rare_curve) Tests Table:\n")
full_tests_rare <- rare_curve@tests %>% mutate(across(c(delta_mean, delta_sd, pvalue), ~ sprintf("%.4f", .)))
print(full_tests_rare, row.names = FALSE)

inverse_simpson <- bcr_data %>% group_by(condition) %>% 
  summarise(inverse_simpson = 1 / sum((table(clone_id) / n())^2)) %>%
  mutate(inverse_simpson = sprintf("%.2f", inverse_simpson))
cat("\nInverse Simpson Index (q=2):\n")
print(inverse_simpson, row.names = FALSE)

shannon <- bcr_data %>% group_by(condition) %>% 
  summarise(shannon = -sum((table(clone_id) / n()) * log(table(clone_id) / n()))) %>%
  mutate(shannon = sprintf("%.3f", shannon))
cat("\nShannon Entropy (q=1 approx):\n")
print(shannon, row.names = FALSE)

non_inverse_simpson <- bcr_data %>% group_by(condition) %>% 
  summarise(non_inverse_simpson = sum((table(clone_id) / n())^2)) %>%
  mutate(non_inverse_simpson = sprintf("%.5f", non_inverse_simpson))
cat("\nNon-Inverse Simpson Index:\n")
print(non_inverse_simpson, row.names = FALSE)

cat("\nDone! All files successfully saved.\n")
