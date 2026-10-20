library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggnewscale)

set.seed(123) # for reproducibility of the GAM smooths fitted in the plot below

# Converts a Mantel test p-value into a conventional significance code
sig_label <- function(p) {
  dplyr::case_when(is.na(p) ~ "",
                   p < 0.001 ~ "***",
                   p < 0.01 ~ "**",
                   p < 0.05 ~ "*",
                   p < 0.1 ~ ".",
                   TRUE ~ "ns")
}

# Turns a set of same-sized square distance matrices into one "long" data
# frame with one row per site pair (the lower triangle of each matrix),
# so their values can be plotted against each other
dist_lower_df <- function(...) {
  mats <- list(...)
  stopifnot(length(mats) > 0)
  rn <- rownames(as.matrix(mats[[1]]))
  idx <- which(lower.tri(as.matrix(mats[[1]])), arr.ind = TRUE)
  out <- tibble(site_a = rn[idx[, 1]], site_b = rn[idx[, 2]])
  for (nm in names(mats)) {
    mat <- as.matrix(mats[[nm]])
    out[[nm]] <- mat[idx]
  }
  out
}

# Runs a Mantel test between two distance matrices and returns a tidy
# one-row summary (correlation, p-value, significance code)
run_mantel_tidy <- function(x, y, label, method = "pearson", permutations = 9999) {
  m <- vegan::mantel(as.dist(x), as.dist(y), method = method, permutations = permutations)
  tibble(test = label, r = unname(m$statistic), p_val = m$signif, sig = sig_label(m$signif))
}


# ---- 1. Load input data --------------------------------------------------

# site x chemical-feature presence/absence matrices
# To recreate these datasets, use the prepare_metabolomics_df.R script. 
leaf.chem.mat <- read.csv("02_processed/leaf_site_mtb_mat.csv", check.names = FALSE)
flower.chem.mat <- read.csv("02_processed/flower_site_mtb_mat.csv", check.names = FALSE)

species.mat <- read.csv("01_raw/flora_mat.csv", check.names = FALSE)   # site x plant-species presence/absence
covariates.df <- read.csv("01_raw/covariates.csv") # per-site covariates

# Precomputed community-weighted dissimilarities (human-relevant traits, ecological traits, phylogeny) for the leaf and flower sample sets. 
# These were built upstream (e.g. with picante / V.PhyloMaker2) and are simply loaded here as ready-made distance objects, then renamed for clarity.
# To recreate this dataset, use the phylo_funct_dist.R script. 
load("02_processed/phyto_funct_dist.RData")

################################################################################
# ---- 2. LEAF analysis -----------------------------------------------------
################################################################################

# --- 2a. Restrict to sites with chemistry, floristic, and covariate data ---
leaf.keep.sites <- Reduce(intersect, list(leaf.chem.mat$ugs, species.mat$ugs, covariates.df$ugs))
leaf.chem.mat <- leaf.chem.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))
leaf.chem.dist <- vegan::vegdist(leaf.chem.mat %>% column_to_rownames("ugs"), method = "jaccard", binary = TRUE)

# --- 2b. Mantel tests: is leaf chemical dissimilarity correlated with functional / phylogenetic dissimilarity? ---
leaf.mantel_results <- bind_rows(run_mantel_tidy(leaf.chem.dist, leaf.all.dist.comdist, "Metabolomic vs functional MPD"),
                                 run_mantel_tidy(leaf.chem.dist, leaf.all.dist.comdistnt, "Metabolomic vs functional MNTD"),
                                 run_mantel_tidy(leaf.chem.dist, leaf.phylo.dist.comdist, "Metabolomic vs phylogenetic MPD"),
                                 run_mantel_tidy(leaf.chem.dist, leaf.phylo.dist.comdistnt, "Metabolomic vs phylogenetic MNTD"))

# --- 2c. Pairwise site-by-site distances, for the scatterplots ---
leaf.pair_df <- dist_lower_df(all_trait_dist = leaf.all.dist.comdist,
                              all_trait_nt_dist = leaf.all.dist.comdistnt,
                              phylo_dist = leaf.phylo.dist.comdist,
                              phylo_nt_dist = leaf.phylo.dist.comdistnt,
                              dist = leaf.chem.dist) %>%
  left_join(covariates.df %>% select(ugs, ugs_type) %>% rename(site_a = ugs, type_a = ugs_type), by = "site_a") %>%
  left_join(covariates.df %>% select(ugs, ugs_type) %>% rename(site_b = ugs, type_b = ugs_type), by = "site_b") %>%
  mutate(pair_type = ifelse(type_a == type_b, type_a, "Between types"))

# reshape into one row per site-pair x predictor combination, for faceting
leaf.pair_plot_df <- bind_rows(leaf.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Functional MPD", x = all_trait_dist, dist),
                               leaf.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Phylogenetic MPD", x = phylo_dist, dist),
                               leaf.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Functional MNTD", x = all_trait_nt_dist, dist),
                               leaf.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Phylogenetic MNTD", x = phylo_nt_dist, dist)) %>%
  mutate(predictor = factor(predictor, levels = c("Functional MPD", "Phylogenetic MPD", "Functional MNTD", "Phylogenetic MNTD")))

################################################################################
# ---- 3. FLOWER analysis ----------------------------------------------------
################################################################################

# --- 3a. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.chem.dist <- vegan::vegdist(flower.chem.mat %>% column_to_rownames("ugs"), method = "jaccard", binary = TRUE)

# --- 3b. Mantel tests: is flower chemical dissimilarity correlated with functional / phylogenetic dissimilarity? ---
flower.mantel_results <- bind_rows(run_mantel_tidy(flower.chem.dist, flower.all.dist.comdist, "Metabolomic vs functional MPD"),
                                   run_mantel_tidy(flower.chem.dist, flower.all.dist.comdistnt, "Metabolomic vs functional MNTD"),
                                   run_mantel_tidy(flower.chem.dist, flower.phylo.dist.comdist, "Metabolomic vs phylogenetic MPD"),
                                   run_mantel_tidy(flower.chem.dist, flower.phylo.dist.comdistnt, "Metabolomic vs phylogenetic MNTD"))

# --- 3c. Pairwise site-by-site distances, for the scatterplots ---
flower.pair_df <- dist_lower_df(all_trait_dist = flower.all.dist.comdist,
                                all_trait_nt_dist = flower.all.dist.comdistnt,
                                phylo_dist = flower.phylo.dist.comdist,
                                phylo_nt_dist = flower.phylo.dist.comdistnt,
                                dist = flower.chem.dist) %>%
  left_join(covariates.df %>% select(ugs, ugs_type) %>% rename(site_a = ugs, type_a = ugs_type), by = "site_a") %>%
  left_join(covariates.df %>% select(ugs, ugs_type) %>% rename(site_b = ugs, type_b = ugs_type), by = "site_b") %>%
  mutate(pair_type = ifelse(type_a == type_b, type_a, "Between types"))

flower.pair_plot_df <- bind_rows(flower.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Functional MPD",    x = all_trait_dist,    dist),
                                 flower.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Phylogenetic MPD",  x = phylo_dist,        dist),
                                 flower.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Functional MNTD",   x = all_trait_nt_dist, dist),
                                 flower.pair_df %>% transmute(site_a, site_b, pair_type, predictor = "Phylogenetic MNTD", x = phylo_nt_dist,     dist)) %>%
  mutate(predictor = factor(predictor, levels = c("Functional MPD", "Phylogenetic MPD", "Functional MNTD", "Phylogenetic MNTD")))


# ---- 4. Combine tissues and build Mantel-test annotations -----------------

pair_plot_df <- bind_rows(leaf.pair_plot_df   %>% mutate(tissue = "Leaf"),
                          flower.pair_plot_df %>% mutate(tissue = "Flower"))

mantel_labels <- bind_rows(leaf.mantel_results   %>% mutate(tissue = "Leaf"),
                           flower.mantel_results %>% mutate(tissue = "Flower")) %>%
  mutate(label = paste0("r = ", round(r, 2),
                        "\np ", ifelse(p_val < 0.001, "< 0.001", paste0("= ", format.pval(p_val, digits = 2, eps = 0.001)))))

# map each Mantel test onto the predictor panel it annotates (case_when's fallback always returns a value here, so every row gets a predictor)
mantel_labels <- mantel_labels %>%
  mutate(predictor = case_when(test == "Metabolomic vs functional MPD" ~ "Functional MPD",
                               test == "Metabolomic vs functional MNTD" ~ "Functional MNTD",
                               test == "Metabolomic vs phylogenetic MPD" ~ "Phylogenetic MPD",
                               test == "Metabolomic vs phylogenetic MNTD" ~ "Phylogenetic MNTD",
                               TRUE ~ test), 
         predictor = factor(predictor, levels = c("Functional MPD", "Phylogenetic MPD", "Functional MNTD", "Phylogenetic MNTD"))) 


# ---- 5. Build and save the figure -----------------------------------------

ggplot(pair_plot_df, aes(x = x, y = dist)) +
  geom_point(aes(color = pair_type), alpha = 0.35, size = 1.2) +
  labs(color = "Site pair") +
  ggnewscale::new_scale_color() +
  geom_smooth(aes(color = "GAM"), method = "gam", formula = y ~ s(x, k = 4), linewidth = 1, se = TRUE) +
  geom_smooth(aes(color = "LM"), method = "lm", linewidth = 1, se = TRUE) +
  geom_text(data = mantel_labels, aes(x = Inf, y = -Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1,
            vjust = -0.5,
            size = 4) +
  facet_grid(factor(tissue, levels = c("Leaf", "Flower")) ~ predictor, scales = "free_x") +
  scale_color_manual(values = c("GAM" = "black", "LM" = "blue")) +
  labs(x = "Pairwise distance", y = "Pairwise feature distance", color = "Model") +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 12),
        panel.spacing.x = unit(15, "pt"),
        legend.direction = "vertical", legend.position = "right",
        legend.text = element_text(size = 14),
        legend.key.spacing.y = unit(2, "pt"),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0), legend.margin = ggplot2::margin(0, 0, 0, 0))

ggsave("figures/fig3.png", height = 5, width = 12)
ggsave("figures/fig3.pdf", height = 5, width = 12)
