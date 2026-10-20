library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)
library(forcats)
library(ggforce)

set.seed(123)

# Ordinates a dissimilarity object (PCoA, with a Lingoes correction for
# negative eigenvalues) and returns its first n_axes as a labelled data
# frame -- used to turn each community-distance matrix into a small set of
# numeric predictor axes for varpart() below
pcoa_to_df <- function(dist_obj, prefix = "Axis", n_axes = 5) {
  pc <- vegan::wcmdscale(as.dist(dist_obj), eig = TRUE, add = "lingoes")
  positive <- sum(pc$eig > 0, na.rm = TRUE)
  keep <- seq_len(min(n_axes, ncol(pc$points), positive))
  axs <- as.data.frame(pc$points[, keep, drop = FALSE])
  colnames(axs) <- paste0(prefix, seq_along(keep))
  tibble::rownames_to_column(axs, "ugs")
}

# Fixed layout for the 3-circle Venn diagram -- identical for every tissue x method panel, since only the R2 numbers inside change
venn_positions.df <- tibble(x = c(0, 1.15, 0.575), y = c(0.35, 0.35, 1.25),
                            group = c("Human traits", "Ecological traits", "Phylogeny"))

# Fixed label positions for the 7 variance-partition fractions -- same reasoning as above
label_positions.df <- tibble(fraction = c("[H] Human traits unique", "[O] Ecological traits unique", "[P] Phylogeny unique",
                                          "[H\u2229O]", "[O\u2229P]", "[H\u2229P]", "[H\u2229O\u2229P]"),
                             x = c(-0.45, 1.60, 0.575, 0.575, 0, 1.15, 0.575),
                             y = c(0.15, 0.15, 1.55, -0.05, 0.92, 0.92, 0.6))

# ---- 1. Load input data --------------------------------------------------

# site x chemical-feature presence/absence matrices (each has a "ugs" site-ID column)
leaf.chem.mat   <- read.csv("02_processed/leaf_site_mtb_mat.csv",   check.names = FALSE)
flower.chem.mat <- read.csv("02_processed/flower_site_mtb_mat.csv", check.names = FALSE)

species.mat   <- read.csv("01_raw/flora_mat.csv", check.names = FALSE) # site x plant-species matrix (only used to identify surveyed sites)
covariates.df <- read.csv("01_raw/covariates.csv")                       # only its "ugs" column is used below, for site alignment

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

# --- 2b. Reduce each of the three predictor dissimilarities (human traits, ecological traits, phylogeny) to their first 5 PCoA axes, for both the MPD- and MNTD-based versions ---
leaf.human_pcoa_mpd <- pcoa_to_df(leaf.human.dist.comdist, prefix = "Human",  n_axes = 5)
leaf.eco_pcoa_mpd <- pcoa_to_df(leaf.eco.dist.comdist,   prefix = "Others", n_axes = 5)
leaf.phylo_pcoa_mpd <- pcoa_to_df(leaf.phylo.dist.comdist, prefix = "Phylo",  n_axes = 5)

leaf.human_pcoa_mntd <- pcoa_to_df(leaf.human.dist.comdistnt, prefix = "Human",  n_axes = 5)
leaf.eco_pcoa_mntd <- pcoa_to_df(leaf.eco.dist.comdistnt,   prefix = "Others", n_axes = 5)
leaf.phylo_pcoa_mntd <- pcoa_to_df(leaf.phylo.dist.comdistnt, prefix = "Phylo",  n_axes = 5)

# --- 2c. Variation partitioning: how much of chemical dissimilarity is uniquely/jointly explained by each predictor block ---
leaf.X_human_mpd <- leaf.human_pcoa_mpd %>% column_to_rownames("ugs")
leaf.X_eco_mpd <- leaf.eco_pcoa_mpd %>% column_to_rownames("ugs")
leaf.X_phylo_mpd <- leaf.phylo_pcoa_mpd %>% column_to_rownames("ugs")

leaf.X_human_mntd <- leaf.human_pcoa_mntd %>% column_to_rownames("ugs")
leaf.X_eco_mntd <- leaf.eco_pcoa_mntd %>% column_to_rownames("ugs")
leaf.X_phylo_mntd <- leaf.phylo_pcoa_mntd %>% column_to_rownames("ugs")

leaf.vp_mpd <- vegan::varpart(leaf.chem.dist, leaf.X_human_mpd,  leaf.X_eco_mpd,  leaf.X_phylo_mpd)
leaf.vp_mntd <- vegan::varpart(leaf.chem.dist, leaf.X_human_mntd, leaf.X_eco_mntd, leaf.X_phylo_mntd)

# vegan::varpart's indfract table lists, in order: unique-to-block-1, unique-to-block-2, unique-to-block-3 (rows 1:3), then the 4 shared fractions (rows 4:7) 
leaf.vp_df_mpd <- tibble(fraction = label_positions.df$fraction,
                         R2_adj = c(leaf.vp_mpd$part$indfract$Adj.R.square[1:3], leaf.vp_mpd$part$indfract$Adj.R.square[4:7])) %>%
  dplyr::mutate(R2_plot = pmax(R2_adj, 0), fraction = forcats::fct_inorder(fraction))

leaf.vp_df_mntd <- tibble(fraction = label_positions.df$fraction,
                          R2_adj = c(leaf.vp_mntd$part$indfract$Adj.R.square[1:3], leaf.vp_mntd$part$indfract$Adj.R.square[4:7])) %>%
  dplyr::mutate(R2_plot = pmax(R2_adj, 0), fraction = forcats::fct_inorder(fraction))

################################################################################
# ---- 3. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 3a. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.chem.dist <- vegan::vegdist(flower.chem.mat %>% column_to_rownames("ugs"), method = "jaccard", binary = TRUE)

# --- 3b. Reduce each predictor dissimilarity to its first 5 PCoA axes ---
flower.human_pcoa_mpd <- pcoa_to_df(flower.human.dist.comdist, prefix = "Human",  n_axes = 5)
flower.eco_pcoa_mpd <- pcoa_to_df(flower.eco.dist.comdist, prefix = "Others", n_axes = 5)
flower.phylo_pcoa_mpd <- pcoa_to_df(flower.phylo.dist.comdist, prefix = "Phylo",  n_axes = 5)

flower.human_pcoa_mntd <- pcoa_to_df(flower.human.dist.comdistnt, prefix = "Human",  n_axes = 5)
flower.eco_pcoa_mntd <- pcoa_to_df(flower.eco.dist.comdistnt, prefix = "Others", n_axes = 5)
flower.phylo_pcoa_mntd <- pcoa_to_df(flower.phylo.dist.comdistnt, prefix = "Phylo",  n_axes = 5)

# --- 3c. Variation partitioning ---
flower.X_human_mpd <- flower.human_pcoa_mpd %>% column_to_rownames("ugs")
flower.X_eco_mpd <- flower.eco_pcoa_mpd   %>% column_to_rownames("ugs")
flower.X_phylo_mpd <- flower.phylo_pcoa_mpd %>% column_to_rownames("ugs")

flower.X_human_mntd <- flower.human_pcoa_mntd %>% column_to_rownames("ugs")
flower.X_eco_mntd <- flower.eco_pcoa_mntd   %>% column_to_rownames("ugs")
flower.X_phylo_mntd <- flower.phylo_pcoa_mntd %>% column_to_rownames("ugs")

flower.vp_mpd  <- vegan::varpart(flower.chem.dist, flower.X_human_mpd,  flower.X_eco_mpd,  flower.X_phylo_mpd)
flower.vp_mntd <- vegan::varpart(flower.chem.dist, flower.X_human_mntd, flower.X_eco_mntd, flower.X_phylo_mntd)

flower.vp_df_mpd <- tibble(fraction = label_positions.df$fraction,
                           R2_adj = c(flower.vp_mpd$part$indfract$Adj.R.square[1:3], flower.vp_mpd$part$indfract$Adj.R.square[4:7])) %>%
  dplyr::mutate(R2_plot = pmax(R2_adj, 0), fraction = forcats::fct_inorder(fraction))

flower.vp_df_mntd <- tibble(fraction = label_positions.df$fraction,
                            R2_adj = c(flower.vp_mntd$part$indfract$Adj.R.square[1:3], flower.vp_mntd$part$indfract$Adj.R.square[4:7])) %>%
  dplyr::mutate(R2_plot = pmax(R2_adj, 0), fraction = forcats::fct_inorder(fraction))

# ---- 4. Assemble the plot data for all 4 panels (tissue x method) ----------

# Venn-circle positions: identical layout repeated for each of the 4 panels
venn_df <- bind_rows(data.frame(venn_positions.df, tissue = "Leaf",   method = "MPD"),
                     data.frame(venn_positions.df, tissue = "Leaf",   method = "MNTD"),
                     data.frame(venn_positions.df, tissue = "Flower", method = "MPD"),
                     data.frame(venn_positions.df, tissue = "Flower", method = "MNTD")) %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")),
         method = factor(method, levels = c("MPD", "MNTD")),
         group = factor(group, levels = c("Human traits", "Ecological traits", "Phylogeny")))

# Label positions + each panel's R2 values, joined onto the shared position template
label_df <- bind_rows(data.frame(label_positions.df %>% left_join(leaf.vp_df_mpd   %>% select(fraction, R2_adj), by = "fraction"), tissue = "Leaf",   method = "MPD"),
                      data.frame(label_positions.df %>% left_join(leaf.vp_df_mntd  %>% select(fraction, R2_adj), by = "fraction"), tissue = "Leaf",   method = "MNTD"),
                      data.frame(label_positions.df %>% left_join(flower.vp_df_mpd %>% select(fraction, R2_adj), by = "fraction"), tissue = "Flower", method = "MPD"),
                      data.frame(label_positions.df %>% left_join(flower.vp_df_mntd %>% select(fraction, R2_adj), by = "fraction"), tissue = "Flower", method = "MNTD")) %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")),
         method = factor(method, levels = c("MPD", "MNTD")))

# total variance explained (sum of all 7 fractions), annotated above each panel
total_var_df <- tibble(tissue = c("Leaf", "Leaf", "Flower", "Flower"),
                       method = c("MPD", "MNTD", "MPD", "MNTD"),
                       total_R2 = c(sum(leaf.vp_mpd$part$indfract$Adj.R.square[1:7], na.rm = TRUE),
                                    sum(leaf.vp_mntd$part$indfract$Adj.R.square[1:7], na.rm = TRUE),
                                    sum(flower.vp_mpd$part$indfract$Adj.R.square[1:7], na.rm = TRUE),
                                    sum(flower.vp_mntd$part$indfract$Adj.R.square[1:7], na.rm = TRUE))) %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")),
         method = factor(method, levels = c("MPD", "MNTD")),
         label = sprintf("Total explained = %.3f", total_R2),
         x = 0.575, y = 2.45)

# ---- 5. Build and save the figure -----------------------------------------

ggplot() +
  ggforce::geom_circle(data = venn_df, aes(x0 = x, y0 = y, r = 1, fill = group), alpha = 0.45, color = "grey20", linewidth = 0.5) +
  geom_text(data = label_df, aes(x = x, y = y, label = round(R2_adj, 3)), vjust = 0, size = 4) +
  geom_text(data = total_var_df, aes(x = x, y = y, label = label), size = 4) +
  facet_grid(tissue ~ method) +
  labs(fill = NULL) +
  paletteer::scale_fill_paletteer_d("ggthemes::Classic_Color_Blind") +
  coord_fixed(clip = "off") +
  theme_void(base_size = 12, paper = scales::alpha("white", 1)) +
  theme(strip.text = element_text(size = 12),
        strip.text.y = element_text(angle = -90),
        legend.position = "bottom")

ggsave("figures/figS5.png", height = 6, width = 8)
