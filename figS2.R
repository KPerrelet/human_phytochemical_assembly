library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(vegan)

set.seed(123) # for reproducibility of the permutation tests (protest) below

# Converts a p-value into a conventional significance code
sig_label <- function(p) {
  dplyr::case_when(is.na(p) ~ "",
                   p < 0.001 ~ "***",
                   p < 0.01 ~ "**",
                   p < 0.05 ~ "*",
                   p < 0.1 ~ ".",
                   TRUE ~ "ns")
}

# Color palette distinguishing the two tissues
tissue.colors <- c("Leaf" = "#007C91", "Flower" = "#8F2D56")


# ---- 1. Load input data --------------------------------------------------

# site x chemical-feature presence/absence matrices
# To recreate these datasets, use the prepare_metabolomics_df.R script. 
leaf.chem.mat <- read.csv("02_processed/leaf_site_mtb_mat.csv",   check.names = FALSE)
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

# --- 2b. PCoA (with a Lingoes correction for negative eigenvalues) on each
# of the four dissimilarity matrices ---
leaf.human.ord <- vegan::wcmdscale(leaf.human.dist.comdistnt, eig = TRUE, add = "lingoes")
leaf.eco.ord   <- vegan::wcmdscale(leaf.eco.dist.comdistnt,   eig = TRUE, add = "lingoes")
leaf.phylo.ord <- vegan::wcmdscale(leaf.phylo.dist.comdistnt, eig = TRUE, add = "lingoes")
leaf.chem.ord  <- vegan::wcmdscale(leaf.chem.dist,  eig = TRUE, add = "lingoes")

# only sites with valid coordinates in all four ordinations can be compared
leaf.common.sites <- Reduce(intersect, list(rownames(leaf.human.ord$points), rownames(leaf.eco.ord$points),
                                            rownames(leaf.phylo.ord$points), rownames(leaf.chem.ord$points)))

leaf.H <- leaf.human.ord$points[leaf.common.sites, 1:2, drop = FALSE]
leaf.O <- leaf.eco.ord$points[leaf.common.sites, 1:2, drop = FALSE]
leaf.P <- leaf.phylo.ord$points[leaf.common.sites, 1:2, drop = FALSE]
leaf.M <- leaf.chem.ord$points[leaf.common.sites, 1:2, drop = FALSE]

# --- 2c. Pairwise Procrustes correlation tests (protest) between every pair of the four ordination spaces (Human traits H, Ecological traits O, Phylogeny P, chemical/Metabolomics M) ---
leaf.proc_HO <- vegan::protest(leaf.H, leaf.O, permutations = 9999)
leaf.proc_HP <- vegan::protest(leaf.H, leaf.P, permutations = 9999)
leaf.proc_OP <- vegan::protest(leaf.O, leaf.P, permutations = 9999)
leaf.proc_PM <- vegan::protest(leaf.P, leaf.M, permutations = 9999)
leaf.proc_HM <- vegan::protest(leaf.H, leaf.M, permutations = 9999)
leaf.proc_OM <- vegan::protest(leaf.O, leaf.M, permutations = 9999)

leaf.procrustes.df <- tibble(comparison = c("Human traits\n-Ecological traits", "Human traits\n-phylogeny", "Ecological traits\n-phylogeny",
                                            "Phylogeny\n-metabolomics", "Human traits\n-metabolomics", "Ecological traits\n-metabolomics"),
                             r = c(leaf.proc_HO$t0, leaf.proc_HP$t0, leaf.proc_OP$t0, leaf.proc_PM$t0, leaf.proc_HM$t0, leaf.proc_OM$t0),
                             p_val = c(leaf.proc_HO$signif, leaf.proc_HP$signif, leaf.proc_OP$signif, leaf.proc_PM$signif, leaf.proc_HM$signif, leaf.proc_OM$signif)) %>%
  mutate(sig = sig_label(p_val))

################################################################################
# ---- 3. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 3a. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.chem.dist <- vegan::vegdist(flower.chem.mat %>% column_to_rownames("ugs"), method = "jaccard", binary = TRUE)

# --- 3b. PCoA on each of the four dissimilarity matrices ---
flower.human.ord <- vegan::wcmdscale(flower.human.dist.comdistnt, eig = TRUE, add = "lingoes")
flower.eco.ord   <- vegan::wcmdscale(flower.eco.dist.comdistnt,   eig = TRUE, add = "lingoes")
flower.phylo.ord <- vegan::wcmdscale(flower.phylo.dist.comdistnt, eig = TRUE, add = "lingoes")
flower.chem.ord  <- vegan::wcmdscale(flower.chem.dist,  eig = TRUE, add = "lingoes")

flower.common.sites <- Reduce(intersect, list(rownames(flower.human.ord$points), rownames(flower.eco.ord$points),
                                              rownames(flower.phylo.ord$points), rownames(flower.chem.ord$points)))

flower.H <- flower.human.ord$points[flower.common.sites, 1:2, drop = FALSE]
flower.O <- flower.eco.ord$points[flower.common.sites, 1:2, drop = FALSE]
flower.P <- flower.phylo.ord$points[flower.common.sites, 1:2, drop = FALSE]
flower.M <- flower.chem.ord$points[flower.common.sites, 1:2, drop = FALSE]

# --- 3c. Pairwise Procrustes correlation tests ---
flower.proc_HO <- vegan::protest(flower.H, flower.O, permutations = 9999)
flower.proc_HP <- vegan::protest(flower.H, flower.P, permutations = 9999)
flower.proc_OP <- vegan::protest(flower.O, flower.P, permutations = 9999)
flower.proc_PM <- vegan::protest(flower.P, flower.M, permutations = 9999)
flower.proc_HM <- vegan::protest(flower.H, flower.M, permutations = 9999)
flower.proc_OM <- vegan::protest(flower.O, flower.M, permutations = 9999)

flower.procrustes.df <- tibble(comparison = c("Human traits\n-Ecological traits", "Human traits\n-phylogeny", "Ecological traits\n-phylogeny",
                                              "Phylogeny\n-metabolomics", "Human traits\n-metabolomics", "Ecological traits\n-metabolomics"),
                               r = c(flower.proc_HO$t0, flower.proc_HP$t0, flower.proc_OP$t0, flower.proc_PM$t0, flower.proc_HM$t0, flower.proc_OM$t0),
                               p_val = c(flower.proc_HO$signif, flower.proc_HP$signif, flower.proc_OP$signif, flower.proc_PM$signif, flower.proc_HM$signif, flower.proc_OM$signif)) %>%
  mutate(sig = sig_label(p_val))


# ---- 4. Combine tissues and build the figure --------------------------------

procrustes_results.df <- rbind(data.frame(leaf.procrustes.df, tissue = "Leaf"),
                               data.frame(flower.procrustes.df, tissue = "Flower"))

ggplot(procrustes_results.df, aes(x = comparison, y = r, color = tissue, fill = tissue)) +
  geom_col(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.2) +
  scale_color_manual(values = tissue.colors) +
  scale_fill_manual(values = tissue.colors) +
  labs(x = NULL, y = "R\u00b2", fill = NULL, color = NULL) +
  theme_pubclean(base_size = 12) +
  theme(legend.position = "top",
        axis.text = element_text(size = 12))

ggsave("figures/figS2.png", height = 5, width = 9)
