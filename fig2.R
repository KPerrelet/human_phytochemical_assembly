library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(vegan)

set.seed(123) # for reproducibility of the ordisurf() surface fits below

# Color palette for the five site land-use categories
site.colors <- c("Allotment" = "#FEC601",
                 "Garden" = "#072AC8",
                 "Park" = "#594F3B",
                 "Real Estate" = "#3DA5D9",
                 "Ruderal" = "#DD403A")

# Point shapes used to distinguish the four ordination "spaces" (chemistry vs. the three spaces that get Procrustes-rotated onto it)
type.shapes <- c("Chemical feature" = 16,
                 "Human traits" = 24,
                 "Ecological traits"  = 25,
                 "Phylogeny" = 15)

# % of total variance explained by each positive PCoA axis of an ordination
pcoa_var <- function(ord, label) {
  eig <- ord$eig[ord$eig > 0]
  tibble(type = label,
         axis = paste0("PCoA", seq_along(eig)),
         variance_explained = eig / sum(eig) * 100)
}

# Builds an axis label like "Procrustes dimension 1 (12.34-56.78%)" summarising
# the range of variance explained by that axis across the four ordinations
axis_label <- function(pcoa.var.df, axis.name, dim.label) {
  vals <- pcoa.var.df$variance_explained[pcoa.var.df$axis == axis.name]
  sprintf("%s (%.2f\u2014%.2f%%)", dim.label, min(vals), max(vals))
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
leaf.chem.dist <- vegdist(leaf.chem.mat %>% column_to_rownames("ugs"), method = "jaccard", binary = TRUE)

# --- 2b. PCoA (with a Lingoes correction for negative eigenvalues) on each of the four dissimilarity matrices ---
leaf.human.ord <- wcmdscale(leaf.human.dist.comdistnt, eig = TRUE, add = "lingoes")
leaf.eco.ord <- wcmdscale(leaf.eco.dist.comdistnt,   eig = TRUE, add = "lingoes")
leaf.phylo.ord <- wcmdscale(leaf.phylo.dist.comdistnt, eig = TRUE, add = "lingoes")
leaf.chem.ord <- wcmdscale(leaf.chem.dist,  eig = TRUE, add = "lingoes")

# only sites with valid coordinates in all four ordinations can be compared
leaf.common.sites <- Reduce(intersect, list(rownames(leaf.human.ord$points), rownames(leaf.eco.ord$points),
                                            rownames(leaf.phylo.ord$points), rownames(leaf.chem.ord$points)))

# --- 2c. Procrustes rotation: rotate/rescale each trait or phylogeny ordination (first 2 axes) onto the chemical ordination, which acts as the common target space ---
leaf.human.proc <- procrustes(leaf.chem.ord$points[leaf.common.sites, 1:2], leaf.human.ord$points[leaf.common.sites, 1:2], symmetric = FALSE)
leaf.eco.proc <- procrustes(leaf.chem.ord$points[leaf.common.sites, 1:2], leaf.eco.ord$points[leaf.common.sites, 1:2],   symmetric = FALSE)
leaf.phylo.proc <- procrustes(leaf.chem.ord$points[leaf.common.sites, 1:2], leaf.phylo.ord$points[leaf.common.sites, 1:2], symmetric = FALSE)

# collect the rotated coordinates from each of the 4 spaces into one data frame, plus the site's land-use category (ugs_type)
leaf.aligned.df <- bind_rows(as.data.frame(leaf.human.proc$Yrot) %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Human traits"),
                             as.data.frame(leaf.eco.proc$Yrot) %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Ecological traits"),
                             as.data.frame(leaf.phylo.proc$Yrot) %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Phylogeny"),
                             as.data.frame(leaf.human.proc$X) %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Chemical feature")) %>% 
  left_join(covariates.df %>% select(ugs, ugs_type), by = "ugs") %>%
  mutate(type = factor(type, levels = c("Human traits", "Ecological traits", "Phylogeny", "Chemical feature")))

# --- 2d. Fit a smoothed response surface (ordisurf, a GAM-based method) of the human-influence index over the chemical-feature ordination space. This produces the background contour lines drawn on the panel. ---
leaf.chem.points <- leaf.aligned.df %>%
  filter(type == "Chemical feature") %>%
  left_join(covariates.df %>% select(ugs, human_index), by = "ugs") %>%
  filter(!is.na(human_index))

leaf.surf.fit <- ordisurf(leaf.chem.points[, c("Dim1", "Dim2")], leaf.chem.points$human_index, plot = FALSE)

leaf.surf.df <- expand.grid(Dim1 = leaf.surf.fit$grid$x, Dim2 = leaf.surf.fit$grid$y) %>%
  mutate(human_index = as.vector(leaf.surf.fit$grid$z), type = "Chemical feature")

# --- 2e. Line segments from each alternative ordination's point to its true position in chemical space, showing how far trait/phylogeny space diverges from chemistry ---
leaf.chem.ref <- leaf.aligned.df %>% filter(type == "Chemical feature") %>% select(ugs, Dim1_ref = Dim1, Dim2_ref = Dim2)
leaf.segments.df <- leaf.aligned.df %>% filter(type != "Chemical feature") %>% left_join(leaf.chem.ref, by = "ugs")

# --- 2f. Variance explained by each axis, for the panel's axis labels ---
leaf.pcoa.var.df <- bind_rows(pcoa_var(leaf.human.ord, "Human traits"),
                              pcoa_var(leaf.eco.ord, "Ecological traits"),
                              pcoa_var(leaf.phylo.ord, "Phylogeny"),
                              pcoa_var(leaf.chem.ord, "Chemical feature")) %>%
  filter(axis %in% c("PCoA1", "PCoA2"))

################################################################################
# ---- 3. FLOWER analysis ----------------------------------------------------
################################################################################

# --- 3a. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.chem.dist <- vegdist(flower.chem.mat %>% column_to_rownames("ugs"), method = "jaccard", binary = TRUE)

# --- 3b. PCoA (with a Lingoes correction for negative eigenvalues) on each of the four dissimilarity matrices ---
flower.human.ord <- wcmdscale(flower.human.dist.comdistnt, eig = TRUE, add = "lingoes")
flower.eco.ord <- wcmdscale(flower.eco.dist.comdistnt,   eig = TRUE, add = "lingoes")
flower.phylo.ord <- wcmdscale(flower.phylo.dist.comdistnt, eig = TRUE, add = "lingoes")
flower.chem.ord <- wcmdscale(flower.chem.dist,  eig = TRUE, add = "lingoes")

# only sites with valid coordinates in all four ordinations can be compared
flower.common.sites <- Reduce(intersect, list(rownames(flower.human.ord$points), rownames(flower.eco.ord$points),
                                              rownames(flower.phylo.ord$points), rownames(flower.chem.ord$points)))

# --- 3c. Procrustes rotation: rotate/rescale each trait or phylogeny ordination (first 2 axes) onto the chemical ordination, which acts as the common target space ---
flower.human.proc <- procrustes(flower.chem.ord$points[flower.common.sites, 1:2], flower.human.ord$points[flower.common.sites, 1:2], symmetric = FALSE)
flower.eco.proc <- procrustes(flower.chem.ord$points[flower.common.sites, 1:2], flower.eco.ord$points[flower.common.sites, 1:2],   symmetric = FALSE)
flower.phylo.proc <- procrustes(flower.chem.ord$points[flower.common.sites, 1:2], flower.phylo.ord$points[flower.common.sites, 1:2], symmetric = FALSE)

# collect the rotated coordinates from each of the 4 spaces into one data frame, plus the site's land-use category (ugs_type)
flower.aligned.df <- bind_rows(as.data.frame(flower.human.proc$Yrot) %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Human traits"),
                               as.data.frame(flower.eco.proc$Yrot)   %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Ecological traits"),
                               as.data.frame(flower.phylo.proc$Yrot) %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Phylogeny"),
                               as.data.frame(flower.human.proc$X)    %>% setNames(c("Dim1", "Dim2")) %>% rownames_to_column("ugs") %>% mutate(type = "Chemical feature")) %>%
  left_join(covariates.df %>% select(ugs, ugs_type), by = "ugs") %>%
  mutate(type = factor(type, levels = c("Human traits", "Ecological traits", "Phylogeny", "Chemical feature")))

# --- 3d. Fit a smoothed response surface of the human-influence index over the chemical-feature ordination space ---
flower.chem.points <- flower.aligned.df %>%
  filter(type == "Chemical feature") %>%
  left_join(covariates.df %>% select(ugs, human_index), by = "ugs") %>%
  filter(!is.na(human_index))

flower.surf.fit <- ordisurf(flower.chem.points[, c("Dim1", "Dim2")], flower.chem.points$human_index, plot = FALSE)

flower.surf.df <- expand.grid(Dim1 = flower.surf.fit$grid$x, Dim2 = flower.surf.fit$grid$y) %>%
  mutate(human_index = as.vector(flower.surf.fit$grid$z), type = "Chemical feature")

# --- 3e. Line segments from each alternative ordination's point to its true position in chemical space ---
flower.chem.ref <- flower.aligned.df %>% filter(type == "Chemical feature") %>% select(ugs, Dim1_ref = Dim1, Dim2_ref = Dim2)
flower.segments.df <- flower.aligned.df %>% filter(type != "Chemical feature") %>% left_join(flower.chem.ref, by = "ugs")

# --- 3f. Variance explained by each axis, for the panel's axis labels ---
flower.pcoa.var.df <- bind_rows(pcoa_var(flower.human.ord, "Human traits"),
                                pcoa_var(flower.eco.ord, "Ecological traits"),
                                pcoa_var(flower.phylo.ord, "Phylogeny"),
                                pcoa_var(flower.chem.ord, "Chemical feature")) %>%
  filter(axis %in% c("PCoA1", "PCoA2"))

################################################################################
# ---- 4. Build the two panels and save the figure --------------------------
################################################################################

# add a "tissue" column so facet_wrap() can draw the "Leaf"/"Flower" strip label
leaf.aligned.df$tissue  <- "Leaf"
leaf.surf.df$tissue <- "Leaf"
leaf.segments.df$tissue <- "Leaf"

flower.aligned.df$tissue  <- "Flower"
flower.surf.df$tissue <- "Flower"
flower.segments.df$tissue <- "Flower"

p1 <- ggplot() +
  metR::geom_contour2(data = leaf.surf.df, aes(x = Dim1, y = Dim2, z = human_index, group = type),
                      inherit.aes = FALSE, color = "grey25", linewidth = 0.3, bins = 8) +
  geom_segment(data = leaf.segments.df, aes(x = Dim1, y = Dim2, xend = Dim1_ref, yend = Dim2_ref, color = ugs_type),
               alpha = 0.45, linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = leaf.aligned.df %>% filter(type == "Chemical feature"), aes(x = Dim1, y = Dim2, color = ugs_type, shape = type), stroke = 0.5, size = 4) +
  geom_point(data = leaf.aligned.df %>% filter(type != "Chemical feature"), aes(x = Dim1, y = Dim2, color = ugs_type, fill = ugs_type, shape = type), stroke = 0, size = 2.5) +
  metR::geom_label_contour(data = leaf.surf.df, aes(x = Dim1, y = Dim2, z = human_index, group = type, label = sprintf("%.1f", ..level..)),
                           label.placer = metR::label_placer_fraction(0.3), label.size = 0, label.padding = unit(1, "pt"), bins = 8, size = 4, color = "black", skip = 1) +
  facet_wrap(~tissue) +
  scale_color_manual(values = site.colors) +
  scale_fill_manual(values = site.colors) +
  scale_shape_manual(values = type.shapes) +
  labs(x = axis_label(leaf.pcoa.var.df, "PCoA1", "Procrustes dimension 1"),
       y = axis_label(leaf.pcoa.var.df, "PCoA2", "Procrustes dimension 2"),
       color = NULL, fill = NULL, shape = NULL) +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14),
        legend.box = "vertical", legend.margin = ggplot2::margin(0, 0, 0, 0)) + 
  guides(shape = guide_legend(override.aes = list(fill = "black", color = "black")))

p2 <- ggplot() +
  metR::geom_contour2(data = flower.surf.df, aes(x = Dim1, y = Dim2, z = human_index, group = type),
                      inherit.aes = FALSE, color = "grey25", linewidth = 0.3, bins = 8) +
  geom_segment(data = flower.segments.df, aes(x = Dim1, y = Dim2, xend = Dim1_ref, yend = Dim2_ref, color = ugs_type),
               alpha = 0.45, linewidth = 0.4, show.legend = FALSE) +
  geom_point(data = flower.aligned.df %>% filter(type == "Chemical feature"), aes(x = Dim1, y = Dim2, color = ugs_type, shape = type), stroke = 0.5, size = 4) +
  geom_point(data = flower.aligned.df %>% filter(type != "Chemical feature"), aes(x = Dim1, y = Dim2, color = ugs_type, fill = ugs_type, shape = type), stroke = 0, size = 2.5) +
  metR::geom_label_contour(data = flower.surf.df, aes(x = Dim1, y = Dim2, z = human_index, group = type, label = sprintf("%.1f", ..level..)),
                           label.placer = metR::label_placer_fraction(0.3), label.size = 0, label.padding = unit(1, "pt"), bins = 8, size = 4, color = "black", skip = 1) +
  facet_wrap(~tissue) +
  scale_color_manual(values = site.colors) +
  scale_fill_manual(values = site.colors) +
  scale_shape_manual(values = type.shapes) +
  labs(x = axis_label(flower.pcoa.var.df, "PCoA1", "Procrustes dimension 1"),
       y = axis_label(flower.pcoa.var.df, "PCoA2", "Procrustes dimension 2"),
       color = NULL, shape = NULL) +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14),
        legend.box = "vertical", legend.margin = ggplot2::margin(0, 0, 0, 0))

ggarrange(p1, p2, common.legend = TRUE)

ggsave("figures/fig2.png", height = 7, width = 14)
ggsave("figures/fig2.pdf", height = 7, width = 14)
