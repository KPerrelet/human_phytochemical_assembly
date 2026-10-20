library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(vegan)

set.seed(123)

# Color palette for the five site land-use categories
site.colors <- c("Allotment" = "#FEC601",
                 "Garden" = "#072AC8",
                 "Park" = "#594F3B",
                 "Real Estate" = "#3DA5D9",
                 "Ruderal" = "#DD403A")

# % of total variance explained by each positive PCoA axis of an ordination
pcoa_var <- function(ord) {
  eig <- ord$eig[ord$eig > 0]
  tibble(axis = paste0("PCoA", seq_along(eig)), variance_explained = eig / sum(eig) * 100)
}


# ---- 1. Load input data --------------------------------------------------

# metadata listing which plant species were actually profiled for each tissue
leaf.meta   <- read.csv("01_raw/leaf_metadata.csv")
flower.meta <- read.csv("01_raw/flower_metadata.csv")

species.mat   <- read.csv("01_raw/flora_mat.csv", check.names = FALSE) # site x plant-species presence/absence
covariates.df <- read.csv("01_raw/covariates.csv") # per-site covariates

leaf.species.mat   <- species.mat %>% select(ugs, leaf.meta$Plants) # restricted to species profiled for leaf metabolomics
flower.species.mat <- species.mat %>% select(ugs, flower.meta$Plants) # restricted to species profiled for flower VOCs

################################################################################
# ---- 2. LEAF analysis -----------------------------------------------------
################################################################################

# PCoA on a Jaccard dissimilarity of plant-community composition across leaf sites
leaf.species.pcoa <- cmdscale(vegdist(leaf.species.mat %>% select(-ugs), method = "jaccard", binary = TRUE), eig = TRUE)

leaf.pcoa.df <- data.frame(dim1 = leaf.species.pcoa$points[, 1],
                           dim2 = leaf.species.pcoa$points[, 2],
                           ugs = leaf.species.mat$ugs,
                           tissue = "Leaf") %>%
  left_join(covariates.df, by = "ugs")

p1 <- ggplot(leaf.pcoa.df, aes(x = dim1, y = dim2, color = ugs_type, fill = ugs_type)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  facet_wrap(~ tissue) +
  scale_color_manual(values = site.colors) +
  scale_fill_manual(values = site.colors) +
  labs(x = pcoa_var(leaf.species.pcoa) %>%
         filter(axis == "PCoA1") %>%
         with(paste0(axis, ": ", round(variance_explained, 2), "%")),
       y = pcoa_var(leaf.species.pcoa) %>%
         filter(axis == "PCoA2") %>%
         with(paste0(axis, ": ", round(variance_explained, 2), "%")),
       color = NULL, fill = NULL) +
  theme_pubr(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))


################################################################################
# ---- 3. FLOWER analysis ------------------------------------------------------
################################################################################

# PCoA on a Jaccard dissimilarity of plant-community composition across flower sites
flower.species.pcoa <- cmdscale(vegdist(flower.species.mat %>% select(-ugs), method = "jaccard", binary = TRUE), eig = TRUE)

flower.pcoa.df <- data.frame(dim1 = flower.species.pcoa$points[, 1],
                             dim2 = flower.species.pcoa$points[, 2],
                             ugs = flower.species.mat$ugs,
                             tissue = "Flower") %>%
  left_join(covariates.df, by = "ugs")

p2 <- ggplot(flower.pcoa.df, aes(x = dim1, y = dim2, color = ugs_type, fill = ugs_type)) +
  geom_point() +
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  facet_wrap(~ tissue) +
  scale_color_manual(values = site.colors) +
  scale_fill_manual(values = site.colors) +
  labs(x = pcoa_var(flower.species.pcoa) %>%
         filter(axis == "PCoA1") %>%
         with(paste0(axis, ": ", round(variance_explained, 2), "%")),
       y = pcoa_var(flower.species.pcoa) %>%
         filter(axis == "PCoA2") %>%
         with(paste0(axis, ": ", round(variance_explained, 2), "%")),
       color = NULL, fill = NULL) +
  theme_pubr(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))


# ---- 4. Assemble and save the figure ---------------------------------------

ggarrange(p1, p2, common.legend = TRUE)

ggsave("figures/figS9.png", height = 5, width = 10)
