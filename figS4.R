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

# ---- 1. Load input data --------------------------------------------------

# site x chemical-feature presence/absence matrices
# To recreate these datasets, use the prepare_metabolomics_df.R script. 
leaf.chem.mat <- read.csv("02_processed/leaf_site_mtb_mat.csv",   check.names = FALSE)
flower.chem.mat <- read.csv("02_processed/flower_site_mtb_mat.csv", check.names = FALSE)

# metadata listing which plant species were actually profiled for each tissue
leaf.meta <- read.csv("01_raw/leaf_metadata.csv")
flower.meta <- read.csv("01_raw/flower_metadata.csv")

species.mat <- read.csv("01_raw/flora_mat.csv", check.names = FALSE)   # site x plant-species presence/absence

covariates.df <- read.csv("01_raw/covariates.csv") # per-site covariates

leaf.species.mat   <- species.mat %>% select(ugs, leaf.meta$Plants)   # restricted to species profiled for leaf metabolomics
flower.species.mat <- species.mat %>% select(ugs, flower.meta$Plants) # restricted to species profiled for flower VOCs

################################################################################
# ---- 2. LEAF analysis -----------------------------------------------------
################################################################################

# site-level chemical feature richness and plant species richness
leaf.alpha.df <- data.frame(ugs = leaf.chem.mat$ugs,
                            feat_rich = specnumber(leaf.chem.mat %>% dplyr::select(-ugs)),
                            plant_rich = specnumber(leaf.species.mat %>% dplyr::select(-ugs)),
                            tissue = "Leaf")

################################################################################
# ---- 3. FLOWER analysis ------------------------------------------------------
################################################################################

flower.alpha.df <- data.frame(ugs = flower.chem.mat$ugs,
                              feat_rich = specnumber(flower.chem.mat %>% dplyr::select(-ugs)),
                              plant_rich = specnumber(flower.species.mat %>% dplyr::select(-ugs)),
                              tissue = "Flower")

# ---- 4. Combine tissues and compute pairwise significance groups ------------

alpha.df <- rbind(leaf.alpha.df, flower.alpha.df) %>%
  left_join(covariates.df, by = "ugs")

# for each tissue: pairwise Wilcoxon tests between all land-use types (BH-
# corrected), collapsed into compact-letter-display groups -- types sharing
# a letter are not significantly different from each other
cld_labels.df <- alpha.df %>%
  group_by(tissue) %>%
  group_modify(~ {
    pw <- pairwise.wilcox.test(.x$feat_rich, .x$ugs_type, p.adjust.method = "BH")
    
    p_mat <- pw$p.value
    p_vec <- as.vector(p_mat)
    
    rnames <- rownames(p_mat)
    cnames <- colnames(p_mat)
    names(p_vec) <- apply(expand.grid(rnames, cnames), 1, paste, collapse = "-")
    
    p_vec <- p_vec[!is.na(p_vec)]
    
    letters_out <- multcompView::multcompLetters(p_vec)$Letters
    
    tibble(ugs_type = names(letters_out),
           cld = letters_out,
           y_position = max(.x$feat_rich))
  }) %>%
  mutate(y_position = 1.05 * y_position) # place the letter just above the tallest point/box in each panel


# ---- 5. Build and save the figure -----------------------------------------

ggplot(alpha.df, aes(x = ugs_type, y = feat_rich, fill = ugs_type, color = ugs_type)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  geom_text(data = cld_labels.df,
            aes(x = ugs_type, y = y_position, label = cld),
            inherit.aes = FALSE,
            size = 4.5) +
  facet_wrap(~ factor(tissue, levels = c("Leaf", "Flower")), scales = "free") +
  scale_color_manual(values = site.colors) +
  scale_fill_manual(values = site.colors) +
  labs(x = NULL,
       y = "Feature richness") +
  theme_pubr(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("figures/figS4.png", height = 5, width = 8)
