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

# ---- 4. Combine tissues and fit per-group regressions -----------------------

alpha.df <- rbind(leaf.alpha.df, flower.alpha.df) %>%
  left_join(covariates.df, by = "ugs")

# axis ranges per tissue panel, used below to position the model-summary labels
tissue_ranges.df <- alpha.df %>%
  group_by(tissue) %>%
  dplyr::summarise(xmax = max(plant_rich, na.rm = TRUE),
                   ymin = min(feat_rich, na.rm = TRUE),
                   yrange = diff(range(feat_rich, na.rm = TRUE)))

# fit a negative-binomial regression of feature richness on (scaled) plant
# richness, separately for every tissue x land-use-type group, and build a
# "beta / pseudo-R2" label positioned in the corner of each facet
glm_labels.df <- alpha.df %>%
  dplyr::group_by(tissue, ugs_type) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    mod = purrr::map(data, ~ MASS::glm.nb(feat_rich ~ scale(plant_rich), data = .x)),
    beta = purrr::map_dbl(mod, ~ coef(.x)[["scale(plant_rich)"]]),
    pseudo_r2 = purrr::map_dbl(mod, ~ 1 - (.x$deviance / .x$null.deviance)),
    label = paste0("\u03b2 = ", round(beta, 2),
                   "; pseudo-R\u00b2 = ", round(pseudo_r2, 2))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(tissue, rev(ugs_type)) %>%
  dplyr::left_join(tissue_ranges.df, by = "tissue") %>%
  dplyr::group_by(tissue) %>%
  dplyr::mutate(x = xmax,
                y = ymin + dplyr::row_number() * 0.1 * yrange) %>% # stack labels vertically, one row per land-use type
  dplyr::ungroup()


# ---- 5. Build and save the figure -----------------------------------------

ggplot(alpha.df, aes(y = feat_rich, x = plant_rich, color = ugs_type, fill = ugs_type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", alpha = 0.2) +
  geom_text(data = glm_labels.df,
            aes(x = x, y = y, label = label, color = ugs_type),
            hjust = 1,
            vjust = 2,
            show.legend = FALSE) +
  facet_wrap(~ factor(tissue, levels = c("Leaf", "Flower")), scales = "free") +
  scale_color_manual(values = site.colors) +
  scale_fill_manual(values = site.colors) +
  labs(x = "Plant richness",
       y = "Feature richness",
       color = NULL, fill = NULL) +
  theme_pubr(base_size = 12) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))

ggsave("figures/figS3.png", height = 5, width = 10)
