library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggpattern)
library(mixOmics)

set.seed(123)

# Color palette for the five site land-use categories
site.colors <- c("Allotment" = "#FEC601",
                 "Garden" = "#072AC8",
                 "Park" = "#594F3B",
                 "Real Estate" = "#3DA5D9",
                 "Ruderal" = "#DD403A")

# Color palette distinguishing the two broad trait groups
trait_group_colors <- c("Human-associated" = "#B4563D",
                        "Ecological" = "#5B8DB8")


# ---- 1. Load input data --------------------------------------------------

# site x chemical-feature presence/absence matrices
# To recreate these datasets, use the prepare_metabolomics_df.R script. 
leaf.chem.mat <- read.csv("02_processed/leaf_site_mtb_mat.csv", check.names = FALSE)
flower.chem.mat <- read.csv("02_processed/flower_site_mtb_mat.csv", check.names = FALSE)

# metadata listing which plant species were actually profiled for each tissue
leaf.meta <- read.csv("01_raw/leaf_metadata.csv")
flower.meta <- read.csv("01_raw/flower_metadata.csv")

species.mat <- read.csv("01_raw/flora_mat.csv", check.names = FALSE)   # site x plant-species presence/absence
species.taxo <- read.csv("01_raw/flora_taxo.csv", check.names = FALSE)  # species-level taxonomy
species.traits.raw <- read.csv("01_raw/flora_traits.csv", check.names = FALSE) # raw species-level trait table

covariates.df <- read.csv("01_raw/covariates.csv") # per-site covariates

# ---- 2. Trait definitions and preprocessing (shared) -----------------------

# Traits considered "human-associated" (visible/experiential plant features) vs. "ecological" (habitat preference / life-history traits)
human_traits <- c("Growth.form", "Plant.height.vegetative..m.",
                  "Life.span", "Life.annual",
                  "Flower.color", "Flower.symmetry", "Flower.bilateral",
                  "Leaf.color", "Leaf.greenness", "Leaf.type", "Leaf.compoundness", "Leaf.elliptic", "Leaf.phenology", "Leaf.deciduous", "Fall.Conspicuous",
                  "is_alien")

eco_traits <- c("SLA..mm2.mg.", "Seed.dry.mass..mg.", "Leaf.water.content.per.leaf.dry.mass..g.W..g.DM.", "LDMC..g.g.",
                "Winter.hardiness", "Temperature", "Light", "Moisture", "Nutrients", "Humus", "Aeration",
                "CRS.life.strategies..competitive.", "CRS.life.strategies..stress_tolerant.", "CRS.life.strategies..ruderal.",
                "Flowering.duration")

all_traits <- c(human_traits, eco_traits)

# Human-readable labels for the derived predictor variables (CWM_* and PROP_*), used in plot legends/axes
trait_labels <- c("PROP_Growth.form__herb" = "Growth form: herb",
                  "PROP_Growth.form__liana" = "Growth form: liana",
                  "PROP_Growth.form__shrub" = "Growth form: shrub",
                  "PROP_Growth.form__tree" = "Growth form: tree",
                  "CWM_Plant.height.vegetative..m." = "Plant height (m)",
                  "CWM_Life.span" = "Life span",
                  "PROP_Flower.color__blue" = "Flower color: blue",
                  "PROP_Flower.color__green" = "Flower color: green",
                  "PROP_Flower.color__red" = "Flower color: red",
                  "PROP_Flower.color__white" = "Flower color: white",
                  "PROP_Flower.color__yellow" = "Flower color: yellow",
                  "CWM_Flower.bilateral" = "Bilateral symmetry",
                  "CWM_Leaf.greenness" = "Leaf morph: greenness",
                  "CWM_Leaf.compoundness" = "Leaf morph: compoundness",
                  "CWM_Leaf.elliptic" = "Leaf morph: ellipticness",
                  "CWM_Leaf.deciduous" = "Leaf morph: deciduousness",
                  "CWM_Fall.Conspicuous" = "Fall conspicuousness",
                  "CWM_SLA..mm2.mg." = "SLA",
                  "CWM_Seed.dry.mass..mg." = "Seed mass",
                  "CWM_Leaf.water.content.per.leaf.dry.mass..g.W..g.DM." = "Leaf water content",
                  "CWM_LDMC..g.g." = "LDMC",
                  "CWM_Winter.hardiness" = "Winter hardiness",
                  "CWM_Temperature" = "Temperature pref.",
                  "CWM_Light" = "Light pref.",
                  "CWM_Moisture" = "Moisture pref.",
                  "CWM_Nutrients" = "Nutrient pref.",
                  "CWM_Humus" = "Humus pref.",
                  "CWM_Aeration" = "Aeration pref.",
                  "CWM_CRS.life.strategies..competitive." = "CRS competitive",
                  "CWM_CRS.life.strategies..stress_tolerant." = "CRS stress-tolerant",
                  "CWM_CRS.life.strategies..ruderal." = "CRS ruderal",
                  "CWM_Flowering.duration" = "Flowering duration",
                  "CWM_is_alien" = "Non-native")

# Ordinal traits are stored in the raw data as ordered text categories; these maps convert each category to its rank (1 = lowest) so the trait can be averaged numerically
ordinal_maps <- list(Temperature = c("alpine and nival" = 1, "lower alpine, supra subalpine and upper subalpine" = 2, "subalpine" = 3,
                                     "lower subalpine and upper montane" = 4, "montane" = 5, "lower montane and upper colline" = 6,
                                     "colline" = 7, "warm colline" = 8, "very warm colline" = 9),
                     Light = c("deep shade" = 1, "shade" = 2, "semi-shade" = 3,
                               "well lit places" = 4, "full light" = 5),
                     Moisture = c("flooded, i.e. submerged" = 1, "wet" = 2, "very moist" = 3,
                                  "moderately moist" = 4, "moist" = 5, "dry" = 6,
                                  "moderately dry" = 7, "very dry" = 8),
                     Nutrients = c("very infertile" = 1, "infertile" = 2,
                                   "medium infertile to medium fertile" = 3,
                                   "fertile" = 4, "very fertile and over-rich" = 5),
                     Humus = c("little or no humus" = 1, "moderate humus content" = 2,
                               "high humus content" = 3),
                     Aeration = c("bad aeration" = 1, "moderate aeration" = 2, "good aeration" = 3))

# Cleans a raw species x trait table: converts character columns to factors, recodes ordinal traits to ranks, and collapses a few sparse/rare categories into a broader neighbouring category
process_traits <- function(df) {
  traits.df <- df %>%
    dplyr::mutate(across(where(is.character), as.factor)) %>%
    column_to_rownames("species")
  
  for (nm in intersect(names(ordinal_maps), colnames(traits.df))) {
    traits.df[[nm]] <- ordinal_maps[[nm]][as.character(traits.df[[nm]])]
    traits.df[[nm]] <- as.ordered(traits.df[[nm]])
  }
  
  traits.df <- traits.df %>%
    dplyr::mutate(Flower.color = case_when(Flower.color %in% ("orange") ~ "red", 
                                           Flower.color %in% c("pink/purple/violet") ~ "blue", 
                                           Flower.color %in% c("brown/beige/gray") ~ "white", 
                                           TRUE ~ Flower.color), 
                  Life.span = case_when(Life.span %in% c("annual", "biennial") ~ 1, 
                                        TRUE ~ 0),
                  Flower.bilateral = ifelse(Flower.bilateral == "bilateral", 1, 0),
                  Leaf.greenness = ifelse(Leaf.greenness == "green", 1, 0),
                  Leaf.elliptic = ifelse(Leaf.elliptic == "elliptic", 1, 0),
                  Leaf.compoundness = ifelse(Leaf.compoundness == "compound", 1, 0),
                  Leaf.deciduous = ifelse(Leaf.deciduous == "deciduous", 1, 0),
                  Fall.Conspicuous = ifelse(Fall.Conspicuous == "Fall Conspicuous", 1, 0)) %>%
    dplyr::mutate(Flower.color = as.factor(Flower.color))
  
  return(traits.df)
}

# Looks up, for every predictor trait, which trait group it belongs to (Human-associated / Ecological) and its display label
trait_group_lookup <- tibble(trait = all_traits,
                             trait_group = case_when(trait %in% c(human_traits) ~ "Human-associated", 
                                                     TRUE ~ "Ecological"),
                             trait_label = dplyr::recode(trait, !!!trait_labels, .default = trait))

# Builds community-weighted site-level trait predictors from a site x species community matrix and a species x trait table. 
# Numeric/ordered traits become site means (CWM_*)
# Factor traits become site-level proportions of species in each level (PROP_*), keeping only levels present at >= min_level_prop of sites.
make_site_trait_predictors <- function(comm, traits_df, min_level_prop = 0.03) {
  
  out <- list()
  meta <- list()
  
  for (tr in colnames(traits_df)) {
    x <- traits_df[[tr]]
    
    if (is.numeric(x) || is.ordered(x) || is.integer(x)) {
      xn <- as.numeric(x)
      vals <- apply(comm, 1, function(site_row) {
        spp <- names(site_row[site_row > 0])
        mean(xn[match(spp, rownames(traits_df))], na.rm = TRUE)
      })
      nm <- paste0("CWM_", make.names(tr))
      out[[nm]] <- vals
      meta[[nm]] <- tibble(variable = nm, trait = tr, level = NA_character_,
                           variable_type = "CWM numeric/ordinal")
    } else {
      xf <- as.factor(x)
      lvls <- levels(droplevels(xf))
      for (lv in lvls) {
        vals <- apply(comm, 1, function(site_row) {
          spp <- names(site_row[site_row > 0])
          trait_vals <- xf[match(spp, rownames(traits_df))]
          mean(trait_vals == lv, na.rm = TRUE)
        })
        if (mean(vals > 0, na.rm = TRUE) >= min_level_prop) {
          nm <- paste0("PROP_", make.names(tr), "__", make.names(lv))
          out[[nm]] <- vals
          meta[[nm]] <- tibble(variable = nm, trait = tr, level = as.character(lv),
                               variable_type = "Assemblage proportion")
        }
      }
    }
  }
  
  X <- as.data.frame(out)
  X$ugs <- rownames(comm)
  meta_df <- bind_rows(meta) %>%
    left_join(trait_group_lookup, by = "trait") %>%
    dplyr::mutate(display = ifelse(is.na(level), trait_label, paste0(trait_label, ": ", level)),
                  display = stringr::str_replace_all(display, "_", " "))
  
  list(X = X, meta = meta_df)
}

################################################################################
# ---- 3. LEAF analysis -------------------------------------------------------
################################################################################

# --- 3a. Restrict species-level data to the species profiled in leaf metabolomics ---
leaf.species.traits.raw <- species.traits.raw %>% filter(species %in% leaf.meta$Plants)
leaf.species.mat <- species.mat %>% dplyr::select(ugs, leaf.meta$Plants)

# --- 3b. Restrict to sites with chemistry, floristic, and covariate data ---
leaf.keep.sites <- Reduce(intersect, list(leaf.chem.mat$ugs, leaf.species.mat$ugs, covariates.df$ugs))
leaf.chem.mat <- leaf.chem.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))
leaf.species.mat <- leaf.species.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))

# --- 3c. Clean species traits and restrict to species present in the trait table ---
leaf.traits.df <- process_traits(leaf.species.traits.raw)

leaf.candidate_species <- intersect(species.taxo$species[species.taxo$class == "Magnoliopsida"],
                                    rownames(leaf.traits.df))
leaf.traits.df <- leaf.traits.df[leaf.candidate_species, , drop = FALSE]

leaf.species.comm <- leaf.species.mat %>%
  dplyr::select(all_of(leaf.candidate_species)) %>%
  as.matrix()
rownames(leaf.species.comm) <- leaf.species.mat$ugs

# --- 3d. Binarize the chemical feature matrix (presence/absence per site) ---
leaf.chem.pa <- leaf.chem.mat %>%
  column_to_rownames("ugs") %>%
  dplyr::mutate(across(everything(), ~ as.numeric(. > 0))) %>%
  as.matrix()

# --- 3e. Build community-weighted trait predictors per site (the PLS X matrix) ---
leaf.trait_obj <- make_site_trait_predictors(leaf.species.comm, leaf.traits.df %>% dplyr::select(any_of(all_traits)))
leaf.trait.mat <- leaf.trait_obj$X

# the set of predictor columns is fixed here (from the leaf data) and then also imposed on the flower predictors below, so both tissues share an identical set of predictor columns for comparability
trait_vars <- setdiff(colnames(leaf.trait.mat), "ugs")

leaf.trait.mat <- leaf.trait.mat %>%
  arrange(match(ugs, rownames(leaf.chem.pa))) %>%
  dplyr::select(all_of(trait_vars)) %>%
  as.matrix()

leaf.trait.mat.sc <- scale(leaf.trait.mat)

# Chemical features with no variation among sites are removed because they cannot contribute to model discrimination.
leaf.chem.pa.filt <- leaf.chem.pa[, apply(leaf.chem.pa, 2, sd, na.rm = TRUE) > 0]

# --- 3f. Fit the PLS regression (X = traits, Y = chemistry) ---
leaf.pls_fit <- mixOmics::pls(X = leaf.trait.mat.sc,
                              Y = leaf.chem.pa.filt,
                              ncomp = 2,
                              scale = FALSE, 
                              mode = "regression")

leaf.vip_mat <- mixOmics::vip(leaf.pls_fit)

# --- 3g. Trait loading "arrows" for the top 10 most important traits (by VIP, component 1) ---
leaf.arrow.df <- data.frame(variable = rownames(leaf.pls_fit$loadings$X),
                            PLS1 = leaf.pls_fit$loadings$X[, 1],
                            PLS2 = leaf.pls_fit$loadings$X[, 2]) %>%
  mutate(trait_label = dplyr::recode(variable, !!!trait_labels),
         label_text  = trait_label) %>%
  left_join(leaf.trait_obj$meta %>% distinct(variable, trait_group), by = "variable") %>%
  filter(variable %in% names(head(sort(leaf.vip_mat[, "comp1"], decreasing = TRUE), n = 10)))

# --- 3h. Site scores in PLS component space ---
leaf.scores.mat <- leaf.pls_fit$variates$Y[, 1:2]
colnames(leaf.scores.mat) <- c("PLS1", "PLS2")

leaf.scores.df <- as.data.frame(leaf.scores.mat) %>%
  rownames_to_column("ugs") %>%
  left_join(covariates.df %>% dplyr::select(ugs, ugs_type), by = "ugs")

# arrow scaling so arrows fill ~65% of the score cloud
leaf.score_range <- max(abs(c(leaf.scores.df$PLS1, leaf.scores.df$PLS2)), na.rm = TRUE) * 2
leaf.arrow_range <- max(sqrt(leaf.arrow.df$PLS1^2 + leaf.arrow.df$PLS2^2), na.rm = TRUE)
leaf.arrow_sc <- 0.65 * leaf.score_range / leaf.arrow_range

################################################################################
# ---- 4. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 4a. Restrict species-level data to the species profiled in floral VOC sampling ---
flower.species.traits.raw <- species.traits.raw %>% filter(species %in% flower.meta$Plants)
flower.species.mat <- species.mat %>% dplyr::select(ugs, flower.meta$Plants)

# --- 4b. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, flower.species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.species.mat <- flower.species.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))

# --- 4c. Clean species traits and restrict to species present in the trait table ---
flower.traits.df <- process_traits(flower.species.traits.raw)

flower.candidate_species <- intersect(species.taxo$species[species.taxo$class == "Magnoliopsida"],
                                      rownames(flower.traits.df))
flower.traits.df <- flower.traits.df[flower.candidate_species, , drop = FALSE]

flower.species.comm <- flower.species.mat %>%
  dplyr::select(all_of(flower.candidate_species)) %>%
  as.matrix()
rownames(flower.species.comm) <- flower.species.mat$ugs

# --- 4d. Binarize the chemical feature matrix (presence/absence per site) ---
flower.chem.pa <- flower.chem.mat %>%
  column_to_rownames("ugs") %>%
  dplyr::mutate(across(everything(), ~ as.numeric(. > 0))) %>%
  as.matrix()

# --- 4e. Build community-weighted trait predictors per site (the PLS X matrix) ---
flower.trait_obj <- make_site_trait_predictors(flower.species.comm, flower.traits.df %>% dplyr::select(any_of(all_traits)))
flower.trait.mat <- flower.trait_obj$X

# reuses trait_vars from the leaf analysis (see section 3e) so both tissues share an identical set of predictor columns
flower.trait.mat <- flower.trait.mat %>%
  arrange(match(ugs, rownames(flower.chem.pa))) %>%
  dplyr::select(all_of(trait_vars)) %>%
  as.matrix()

flower.trait.mat.sc <- scale(flower.trait.mat)

# Chemical features with no variation among sites are removed because they cannot contribute to model discrimination.
flower.chem.pa.filt <- flower.chem.pa[, apply(flower.chem.pa, 2, sd, na.rm = TRUE) > 0]

# --- 4f. Fit the PLS regression (X = traits, Y = chemistry) ---
flower.pls_fit <- mixOmics::pls(X = flower.trait.mat.sc,
                                Y = flower.chem.pa.filt,
                                ncomp = 2,
                                scale = FALSE, 
                                mode = "regression")

flower.vip_mat <- mixOmics::vip(flower.pls_fit)

# --- 4g. Trait loading "arrows" for the top 10 most important traits (by VIP, component 1) ---
flower.arrow.df <- data.frame(variable = rownames(flower.pls_fit$loadings$X),
                              PLS1 = flower.pls_fit$loadings$X[, 1],
                              PLS2 = flower.pls_fit$loadings$X[, 2]) %>%
  mutate(trait_label = dplyr::recode(variable, !!!trait_labels),
         label_text = trait_label) %>%
  left_join(flower.trait_obj$meta %>% distinct(variable, trait_group), by = "variable") %>%
  filter(variable %in% names(head(sort(flower.vip_mat[, "comp1"], decreasing = TRUE), n = 10)))

# --- 4h. Site scores in PLS component space ---
flower.scores.mat <- flower.pls_fit$variates$Y[, 1:2]
colnames(flower.scores.mat) <- c("PLS1", "PLS2")

flower.scores.df <- as.data.frame(flower.scores.mat) %>%
  rownames_to_column("ugs") %>%
  left_join(covariates.df %>% dplyr::select(ugs, ugs_type), by = "ugs")

flower.score_range <- max(abs(c(flower.scores.df$PLS1, flower.scores.df$PLS2)), na.rm = TRUE) * 1.5
flower.arrow_range <- max(sqrt(flower.arrow.df$PLS1^2 + flower.arrow.df$PLS2^2), na.rm = TRUE)
flower.arrow_sc <- 0.65 * flower.score_range / flower.arrow_range

# ---- 5. Combine tissues and build the VIP long-format data ----------------

arrow.df <- rbind(data.frame(leaf.arrow.df, arrow_sc = leaf.arrow_sc, tissue = "Leaf"),
                  data.frame(flower.arrow.df, arrow_sc = flower.arrow_sc, tissue = "Flower")) %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")))

scores.df <- rbind(data.frame(leaf.scores.df, tissue = "Leaf"),
                   data.frame(flower.scores.df, tissue = "Flower")) %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")))

# VIP scores for both PLS components, both tissues, in one long data frame
vip.df <- bind_rows(data.frame(variable = rownames(leaf.vip_mat),
                               comp1 = leaf.vip_mat[, 1],
                               comp2 = leaf.vip_mat[, 2],
                               tissue = "Leaf"),
                    data.frame(variable = rownames(flower.vip_mat),
                               comp1 = flower.vip_mat[, 1],
                               comp2 = flower.vip_mat[, 2],
                               tissue = "Flower"))  %>%
  left_join(bind_rows(leaf.trait_obj$meta,
                      flower.trait_obj$meta) %>%
              distinct(variable, trait, trait_group, display),
            by = "variable") %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")),
         trait_label = recode(variable, !!!trait_labels)) 

# order traits by their leaf VIP score (component 1), top 20, for the y-axis
trait_order <- vip.df %>%
  filter(tissue == "Leaf") %>%
  arrange(desc(comp1)) %>%
  slice_head(n = 20) %>%
  pull(trait_label) 

vip.df <- vip.df %>%
  mutate(trait_label = factor(trait_label, levels = rev(trait_order))) %>%
  pivot_longer(cols = c(comp1, comp2), names_to = "comp") %>%
  mutate(comp = factor(comp, levels = c("comp2", "comp1"))) %>%
  mutate(comp = recode(comp, "comp2" = "Component 2", "comp1" = "Component 1"),
         trait_group = factor(trait_group, levels = c("Human-associated", "Ecological"))) %>%
  na.omit()

# ---- 6. Build and save the figure -----------------------------------------

# Panel A: site scores + trait loading arrows, one facet per tissue
p1 <- ggplot() +
  stat_ellipse(data = scores.df,
               aes(x = PLS1, y = PLS2, color = ugs_type, fill = ugs_type),
               geom = "polygon", alpha = 0.08, linewidth = 0.5) +
  geom_point(data = scores.df,
             aes(x = PLS1, y = PLS2, color = ugs_type),
             size = 2.4, alpha = 0.9) +
  scale_color_manual(name = NULL, values = site.colors) +
  scale_fill_manual(name  = NULL, values = site.colors) +
  ggnewscale::new_scale_color() + # trait arrows use their own color legend, separate from ugs_type
  geom_segment(data = arrow.df,
               aes(x = 0, y = 0, xend  = PLS1 * arrow_sc, yend  = PLS2 * arrow_sc, color = trait_group),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.9,
               show.legend = c(color = TRUE, alpha = FALSE)) +
  ggrepel::geom_label_repel(data = arrow.df,
                            aes(x = PLS1 * arrow_sc, y = PLS2 * arrow_sc, label = stringr::str_wrap(label_text, 20), color = trait_group),
                            linewidth = 0, label.padding = 0.15,
                            box.padding = 0.3, point.padding = 0,
                            force = 2, show.legend = FALSE, size = 3.7) +
  facet_wrap(~ tissue, scales = "free") + 
  scale_color_manual(name = "Trait group", values = trait_group_colors) +
  scale_alpha_identity() +
  labs(x = "Component 1",
       y = "Component 2") +
  theme_pubr(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing = unit(0, "pt"),
        legend.box = "vertical", 
        legend.margin = ggplot2::margin(0, 0, 0, 0)) + 
  guides(color = "none")

# Panel B: VIP scores for the top 20 (leaf-ranked) traits, both components stacked
p2 <- ggplot(vip.df, aes(y = trait_label, x = value, fill = trait_group, color = trait_group, pattern = comp, pattern_color = trait_group, pattern_fill = trait_group, alpha = value < 1)) +
  geom_col_pattern(width = 0.7, position = "stack", linewidth = 0.5,
                   pattern_size = 0.35, pattern_angle = 45,
                   pattern_density = 0.1, pattern_spacing = 0.025,
                   pattern_key_scale_factor = .1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~ tissue) + 
  scale_color_manual(values = trait_group_colors) +
  scale_fill_manual(values = trait_group_colors) +
  scale_alpha_manual(values = c("FALSE" = 0.75, "TRUE" = 0.2)) +
  scale_pattern_manual(values = c("Component 1" = "none", "Component 2" = "stripe"), breaks = c("Component 1", "Component 2")) +
  scale_pattern_color_manual(values = trait_group_colors) +
  scale_pattern_fill_manual(values = trait_group_colors) +
  labs(y = NULL, 
       x = "VIP score", 
       fill = "Trait group:", color = "Trait group:", pattern = NULL) +
  theme_pubclean(base_size = 14) + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 14), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 14), 
        legend.justification = c(1.1, 0.5),
        legend.text = element_text(size = 14),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing = unit(20, "pt")) + 
  guides(alpha = "none", 
         pattern_color = "none",
         pattern_fill = "none", 
         fill = guide_legend(override.aes = list(pattern_alpha = 0)),
         pattern = guide_legend(override.aes = list(fill = "white", 
                                                    color = "black", 
                                                    pattern_fill = "black",
                                                    pattern_color = "black", 
                                                    pattern_density = 0.01, 
                                                    pattern_spacing = 0.08)))

ggarrange(p1, p2, 
          nrow = 2, 
          heights = c(0.9, 1),
          labels = c("A", "B"), font.label = list(size = 17))

ggsave("figures/fig4.png", height = 10, width = 10)
ggsave("figures/fig4.pdf", height = 10, width = 10)
