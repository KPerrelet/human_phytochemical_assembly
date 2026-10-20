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
leaf.chem.mat <- read.csv("02_processed/leaf_site_mtb_mat.csv",   check.names = FALSE)
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
                  "CWM_Leaf.greenness" = "Leaf greenness",
                  "CWM_Leaf.compoundness" = "Leaf compoundness",
                  "CWM_Leaf.elliptic" = "Leaf ellipticness",
                  "CWM_Leaf.deciduous" = "Leaf deciduousness",
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
leaf.chem.mat    <- leaf.chem.mat    %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))
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

################################################################################
# ---- 4. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 4a. Restrict species-level data to the species profiled in floral VOC sampling ---
flower.species.traits.raw <- species.traits.raw %>% filter(species %in% flower.meta$Plants)
flower.species.mat <- species.mat %>% dplyr::select(ugs, flower.meta$Plants)

# --- 4b. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, flower.species.mat$ugs, covariates.df$ugs))
flower.chem.mat    <- flower.chem.mat    %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
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

# ---- 5. Fit PLS models and cross-validate component count -----------------

# fit with more components than we intend to use (5, vs. the 2 used in Figure 4) so that perf() below can show us how much each additional component actually helps predict out-of-sample data
leaf.pls_fit <- mixOmics::pls(X = leaf.trait.mat.sc,
                              Y = leaf.chem.pa.filt,
                              ncomp = 5,
                              scale = FALSE, 
                              mode = "regression")

flower.pls_fit <- mixOmics::pls(X = flower.trait.mat.sc,
                                Y = flower.chem.pa.filt,
                                ncomp  = 5,
                                scale = FALSE, 
                                mode = "regression")

# repeated 10-fold cross-validation, computing Q^2 (total) per component --
leaf.perf_pls <- perf(object = leaf.pls_fit,
                      validation = "Mfold",
                      folds = 10,
                      nrepeat = 100,
                      BPPARAM = BiocParallel::SnowParam(workers = 10), 
                      progressBar = TRUE)

flower.perf_pls <- perf(object = flower.pls_fit,
                        validation = "Mfold",
                        folds = 10,
                        nrepeat = 100,
                        BPPARAM = BiocParallel::SnowParam(workers = 10), 
                        progressBar = TRUE)

# ---- 6. Build and save the figure -----------------------------------------

q2.df <- rbind(data.frame(leaf.perf_pls$measures$Q2.total$summary, tissue = "Leaf"), 
               data.frame(flower.perf_pls$measures$Q2.total$summary, tissue = "Flower")) %>%
  dplyr::rename(comp = comp, Q2_mean = mean, Q2_sd = sd) %>%
  dplyr::mutate(comp = factor(comp), 
                tissue = factor(tissue, levels = c("Leaf", "Flower")))

ggplot(q2.df, aes(x = comp, y = Q2_mean, fill = Q2_mean > 0)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_errorbar(aes(ymin = Q2_mean - Q2_sd, ymax = Q2_mean + Q2_sd), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = 0.0975, linetype = "dotted", color = "#072AC8", linewidth = 0.8) +
  annotate("text", x = Inf, y = 0.0975, hjust = 1.1, vjust = -0.4,
           label = "Q\u00b2 = 0.0975", color = "black", size = 5) +
  facet_wrap(~ tissue) + 
  scale_fill_manual(values = c("TRUE" = "#3DA5D9", "FALSE" = "#DD403A"), guide = "none") +
  labs(x = "Component", y = expression(Q^2~(total))) +
  theme_pubr(base_size = 14) + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 14), 
        axis.text.y = element_text(size = 12))

ggsave("figures/figS10.png", height = 7, width = 10)
