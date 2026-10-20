library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(picante)
library(ape)
library(V.PhyloMaker2)
library(purrr)
library(mFD)
library(broom)
library(pscl)
library(ggnewscale)

set.seed(123)

# Standardizes a numeric vector (z-score), used as a predictor transform in
# the regression models below so coefficients are comparable across metrics
std <- function(x) as.numeric(scale(x))

# Converts a p-value into a conventional significance code
sig_label <- function(p) {
  dplyr::case_when(is.na(p) ~ "",
                   p < 0.001 ~ "***",
                   p < 0.01 ~ "**",
                   p < 0.05 ~ "*",
                   p < 0.1 ~ ".",
                   TRUE ~ "ns")
}

# Color palette for the five site land-use categories
site.colors <- c("Allotment" = "#FEC601",
                 "Garden" = "#072AC8",
                 "Park" = "#594F3B",
                 "Real Estate" = "#3DA5D9",
                 "Ruderal" = "#DD403A")

# Color palette distinguishing the two tissues
tissue.colors <- c("Leaf" = "#007C91", 
                   "Flower" = "#8F2D56")

# Computes site-level functional diversity (trait richness "fric" and trait dispersion "fdis") for a given subset of traits, using the mFD package.
# See Figure 5's script for the fuller step-by-step commentary on this function's internals.
compute_site_mfd <- function(trait_subset, comm, traits_df, label_prefix,
                             max_axes = 6, ind_vect = c("fric", "feve", "fdiv", "fdis")) {
  
  tr <- traits_df[intersect(rownames(traits_df), colnames(comm)),
                  intersect(trait_subset, colnames(traits_df)), drop = FALSE]
  
  tr_cat <- data.frame(trait_name = colnames(tr),
                       trait_type = dplyr::case_when(sapply(tr, is.ordered) ~ "O",
                                                     sapply(tr, is.factor)  ~ "N",
                                                     TRUE ~ "Q"))
  
  sp_ok <- intersect(rownames(tr), colnames(comm))
  comm_sub <- comm[, sp_ok, drop = FALSE]
  comm_sub <- comm_sub[rowSums(comm_sub) > 0, , drop = FALSE]
  tr <- tr[sp_ok, , drop = FALSE]
  
  funct_dist <- mFD::funct.dist(sp_tr = tr, tr_cat = tr_cat, metric = "gower",
                                scale_euclid = "scale_center", ordinal_var = "classic",
                                weight_type = "equal", stop_if_NA = FALSE)
  
  n_max <- min(max_axes, nrow(tr) - 1L)
  fspaces_q <- mFD::quality.fspaces(sp_dist = funct_dist, maxdim_pcoa = n_max,
                                    deviation_weighting = "absolute", fdendro = "average")
  mad_tbl <- fspaces_q$quality_fspaces
  pcoa_rows <- grep("^pcoa", rownames(mad_tbl))
  n_axes <- max(2L, as.integer(which.min(mad_tbl[pcoa_rows, "mad"])))
  
  message(label_prefix, ": ", n_axes, " PCoA axes (MAD = ",
          round(min(mad_tbl[pcoa_rows, "mad"]), 4), ")")
  
  sp_faxes <- fspaces_q$details_fspaces$sp_pc_coord[, seq_len(n_axes), drop = FALSE]
  
  mFD::alpha.fd.multidim(sp_faxes_coord = sp_faxes, asb_sp_w = comm_sub, ind_vect = ind_vect,
                         scaling = TRUE, check_input = TRUE,
                         details_returned = FALSE)$functional_diversity_indices %>%
    tibble::rownames_to_column("ugs") %>%
    dplyr::rename_with(~ paste0(label_prefix, "_", .x), .cols = -ugs)
}

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

ecological_traits <- c("SLA..mm2.mg.", "Seed.dry.mass..mg.", "Leaf.water.content.per.leaf.dry.mass..g.W..g.DM.", "LDMC..g.g.",
                       "Winter.hardiness", "Temperature", "Light", "Moisture", "Nutrients", "Humus", "Aeration",
                       "CRS.life.strategies..competitive.", "CRS.life.strategies..stress_tolerant.", "CRS.life.strategies..ruderal.",
                       "Flowering.duration")

# Ordinal traits are stored in the raw data as ordered text categories; these maps convert each category to its rank (1 = lowest) so the trait can be used in the Gower dissimilarity as ordered/numeric rather than nominal
ordinal_maps <- list(Temperature = c("alpine and nival" = 1, "lower alpine, supra subalpine and upper subalpine" = 2, "subalpine" = 3,
                                     "lower subalpine and upper montane" = 4, "montane" = 5, "lower montane and upper colline" = 6,
                                     "colline" = 7, "warm colline" = 8, "very warm colline" = 9),
                     Light = c("deep shade" = 1, "shade" = 2, "semi-shade" = 3, "well lit places" = 4, "full light" = 5),
                     Moisture = c("flooded, i.e. submerged" = 1, "wet" = 2, "very moist" = 3, "moderately moist" = 4,
                                  "moist" = 5, "dry" = 6, "moderately dry" = 7, "very dry" = 8),
                     Nutrients = c("very infertile" = 1, "infertile" = 2, "medium infertile to medium fertile" = 3,
                                   "fertile" = 4, "very fertile and over-rich" = 5),
                     Humus = c("little or no humus" = 1, "moderate humus content" = 2, "high humus content" = 3),
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


################################################################################
# ---- 3. LEAF analysis -------------------------------------------------------
################################################################################

# --- 3a. Restrict species-level data to the species profiled in leaf metabolomics ---
leaf.species.traits.raw <- species.traits.raw %>% filter(species %in% leaf.meta$Plants)
leaf.species.mat <- species.mat %>% select(ugs, leaf.meta$Plants)
leaf.species.taxo <- species.taxo %>% filter(species %in% leaf.meta$Plants)

# --- 3b. Restrict to sites with chemistry, floristic, and covariate data ---
leaf.keep.sites <- Reduce(intersect, list(leaf.chem.mat$ugs, leaf.species.mat$ugs, covariates.df$ugs))
leaf.chem.mat <- leaf.chem.mat    %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))
leaf.species.mat <- leaf.species.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))
leaf.covariates.df <- covariates.df %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))

# --- 3c. Clean species traits and build the site x species community matrix ---
leaf.traits.df <- process_traits(leaf.species.traits.raw)

leaf.candidate_species <- intersect(leaf.species.taxo$species[leaf.species.taxo$class == "Magnoliopsida"],
                                    intersect(colnames(leaf.species.mat), rownames(leaf.traits.df)))
leaf.traits.df <- leaf.traits.df[leaf.candidate_species, , drop = FALSE]

leaf.species.comm <- leaf.species.mat %>%
  dplyr::select(all_of(leaf.candidate_species)) %>%
  dplyr::mutate(across(everything(), ~ as.numeric(. > 0))) %>%
  as.matrix()
rownames(leaf.species.comm) <- leaf.species.mat$ugs

# --- 3d. Build a phylogeny and compute site-level phylogenetic community structure (MPD, MNTD; picante) ---
leaf.tree.scenarios <- V.PhyloMaker2::phylo.maker(sp.list = leaf.species.taxo %>% dplyr::select(species, genus, family),
                                                  scenarios = c("S1", "S2", "S3"))

leaf.tree <- leaf.tree.scenarios$scenario.3
leaf.tree$tip.label <- gsub("_", " ", leaf.tree$tip.label)
leaf.tree$tip.label <- gsub("  ", " _", leaf.tree$tip.label)

leaf.tree.species <- intersect(colnames(leaf.species.comm), leaf.tree$tip.label)
leaf.tree <- ape::keep.tip(leaf.tree, leaf.tree.species)
leaf.species.comm <- leaf.species.comm[, leaf.tree$tip.label, drop = FALSE]
leaf.traits.df <- leaf.traits.df[leaf.tree$tip.label, , drop = FALSE]

leaf.coph.mat <- ape::cophenetic.phylo(leaf.tree)
leaf.coph.mat <- leaf.coph.mat[colnames(leaf.species.comm), colnames(leaf.species.comm)]

leaf.mpd.obs <- picante::mpd(leaf.species.comm, leaf.coph.mat, abundance.weighted = FALSE)
leaf.mntd.obs <- picante::mntd(leaf.species.comm, leaf.coph.mat, abundance.weighted = FALSE)

leaf.phylo.alpha.df <- tibble(ugs = rownames(leaf.species.comm), MPD = leaf.mpd.obs, MNTD = leaf.mntd.obs)

# --- 3e. Community-level functional diversity for each trait group ---
leaf.human_fd <- compute_site_mfd(human_traits, leaf.species.comm, leaf.traits.df, "human_trait")
leaf.eco_fd   <- compute_site_mfd(ecological_traits, leaf.species.comm, leaf.traits.df, "eco_trait")

leaf.human_fric <- leaf.human_fd %>% dplyr::select(ugs, human_trait_fric)
leaf.eco_fric <- leaf.eco_fd   %>% dplyr::select(ugs, eco_trait_fric)
leaf.human_fdis <- leaf.human_fd %>% dplyr::select(ugs, human_trait_fdis)
leaf.eco_fdis <- leaf.eco_fd   %>% dplyr::select(ugs, eco_trait_fdis)

# --- 3f. Site-level alpha-diversity table for every leaf site ---
leaf.alpha.df <- tibble(ugs = leaf.chem.mat$ugs,
                        mtb_rich = vegan::specnumber(leaf.chem.mat %>% select(-ugs)),
                        plant_rich = vegan::specnumber(leaf.species.comm)) %>%
  left_join(leaf.covariates.df, by = "ugs") %>%
  left_join(leaf.human_fdis, by = "ugs") %>%
  left_join(leaf.eco_fdis, by = "ugs") %>%
  left_join(leaf.human_fric, by = "ugs") %>%
  left_join(leaf.eco_fric, by = "ugs") %>%
  left_join(leaf.phylo.alpha.df, by = "ugs")

# --- 3g. Reshape to long format and regress each metric on the
# human-influence index ---
leaf.alpha.long <- leaf.alpha.df %>%
  dplyr::select(ugs, ugs_type, human_index,
                plant_rich, MPD, MNTD,
                human_trait_fdis, eco_trait_fdis, human_trait_fric, eco_trait_fric,
                mtb_rich) %>%
  pivot_longer(cols = -c(ugs, ugs_type, human_index), names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("plant_rich", "MPD", "MNTD", "mtb_rich",
                                                   "human_trait_fric", "human_trait_fdis", "eco_trait_fric", "eco_trait_fdis"))) %>%
  dplyr::mutate(metric = recode(metric,
                                plant_rich = "Plant\nrichness",
                                mtb_rich = "Feature\nrichness",
                                human_trait_fric = "Human trait\nrichness",
                                human_trait_fdis = "Human trait dispersion",
                                eco_trait_fric = "Ecological trait\nrichness",
                                eco_trait_fdis = "Ecological trait dispersion"))

# fit a model per metric: negative-binomial GLM for the two count metrics
# (plant and feature richness), linear model for the rest
leaf.alpha.reg <- leaf.alpha.long %>%
  group_by(metric) %>%
  dplyr::summarise(model = list(if (unique(metric) %in% c("Plant\nrichness", "Feature\nrichness")) {
    MASS::glm.nb(value ~ std(human_index), data = cur_data()) } else { lm(value ~ std(human_index), data = cur_data())
    }), .groups = "drop") %>%
  dplyr::mutate(beta = map_dbl(model, ~coef(.x)["std(human_index)"]),
                p.value = map_dbl(model, ~broom::tidy(.x) %>% filter(term == "std(human_index)") %>% pull(p.value)),
                r2 = map2_dbl(model, metric, ~{
                  if (.y %in% c("Plant\nrichness", "Feature\nrichness")) { pscl::pR2(.x)["McFadden"] } else { summary(.x)$r.squared }
                })) %>%
  dplyr::select(-model) %>%
  dplyr::mutate(signif = p.value < 0.05) %>%
  dplyr::mutate(sig = sig_label(!signif),
                sig_star = case_when(p.value < 0.001 ~ "< 0.001", p.value < 0.01 ~ "< 0.01",
                                     p.value < 0.05 ~ "< 0.05", TRUE ~ "> 0.05"),
                sig_label = ifelse(signif == TRUE, "Significant", "Non-significant"))

leaf.alpha.long <- leaf.alpha.long %>% left_join(leaf.alpha.reg, by = "metric")

################################################################################
# ---- 4. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 4a. Restrict species-level data to the species profiled in floral VOC sampling ---
flower.species.traits.raw <- species.traits.raw %>% filter(species %in% flower.meta$Plants)
flower.species.mat <- species.mat %>% select(ugs, flower.meta$Plants)
flower.species.taxo <- species.taxo %>% filter(species %in% flower.meta$Plants)

# --- 4b. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, flower.species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.species.mat <- flower.species.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.covariates.df <- covariates.df %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))

# --- 4c. Clean species traits and build the site x species community matrix ---
flower.traits.df <- process_traits(flower.species.traits.raw)

flower.candidate_species <- intersect(flower.species.taxo$species[flower.species.taxo$class == "Magnoliopsida"],
                                      intersect(colnames(flower.species.mat), rownames(flower.traits.df)))
flower.traits.df <- flower.traits.df[flower.candidate_species, , drop = FALSE]

flower.species.comm <- flower.species.mat %>%
  dplyr::select(all_of(flower.candidate_species)) %>%
  dplyr::mutate(across(everything(), ~ as.numeric(. > 0))) %>%
  as.matrix()
rownames(flower.species.comm) <- flower.species.mat$ugs

# --- 4d. Build a phylogeny and compute MPD / MNTD ---
flower.tree.scenarios <- V.PhyloMaker2::phylo.maker(sp.list = flower.species.taxo %>% dplyr::select(species, genus, family),
                                                    scenarios = c("S1", "S2", "S3"))

flower.tree <- flower.tree.scenarios$scenario.3
flower.tree$tip.label <- gsub("_", " ", flower.tree$tip.label)
flower.tree$tip.label <- gsub("  ", " _", flower.tree$tip.label)

flower.tree.species <- intersect(colnames(flower.species.comm), flower.tree$tip.label)
flower.tree <- ape::keep.tip(flower.tree, flower.tree.species)
flower.species.comm <- flower.species.comm[, flower.tree$tip.label, drop = FALSE]
flower.traits.df <- flower.traits.df[flower.tree$tip.label, , drop = FALSE]

flower.coph.mat <- ape::cophenetic.phylo(flower.tree)
flower.coph.mat <- flower.coph.mat[colnames(flower.species.comm), colnames(flower.species.comm)]

flower.mpd.obs <- picante::mpd(flower.species.comm, flower.coph.mat, abundance.weighted = FALSE)
flower.mntd.obs <- picante::mntd(flower.species.comm, flower.coph.mat, abundance.weighted = FALSE)

flower.phylo.alpha.df <- tibble(ugs = rownames(flower.species.comm), MPD = flower.mpd.obs, MNTD = flower.mntd.obs)

# --- 4e. Community-level functional diversity for each trait group ---
flower.human_fd <- compute_site_mfd(human_traits, flower.species.comm, flower.traits.df, "human_trait")
flower.eco_fd <- compute_site_mfd(ecological_traits, flower.species.comm, flower.traits.df, "eco_trait")

flower.human_fric <- flower.human_fd %>% dplyr::select(ugs, human_trait_fric)
flower.eco_fric <- flower.eco_fd %>% dplyr::select(ugs, eco_trait_fric)
flower.human_fdis <- flower.human_fd %>% dplyr::select(ugs, human_trait_fdis)
flower.eco_fdis <- flower.eco_fd %>% dplyr::select(ugs, eco_trait_fdis)

# --- 4f. Site-level alpha-diversity table for every flower site ---
flower.alpha.df <- tibble(ugs = flower.chem.mat$ugs,
                          voc_rich = vegan::specnumber(flower.chem.mat %>% select(-ugs)),
                          plant_rich = vegan::specnumber(flower.species.comm)) %>%
  left_join(flower.covariates.df, by = "ugs") %>%
  left_join(flower.human_fdis, by = "ugs") %>%
  left_join(flower.eco_fdis, by = "ugs") %>%
  left_join(flower.human_fric, by = "ugs") %>%
  left_join(flower.eco_fric, by = "ugs") %>%
  left_join(flower.phylo.alpha.df, by = "ugs")

# --- 4g. Reshape to long format and regress each metric on human influence ---
flower.alpha.long <- flower.alpha.df %>%
  dplyr::select(ugs, ugs_type, human_index,
                plant_rich, MPD, MNTD,
                human_trait_fdis, eco_trait_fdis, human_trait_fric, eco_trait_fric,
                voc_rich) %>%
  pivot_longer(cols = -c(ugs, ugs_type, human_index), names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("plant_rich", "MPD", "MNTD", "voc_rich",
                                                   "human_trait_fric", "human_trait_fdis", "eco_trait_fric", "eco_trait_fdis"))) %>%
  dplyr::mutate(metric = recode(metric,
                                plant_rich = "Plant\nrichness",
                                voc_rich = "Feature\nrichness",
                                human_trait_fric = "Human trait\nrichness",
                                human_trait_fdis = "Human trait dispersion",
                                eco_trait_fric = "Ecological trait\nrichness",
                                eco_trait_fdis = "Ecological trait dispersion"))

flower.alpha.reg <- flower.alpha.long %>%
  group_by(metric) %>%
  dplyr::summarise(model = list(if (unique(metric) %in% c("Plant\nrichness", "Feature\nrichness")) {
    MASS::glm.nb(value ~ std(human_index), data = cur_data())
    } else { 
      lm(value ~ std(human_index), data = cur_data())
    }), .groups = "drop") %>%
  dplyr::mutate(beta = map_dbl(model, ~coef(.x)["std(human_index)"]),
                p.value = map_dbl(model, ~broom::tidy(.x) %>% filter(term == "std(human_index)") %>% pull(p.value)),
                r2 = map2_dbl(model, metric, ~{
                  if (.y %in% c("Plant\nrichness", "Feature\nrichness")) { pscl::pR2(.x)["McFadden"] } else { summary(.x)$r.squared }
                })) %>%
  dplyr::select(-model) %>%
  dplyr::mutate(signif = p.value < 0.05) %>%
  dplyr::mutate(sig = sig_label(!signif),
                sig_star = case_when(p.value < 0.001 ~ "< 0.001", p.value < 0.01 ~ "< 0.01",
                                     p.value < 0.05 ~ "< 0.05", TRUE ~ "> 0.05"),
                sig_label = ifelse(signif == TRUE, "Significant", "Non-significant"))

flower.alpha.long <- flower.alpha.long %>% left_join(flower.alpha.reg, by = "metric")

# ---- 5. Combine tissues, scale metrics, and build the figure ---------------

alpha.long <- rbind(data.frame(leaf.alpha.long %>%
                                 dplyr::filter(metric %in% c("MPD", "MNTD", "Feature\nrichness", "Plant\nrichness", "Human trait\nrichness", "Ecological trait\nrichness")) %>%
                                 dplyr::group_by(metric) %>%
                                 dplyr::mutate(value = bruceR::scaler(value)),
                               tissue = "Leaf"),
                    data.frame(flower.alpha.long %>%
                                 dplyr::filter(metric %in% c("MPD", "MNTD", "Feature\nrichness", "Plant\nrichness", "Human trait\nrichness", "Ecological trait\nrichness")) %>%
                                 dplyr::group_by(metric) %>%
                                 dplyr::mutate(value = bruceR::scaler(value)),
                               tissue = "Flower")) %>%
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")),
         metric = factor(metric, levels = c("Feature\nrichness", "Plant\nrichness", "Human trait\nrichness", "Ecological trait\nrichness", "MPD", "MNTD")))

alpha.reg <- rbind(data.frame(leaf.alpha.reg   %>% dplyr::filter(metric %in% c("MPD", "MNTD", "Feature\nrichness", "Plant\nrichness", "Human trait\nrichness", "Ecological trait\nrichness")), tissue = "Flower"),
                   data.frame(flower.alpha.reg %>% dplyr::filter(metric %in% c("MPD", "MNTD", "Feature\nrichness", "Plant\nrichness", "Human trait\nrichness", "Ecological trait\nrichness")), tissue = "Leaf"))

ggplot(alpha.long, aes(x = human_index, y = value, linetype = sig_label, alpha = sig_label)) +
  geom_point(aes(color = ugs_type), size = 1.4, alpha = 0.25) +
  scale_color_manual(values = site.colors, name = NULL) +
  ggnewscale::new_scale_color() + # the smooths below use their own color legend (by tissue), separate from ugs_type
  geom_smooth(aes(color = tissue, fill = tissue), method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
  geom_text(data = alpha.reg,
            aes(x = -Inf, y = c(1, 1, 1, 1, 1, 1, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9), label = sprintf("R\u00B2 = %.2f", r2), color = tissue),
            inherit.aes = FALSE, hjust = -0.1, vjust = 0, size = 5, show.legend = FALSE) +
  facet_wrap(~ metric, ncol = 6, scales = "free_y") +
  scale_color_manual(values = tissue.colors, name = NULL) +
  scale_fill_manual(values = tissue.colors, name = NULL) +
  scale_linetype_manual(values = c("Significant" = "solid", "Non-significant" = "dotted")) +
  scale_alpha_manual(values = c("Significant" = 0.4, "Non-significant" = 0.1)) +
  labs(x = "Human-influence index", y = NULL, alpha = NULL, linetype = NULL) +
  theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0),
        axis.text = element_text(size = 14), legend.text = element_text(size = 14))

ggsave("figures/figS7.png", height = 5, width = 17)
