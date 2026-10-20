library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(picante)
library(ape)
library(V.PhyloMaker2)
library(purrr)
library(piecewiseSEM)
library(semEff)
library(mFD)
library(vegan)

set.seed(123) # for reproducibility of the SEM bootstrap (bootEff) below

# Standardizes a numeric vector (z-score), used as a predictor transform in the regression/SEM models below so coefficients are comparable across metrics
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

# Computes site-level functional diversity (func. richness "fric" and trait dispersion "fdis") for a given subset of traits, using the mFD package:
#   1. build a Gower trait-dissimilarity matrix between species
#   2. choose the number of PCoA axes that best represents that dissimilarity (lowest mean absolute deviation, "mad")
#   3. compute functional-diversity indices per site from those axes
# comm: site x species community matrix. traits_df: species x trait table.
# label_prefix: prefixes the output column names (e.g. "human_trait_fric").
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

  sp_faxes <- fspaces_q$details_fspaces$sp_pc_coord[, seq_len(n_axes), drop = FALSE]
  
  mFD::alpha.fd.multidim(sp_faxes_coord = sp_faxes, asb_sp_w = comm_sub, ind_vect = ind_vect,
                         scaling = TRUE, check_input = TRUE,
                         details_returned = FALSE)$functional_diversity_indices %>%
    tibble::rownames_to_column("ugs") %>%
    dplyr::rename_with(~ paste0(label_prefix, "_", .x), .cols = -ugs)
}

# Calculate deviance explained (alternative to R2)
deviance_explained <- function(model) {
  1 - model$deviance / model$null.deviance
}

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

ecological_traits <- c("SLA..mm2.mg.", "Seed.dry.mass..mg.", "Leaf.water.content.per.leaf.dry.mass..g.W..g.DM.", "LDMC..g.g.",
                       "Winter.hardiness", "Temperature", "Light", "Moisture", "Nutrients", "Humus", "Aeration",
                       "CRS.life.strategies..competitive.", "CRS.life.strategies..stress_tolerant.", "CRS.life.strategies..ruderal.",
                       "Flowering.duration")

all_traits <- c(human_traits, ecological_traits)

# Ordinal traits are stored in the raw data as ordered text categories; these maps convert each category to its rank (1 = lowest) so the trait can be treated as ordered/numeric
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

# --- 3b. Restrict to sites with chemistry, floristic, and covariate data, and
# align row order across the chemistry and floristic matrices (needed below: the site-level alpha-diversity table in 3f is built by row position, not by join, so the two matrices must stay in the same site order) ---
leaf.keep.sites <- Reduce(intersect, list(leaf.chem.mat$ugs, leaf.species.mat$ugs, covariates.df$ugs))
leaf.chem.mat <- leaf.chem.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))
leaf.species.mat <- leaf.species.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))

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

# --- 3d. Build a phylogeny for the candidate species (V.PhyloMaker2), then compute site-level phylogenetic community structure (MPD, MNTD; picante) ---
leaf.tree.scenarios <- V.PhyloMaker2::phylo.maker(sp.list = leaf.species.taxo %>% dplyr::select(species, genus, family),
                                                  scenarios = c("S1", "S2", "S3"))

# scenario 3: species missing from the mega-tree are bound in as polytomies at the genus (or family) level 
leaf.tree <- leaf.tree.scenarios$scenario.3
leaf.tree$tip.label <- gsub("_", " ", leaf.tree$tip.label)
leaf.tree$tip.label <- gsub("  ", " _", leaf.tree$tip.label)

leaf.tree.species  <- intersect(colnames(leaf.species.comm), leaf.tree$tip.label)
leaf.tree <- ape::keep.tip(leaf.tree, leaf.tree.species)
leaf.species.comm <- leaf.species.comm[, leaf.tree$tip.label, drop = FALSE]
leaf.traits.df <- leaf.traits.df[leaf.tree$tip.label, , drop = FALSE]

leaf.coph.mat <- ape::cophenetic.phylo(leaf.tree)
leaf.coph.mat <- leaf.coph.mat[colnames(leaf.species.comm), colnames(leaf.species.comm)]

leaf.mpd.obs <- picante::mpd(leaf.species.comm, leaf.coph.mat, abundance.weighted = FALSE)
leaf.mntd.obs <- picante::mntd(leaf.species.comm, leaf.coph.mat, abundance.weighted = FALSE)

leaf.phylo.alpha.df <- tibble(ugs = rownames(leaf.species.comm), MPD = leaf.mpd.obs, MNTD = leaf.mntd.obs)

# --- 3e. Community-level functional diversity for each trait group (mFD package: func. richness "fric" and trait dispersion "fdis") ---
leaf.human_fd <- compute_site_mfd(human_traits, leaf.species.comm, leaf.traits.df, "human_trait")
leaf.eco_fd <- compute_site_mfd(ecological_traits, leaf.species.comm, leaf.traits.df, "ecological_trait")

leaf.human_fric <- leaf.human_fd %>% dplyr::select(ugs, human_trait_fric)
leaf.eco_fric <- leaf.eco_fd %>% dplyr::select(ugs, ecological_trait_fric)
leaf.human_fdis <- leaf.human_fd %>% dplyr::select(ugs, human_trait_fdis)
leaf.eco_fdis <- leaf.eco_fd %>% dplyr::select(ugs, ecological_trait_fdis)

# --- 3f. Site-level alpha-diversity table: chemical, floristic, functional, and phylogenetic diversity, plus covariates, for every leaf site ---
leaf.alpha.df <- tibble(ugs = leaf.chem.mat$ugs,
                        mtb_rich = vegan::specnumber(leaf.chem.mat %>% select(-ugs)),
                        plant_rich = vegan::specnumber(leaf.species.comm)) %>%
  left_join(covariates.df, by = "ugs") %>%
  left_join(leaf.human_fdis, by = "ugs") %>%
  left_join(leaf.eco_fdis, by = "ugs") %>%
  left_join(leaf.human_fric, by = "ugs") %>%
  left_join(leaf.eco_fric, by = "ugs") %>%
  left_join(leaf.phylo.alpha.df, by = "ugs")

# --- 3g. Reshape to long format and regress each diversity metric on the human-influence index, to visualise the human-influence gradient (panel A) ---
leaf.alpha.long <- leaf.alpha.df %>%
  dplyr::select(ugs, ugs_type, human_index,
                plant_rich, MPD, MNTD,
                human_trait_fdis, ecological_trait_fdis, human_trait_fric, ecological_trait_fric,
                mtb_rich) %>%
  pivot_longer(cols = -c(ugs, ugs_type, human_index), names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("plant_rich", "MPD", "MNTD", "mtb_rich",
                                                   "human_trait_fric", "human_trait_fdis", "ecological_trait_fric", "ecological_trait_fdis"))) %>%
  dplyr::mutate(metric = recode(metric,
                                plant_rich = "Plant richness",
                                mtb_rich = "Feature richness",
                                human_trait_fric = "Human func. richness",
                                human_trait_fdis = "Human func. dispersion",
                                ecological_trait_fric = "Ecological func. richness",
                                ecological_trait_fdis = "Ecological func. dispersion"))

# fit a model per metric: negative-binomial GLM for the (count) plant and feature richness metrics, linear model for the rest
leaf.alpha.reg <- leaf.alpha.long %>%
  group_by(metric) %>%
  dplyr::summarise(model = list(if (unique(metric) %in% c("Plant richness", "Feature richness")) {
    MASS::glm.nb(value ~ std(human_index), data = cur_data())
    } else { 
      lm(value ~ std(human_index), data = cur_data())
      }), 
    .groups = "drop") %>%
  dplyr::mutate(beta = map_dbl(model, ~coef(.x)["std(human_index)"]),
                p.value = map_dbl(model, ~broom::tidy(.x) %>% filter(term == "std(human_index)") %>% pull(p.value)),
                r2 = map2_dbl(model, metric, ~{
                  if (.y %in% c("Plant richness", "Feature richness")) {deviance_explained(.x)} else {summary(.x)$r.squared}
                })) %>%
  dplyr::select(-model) %>%
  dplyr::mutate(signif = p.value < 0.05) %>%
  dplyr::mutate(sig = sig_label(!signif),
                sig_star = case_when(p.value < 0.001 ~ "< 0.001", p.value < 0.01 ~ "< 0.01",
                                     p.value < 0.05 ~ "< 0.05", TRUE ~ "> 0.05"),
                sig_label = ifelse(signif == TRUE, "Significant", "Non-significant"))

leaf.alpha.long <- leaf.alpha.long %>% left_join(leaf.alpha.reg, by = "metric")

# --- 3h. Piecewise SEM: human influence -> trait diversity -> phylogenetic diversity -> chemical richness (panels B and C) ---
leaf.sem <- piecewiseSEM::psem(lm(human_trait_fric ~ human_index, data = na.omit(leaf.alpha.df)),
                               lm(ecological_trait_fric ~ human_index, data = na.omit(leaf.alpha.df)),
                               lm(MPD ~ human_trait_fric + ecological_trait_fric + human_index, data = na.omit(leaf.alpha.df)),
                               lm(MNTD ~ human_trait_fric + ecological_trait_fric, data = na.omit(leaf.alpha.df)),
                               MASS::glm.nb(mtb_rich ~ MPD + MNTD + human_trait_fric + ecological_trait_fric + human_index, data = na.omit(leaf.alpha.df)),
                               human_trait_fric %~~% ecological_trait_fric,
                               MPD %~~% MNTD)
summary(leaf.sem)

leaf.sem.df <- summary(leaf.sem, .progressBar = FALSE)$coefficients
names(leaf.sem.df)[names(leaf.sem.df) == ""] <- "sig"

leaf.sem.df <- leaf.sem.df %>%
  dplyr::mutate(style = case_when(P.Value < 0.05 ~ "solid", P.Value < 0.1 ~ "dashed", TRUE ~ "solid"),
                color = case_when(P.Value < 0.1 ~ "black", TRUE ~ "grey"),
                Std.Estimate = case_when(color == "grey" ~ NA_real_, TRUE ~ round(Std.Estimate, 3)),
                color = case_when(color == "grey" ~ "grey", Std.Estimate < 0 ~ "#d7301f", TRUE ~ "black"),
                color_bg = case_when(color == "#d7301f" ~ "#fddbc7", color == "black" ~ "grey80", TRUE ~ NA_character_)) %>%
  dplyr::mutate(Predictor = recode(Predictor,
                                   human_index = "Human\ninfluence\nindex",
                                   human_trait_fric = "Human\nfunc. richness",
                                   ecological_trait_fric = "Ecological\nfunc. richness"),
                Response = recode(Response,
                                  human_trait_fric = "Human\nfunc. richness",
                                  ecological_trait_fric = "Ecological\nfunc. richness",
                                  mtb_rich = "Leaf chemical\nfeature richness"))

leaf.dag <- igraph::graph_from_adjacency_matrix(getDAG(leaf.sem))
leaf.edge.list <- igraph::as_edgelist(leaf.dag)

leaf.edge.list[leaf.edge.list == "human_index"] <- "Human\ninfluence\nindex"
leaf.edge.list[leaf.edge.list == "human_trait_fric"] <- "Human\nfunc. richness"
leaf.edge.list[leaf.edge.list == "ecological_trait_fric"] <- "Ecological\nfunc. richness"
leaf.edge.list[leaf.edge.list == "mtb_rich"] <- "Leaf chemical\nfeature richness"

leaf.sem.df <- leaf.sem.df[match(paste(leaf.edge.list[,1], leaf.edge.list[,2]), paste(leaf.sem.df$Predictor, leaf.sem.df$Response)), ]

# bootstrap the model to get direct/indirect/total standardized effects at several confidence levels, used below to choose significance stars
leaf.sem.boot <- bootEff(leaf.sem, R = 1000, parallel = "multicore", type = "nonparametric")

leaf.sem.eff90 <- semEff(leaf.sem.boot, ci.type = "norm", ci.conf = 0.90)$Effects$Table
leaf.sem.eff95 <- semEff(leaf.sem.boot, ci.type = "norm", ci.conf = 0.95)$Effects$Table
leaf.sem.eff99 <- semEff(leaf.sem.boot, ci.type = "norm", ci.conf = 0.99)$Effects$Table
leaf.sem.eff999 <- semEff(leaf.sem.boot, ci.type = "norm", ci.conf = 0.999)$Effects$Table

leaf.sem.eff95$lower90 <- leaf.sem.eff90$lower_ci
leaf.sem.eff95$upper90 <- leaf.sem.eff90$upper_ci
leaf.sem.eff95$lower99 <- leaf.sem.eff99$lower_ci
leaf.sem.eff95$upper99 <- leaf.sem.eff99$upper_ci
leaf.sem.eff95$lower999 <- leaf.sem.eff999$lower_ci
leaf.sem.eff95$upper999 <- leaf.sem.eff999$upper_ci

# node positions for the path diagram (a manually laid-out DAG)
leaf.nodes <- tibble(name = c("Human\ninfluence\nindex", "Human\nfunc. richness", "Ecological\nfunc. richness", "MPD", "MNTD", "Leaf chemical\nfeature richness"),
                     x = c(0, 1, 1, 2.5, 2.5, 3.5), 
                     y = c(1, 2, 0, 0, 2, 1))

leaf.edges <- cbind(data.frame(from = leaf.edge.list[,1], to = leaf.edge.list[,2]), leaf.sem.df) %>%
  left_join(leaf.nodes, by = c("from" = "name")) %>%
  left_join(leaf.nodes, by = c("to" = "name")) %>%
  dplyr::rename(x = x.x, y = y.x, xend = x.y, yend = y.y) %>%
  mutate(label_x = x + 0.45 * (xend - x), label_y = y + 0.45 * (yend - y),
         linewidth = ifelse(is.na(Std.Estimate), 0.2, abs(Std.Estimate) * 2))

# semi-axes of the elliptical node shapes drawn on the diagram
leaf.ellipse.a <- 0.42
leaf.ellipse.b <- 0.25

# manually-specified dotted connector lines for the two correlated-error terms (human_trait_fric <-> ecological_trait_fric, MPD <-> MNTD), drawn between the fixed node coordinates above
leaf.corr_edges <- tribble(~x, ~y, ~xend, ~yend, 1, 2, 1, 0, 2.5, 2, 2.5, 0)

# shorten each arrow at both ends so it starts/ends at the node's ellipse boundary rather than at its center (r = ellipse radius in the arrow's direction; pad adds a small extra gap)
leaf.edges <- leaf.edges %>%
  mutate(dx = xend - x, dy = yend - y, theta = atan2(dy, dx), L = sqrt(dx^2 + dy^2),
         r = 1 / sqrt((cos(theta)^2 / leaf.ellipse.a^2) + (sin(theta)^2 / leaf.ellipse.b^2)),
         pad = r + 0.2,
         x2    = x    + pad * dx / L, y2    = y    + pad * dy / L,
         xend2 = xend - pad * dx / L, yend2 = yend - pad * dy / L)

leaf.sem.plot <- ggplot() +
  geom_curve(data = leaf.corr_edges, aes(x, y, xend = xend, yend = yend),
             curvature = 0, linetype = "dotted", linewidth = 0.6, colour = "grey40", arrow = NULL) +
  geom_segment(data = leaf.edges, aes(x = x2, y = y2, xend = xend2, yend = yend2, color = color, linetype = style, linewidth = linewidth),
               arrow = arrow(length = unit(2.5, "mm")), lineend = "round") + 
  geom_label(data = filter(leaf.edges, !is.na(Std.Estimate)),
             aes(label_x, label_y, label = Std.Estimate, fill = color_bg), colour = "black", label.size = 0, size = 3.5) +
  ggforce::geom_ellipse(data = leaf.nodes, aes(x0 = x, y0 = y, a = leaf.ellipse.a, b = leaf.ellipse.b, angle = 0), fill = "white", colour = NA) +
  geom_text(data = leaf.nodes, aes(x, y, label = name), lineheight = 1, size = 5) +
  scale_colour_identity() + scale_fill_identity() + scale_linetype_identity() + scale_linewidth_identity() +
  coord_cartesian(clip = "off") +
  theme_void() + 
  theme(plot.margin = margin(0, 0, 0, 0))

leaf.sem.eff <- leaf.sem.eff95 %>%
  filter(response == "mtb.rich") %>%
  dplyr::mutate(predictor = factor(predictor, levels = c("MNTD", "MPD", "ecological.trait.fric", "human.trait.fric", "human.index"))) %>%
  dplyr::mutate(predictor = recode(predictor, human.index = "Human\ninfluence index",
                                   human.trait.fric = "Human\nfunc. richness",
                                   ecological.trait.fric = "Ecological\nfunc. richness")) %>%
  mutate(stars = case_when(lower999 > 0 | upper999 < 0 ~ "***", lower99 > 0 | upper99 < 0 ~ "**",
                           lower_ci > 0 | upper_ci < 0 ~ "*", lower90 > 0 | upper90 < 0 ~ "\u00b7", TRUE ~ ""),
         signif_cat = case_when(stars != "" ~ "significant", TRUE ~ "ns")) %>%
  mutate(ci_width = upper_ci - lower_ci, star_x = ifelse(effect > 0, upper_ci + 0.06, lower_ci - 0.06)) 


################################################################################
# ---- 4. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 4a. Restrict species-level data to the species profiled in floral VOC sampling ---
flower.species.traits.raw <- species.traits.raw %>% filter(species %in% flower.meta$Plants)
flower.species.mat <- species.mat %>% select(ugs, flower.meta$Plants)
flower.species.taxo <- species.taxo %>% filter(species %in% flower.meta$Plants)

# --- 4b. Restrict to sites with chemistry, floristic, and covariate data, and align row order (see the corresponding note in 3b) ---
flower.keep.sites <- Reduce(intersect, list(flower.chem.mat$ugs, flower.species.mat$ugs, covariates.df$ugs))
flower.chem.mat <- flower.chem.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))
flower.species.mat <- flower.species.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))

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
flower.eco_fd <- compute_site_mfd(ecological_traits, flower.species.comm, flower.traits.df, "ecological_trait")

flower.human_fric <- flower.human_fd %>% dplyr::select(ugs, human_trait_fric)
flower.eco_fric <- flower.eco_fd %>% dplyr::select(ugs, ecological_trait_fric)
flower.human_fdis <- flower.human_fd %>% dplyr::select(ugs, human_trait_fdis)
flower.eco_fdis <- flower.eco_fd %>% dplyr::select(ugs, ecological_trait_fdis)

# --- 4f. Site-level alpha-diversity table for every flower site ---
flower.alpha.df <- tibble(ugs = flower.chem.mat$ugs,
                          voc_rich = vegan::specnumber(flower.chem.mat %>% select(-ugs)),
                          plant_rich = vegan::specnumber(flower.species.comm)) %>%
  left_join(covariates.df, by = "ugs") %>%
  left_join(flower.human_fdis, by = "ugs") %>%
  left_join(flower.eco_fdis, by = "ugs") %>%
  left_join(flower.human_fric, by = "ugs") %>%
  left_join(flower.eco_fric, by = "ugs") %>%
  left_join(flower.phylo.alpha.df, by = "ugs")

# --- 4g. Reshape to long format and regress each metric on human influence ---
flower.alpha.long <- flower.alpha.df %>%
  dplyr::select(ugs, ugs_type, human_index,
                plant_rich, MPD, MNTD,
                human_trait_fdis, ecological_trait_fdis, human_trait_fric, ecological_trait_fric,
                voc_rich) %>%
  pivot_longer(cols = -c(ugs, ugs_type, human_index), names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels = c("plant_rich", "MPD", "MNTD", "voc_rich",
                                                   "human_trait_fric", "human_trait_fdis", "ecological_trait_fric", "ecological_trait_fdis"))) %>%
  dplyr::mutate(metric = recode(metric,
                                plant_rich = "Plant richness",
                                voc_rich = "Feature richness",
                                human_trait_fric = "Human func. richness",
                                human_trait_fdis = "Human func. dispersion",
                                ecological_trait_fric = "Ecological func. richness",
                                ecological_trait_fdis = "Ecological func. dispersion"))

flower.alpha.reg <- flower.alpha.long %>%
  group_by(metric) %>%
  dplyr::summarise(model = list(if (unique(metric) %in% c("Plant richness", "Feature richness")) {
    MASS::glm.nb(value ~ std(human_index), data = cur_data())
    } else { 
      lm(value ~ std(human_index), data = cur_data())
      }), 
    .groups = "drop") %>%
  dplyr::mutate(beta = map_dbl(model, ~coef(.x)["std(human_index)"]),
                p.value = map_dbl(model, ~broom::tidy(.x) %>% filter(term == "std(human_index)") %>% pull(p.value)),
                r2 = map2_dbl(model, metric, ~{
                  if (.y %in% c("Plant richness", "Feature richness")) {deviance_explained(.x)} else {summary(.x)$r.squared}
                })) %>%
  dplyr::select(-model) %>%
  dplyr::mutate(signif = p.value < 0.05) %>%
  dplyr::mutate(sig = sig_label(!signif),
                sig_star = case_when(p.value < 0.001 ~ "< 0.001", p.value < 0.01 ~ "< 0.01",
                                     p.value < 0.05 ~ "< 0.05", TRUE ~ "> 0.05"),
                sig_label = ifelse(signif == TRUE, "Significant", "Non-significant"))

flower.alpha.long <- flower.alpha.long %>% left_join(flower.alpha.reg, by = "metric")

# --- 4h. Piecewise SEM for flower ---
flower.sem <- piecewiseSEM::psem(lm(human_trait_fric ~ human_index, data = na.omit(flower.alpha.df)),
                                 lm(ecological_trait_fric ~ human_index, data = na.omit(flower.alpha.df)),
                                 lm(MPD ~ human_trait_fric + ecological_trait_fric + human_index, data = na.omit(flower.alpha.df)),
                                 lm(MNTD ~ human_trait_fric + ecological_trait_fric, data = na.omit(flower.alpha.df)),
                                 MASS::glm.nb(voc_rich ~ MPD + MNTD + human_trait_fric + ecological_trait_fric + human_index, data = na.omit(flower.alpha.df)),
                                 human_trait_fric %~~% ecological_trait_fric,
                                 MPD %~~% MNTD)
summary(flower.sem)

flower.sem.df <- summary(flower.sem, .progressBar = FALSE)$coefficients
names(flower.sem.df)[names(flower.sem.df) == ""] <- "sig"

flower.sem.df <- flower.sem.df %>%
  dplyr::mutate(style = case_when(P.Value < 0.05 ~ "solid", P.Value < 0.1 ~ "dashed", TRUE ~ "solid"),
                color = case_when(P.Value < 0.1 ~ "black", TRUE ~ "grey"),
                Std.Estimate = case_when(color == "grey" ~ NA_real_, TRUE ~ round(Std.Estimate, 3)),
                color = case_when(color == "grey" ~ "grey", Std.Estimate < 0 ~ "#d7301f", TRUE ~ "black"),
                color_bg = case_when(color == "#d7301f" ~ "#fddbc7", color == "black" ~ "grey80", TRUE ~ NA_character_)) %>%
  dplyr::mutate(Predictor = recode(Predictor,
                                   human_index = "Human\ninfluence\nindex",
                                   human_trait_fric = "Human\nfunc. richness",
                                   ecological_trait_fric = "Ecological\nfunc. richness"),
                Response = recode(Response,
                                  human_trait_fric = "Human\nfunc. richness",
                                  ecological_trait_fric = "Ecological\nfunc. richness",
                                  voc_rich = "Flower chemical\nfeature richness"))

flower.dag <- igraph::graph_from_adjacency_matrix(getDAG(flower.sem))
flower.edge.list <- igraph::as_edgelist(flower.dag)

flower.edge.list[flower.edge.list == "human_index"] <- "Human\ninfluence\nindex"
flower.edge.list[flower.edge.list == "human_trait_fric"] <- "Human\nfunc. richness"
flower.edge.list[flower.edge.list == "ecological_trait_fric"] <- "Ecological\nfunc. richness"
flower.edge.list[flower.edge.list == "voc_rich"] <- "Flower chemical\nfeature richness"

flower.sem.df <- flower.sem.df[match(paste(flower.edge.list[,1], flower.edge.list[,2]), paste(flower.sem.df$Predictor, flower.sem.df$Response)), ]

flower.sem.boot <- bootEff(flower.sem, R = 1000, parallel = "multicore", type = "nonparametric")

flower.sem.eff90 <- semEff(flower.sem.boot, ci.type = "norm", ci.conf = 0.90)$Effects$Table
flower.sem.eff95 <- semEff(flower.sem.boot, ci.type = "norm", ci.conf = 0.95)$Effects$Table
flower.sem.eff99 <- semEff(flower.sem.boot, ci.type = "norm", ci.conf = 0.99)$Effects$Table
flower.sem.eff999 <- semEff(flower.sem.boot, ci.type = "norm", ci.conf = 0.999)$Effects$Table

flower.sem.eff95$lower90 <- flower.sem.eff90$lower_ci
flower.sem.eff95$upper90 <- flower.sem.eff90$upper_ci
flower.sem.eff95$lower99 <- flower.sem.eff99$lower_ci
flower.sem.eff95$upper99 <- flower.sem.eff99$upper_ci
flower.sem.eff95$lower999 <- flower.sem.eff999$lower_ci
flower.sem.eff95$upper999 <- flower.sem.eff999$upper_ci

flower.nodes <- tibble(name = c("Human\ninfluence\nindex", "Human\nfunc. richness", "Ecological\nfunc. richness", "MPD", "MNTD", "Flower chemical\nfeature richness"),
                       x = c(0, 1, 1, 2.5, 2.5, 3.5), 
                       y = c(1, 2, 0, 0, 2, 1))

flower.edges <- cbind(data.frame(from = flower.edge.list[,1], to = flower.edge.list[,2]), flower.sem.df) %>%
  left_join(flower.nodes, by = c("from" = "name")) %>%
  left_join(flower.nodes, by = c("to" = "name")) %>%
  dplyr::rename(x = x.x, y = y.x, xend = x.y, yend = y.y) %>%
  mutate(label_x = x + 0.45 * (xend - x), label_y = y + 0.45 * (yend - y),
         linewidth = ifelse(is.na(Std.Estimate), 0.2, abs(Std.Estimate) * 2))

flower.ellipse.a <- 0.42
flower.ellipse.b <- 0.25

flower.corr_edges <- tribble(~x, ~y, ~xend, ~yend, 1, 2, 1, 0, 2.5, 2, 2.5, 0)

flower.edges <- flower.edges %>%
  mutate(dx = xend - x, dy = yend - y, theta = atan2(dy, dx), L = sqrt(dx^2 + dy^2),
         r = 1 / sqrt((cos(theta)^2 / leaf.ellipse.a^2) + (sin(theta)^2 / flower.ellipse.b^2)),
         pad = r + 0.2,
         x2    = x    + pad * dx / L, y2    = y    + pad * dy / L,
         xend2 = xend - pad * dx / L, yend2 = yend - pad * dy / L)

flower.sem.plot <- ggplot() +
  geom_curve(data = flower.corr_edges, aes(x, y, xend = xend, yend = yend),
             curvature = 0, linetype = "dotted", linewidth = 0.6, colour = "grey40", arrow = NULL) +
  geom_segment(data = flower.edges, aes(x = x2, y = y2, xend = xend2, yend = yend2, color = color, linetype = style, linewidth = linewidth),
               arrow = arrow(length = unit(2.5, "mm")), lineend = "round") + 
  geom_label(data = filter(flower.edges, !is.na(Std.Estimate)),
             aes(label_x, label_y, label = Std.Estimate, fill = color_bg), colour = "black", label.size = 0, size = 3.5) +
  ggforce::geom_ellipse(data = flower.nodes, aes(x0 = x, y0 = y, a = flower.ellipse.a, b = flower.ellipse.b, angle = 0), fill = "white", colour = NA) +
  geom_text(data = flower.nodes, aes(x, y, label = name), lineheight = 1, size = 5) +
  scale_colour_identity() + scale_fill_identity() + scale_linetype_identity() + scale_linewidth_identity() +
  coord_cartesian(clip = "off") +
  theme_void() + 
  theme(plot.margin = margin(0, 0, 0, 0))

flower.sem.eff <- flower.sem.eff95 %>%
  filter(response == "voc.rich") %>%
  dplyr::mutate(predictor = factor(predictor, levels = c("MNTD", "MPD", "ecological.trait.fric", "human.trait.fric", "human.index"))) %>%
  dplyr::mutate(predictor = recode(predictor, human.index = "Human\ninfluence index",
                                   human.trait.fric = "Human\nfunc. richness",
                                   ecological.trait.fric = "Ecological\nfunc. richness")) %>%
  mutate(stars = case_when(lower999 > 0 | upper999 < 0 ~ "***", lower99 > 0 | upper99 < 0 ~ "**",
                           lower_ci > 0 | upper_ci < 0 ~ "*", lower90 > 0 | upper90 < 0 ~ "\u00b7", TRUE ~ ""),
         signif_cat = case_when(stars != "" ~ "significant", TRUE ~ "ns")) %>%
  mutate(ci_width = upper_ci - lower_ci, star_x = ifelse(effect > 0, upper_ci + 0.06, lower_ci - 0.06)) 


# ---- 5. Combine tissues ----------------------------------------------------

alpha.long <- rbind(data.frame(leaf.alpha.long %>%
                                 dplyr::filter(metric %in% c("MPD", "MNTD", "Feature richness", "Human func. richness", "Ecological func. richness")) %>%
                                 dplyr::group_by(metric) %>%
                                 dplyr::mutate(value_scaled = bruceR::scaler(value)),
                               tissue = "leaf"),
                    data.frame(flower.alpha.long %>%
                                 dplyr::filter(metric %in% c("MPD", "MNTD", "Feature richness", "Human func. richness", "Ecological func. richness")) %>%
                                 dplyr::group_by(metric) %>%
                                 dplyr::mutate(value_scaled = bruceR::scaler(value)),
                               tissue = "flower")) %>%
  mutate(tissue = factor(tissue, levels = c("leaf", "flower")),
         metric = factor(metric, levels = c("Feature richness", "Human func. richness", "Ecological func. richness", "MPD", "MNTD"))) %>%
  mutate(tissue = recode(tissue, "leaf" = "Leaf", "flower" = "Flower"))

alpha.reg <- rbind(data.frame(leaf.alpha.reg   %>% dplyr::filter(metric %in% c("MPD", "MNTD", "Feature richness", "Human func. richness", "Ecological func. richness")), tissue = "Flower"),
                   data.frame(flower.alpha.reg %>% dplyr::filter(metric %in% c("MPD", "MNTD", "Feature richness", "Human func. richness", "Ecological func. richness")), tissue = "Leaf"))

sem.eff <- rbind(data.frame(leaf.sem.eff, tissue = "leaf"), 
                 data.frame(flower.sem.eff, tissue = "flower")) %>%
  mutate(tissue = factor(tissue, levels = c("leaf", "flower")))

# ---- 6. Build the figure's three panels ------------------------------------

# Panel A: chemical feature richness vs. human-influence index, one facet per tissue
richness.plot <- ggplot(alpha.long %>% dplyr::filter(metric == "Feature richness"), aes(x = human_index, y = value, shape = tissue)) +
  geom_point(aes(color = ugs_type), size = 2, alpha = 0.35) +
  scale_color_manual(values = site.colors, name = NULL) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 1, color = "black") +
  geom_text(data = alpha.reg %>% dplyr::filter(metric == "Feature richness"),
            aes(x = Inf, y = -Inf, label = paste0(sprintf("Deviance explained = %.2f", r2), "; p < 0.0001 ")),
            inherit.aes = FALSE, hjust = 1, vjust = -1, size = 4.5, show.legend = FALSE) +
  facet_wrap(~ factor(tissue, levels = c("Leaf", "Flower")), scales = "free_y", ncol = 2, strip.position = "top") +
  labs(x = "Human-influence index",
       y = "Chemical\nfeature richness",
       alpha = NULL,
       linetype = NULL,
       shape = NULL) +
  theme_pubr(base_size = 15) +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        legend.position = "top",
        legend.box = "vertical",
        legend.spacing = unit(0, "pt")) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         fill = guide_legend(nrow = 1, byrow = TRUE)) + 
  guides(shape = "none")

# Panel B: leaf and flower SEM path diagrams, stacked
diagram.plot <- ggarrange(leaf.sem.plot, flower.sem.plot, 
                          ncol = 1, 
                          font.label = list(size = 17))

# Panel C: standardized direct/indirect/total effects on chemical richness, one row per predictor, split by tissue
predictor_levels <- levels(factor(sem.eff$predictor)) # row order
tissue_levels <- c("leaf", "flower") # keep in sync with site.colors' tissue order
n_tissue <- length(tissue_levels)

group_width <- 0.8
bar_height  <- 0.3

sem.eff <- sem.eff %>%
  mutate(predictor = factor(predictor, levels = predictor_levels),
         tissue = factor(tissue, levels = tissue_levels),
         y_pos = as.numeric(predictor) + group_width * ((as.numeric(tissue) - 0.5) / n_tissue - 0.5))

effects.plot <- ggplot(sem.eff, aes(x = effect, shape = tissue)) +
  ggpattern::geom_col_pattern(data = sem.eff %>% filter(effect_type %in% c("direct", "indirect")),
                              aes(x = effect, y = y_pos, pattern = effect_type, group = interaction(tissue, effect_type)),
                              orientation = "y",            
                              position = "stack", 
                              color = "black", 
                              fill = "grey85",
                              alpha = 0.2, width = bar_height,
                              pattern_size = 0.35, pattern_angle = 45,
                              pattern_density = 0.01, pattern_spacing = 0.05,
                              pattern_key_scale_factor = .1, 
                              show.legend = c("pattern" = TRUE, "pattern_colour" = FALSE, "color" = FALSE, "fill" = FALSE)) +
  geom_point(data = sem.eff %>% filter(effect_type == "total"),
             aes(x = effect, y = y_pos, alpha = signif_cat), size = 4) +
  geom_errorbar(data = sem.eff %>% filter(effect_type == "total"),
                aes(xmin = lower_ci, xmax = upper_ci, y = y_pos, alpha = signif_cat),
                width = 0, linewidth = 1.2) +
  geom_text(data = sem.eff %>% filter(effect_type == "total"),
            aes(x = star_x, y = y_pos, label = stars),
            vjust = 0.75, size = 6, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40") +
  scale_y_continuous(breaks = seq_along(predictor_levels),
                     labels = predictor_levels,
                     expand = expansion(add = 0.1)) +
  ggpattern::scale_pattern_manual(values = c(direct = "none", indirect = "stripe")) +
  scale_alpha_manual(values = c(significant = 1, marginal = 1, ns = 0.3), guide = "none") +
  labs(x = "Standardized coefficients\nand 95% CI", y = "",
       alpha = NULL, 
       shape = "Feature type:", 
       color = "Tissue",
       fill = "Tissue",
       pattern = "Effect type:") +
  theme_pubr(base_size = 15) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 14),
        plot.margin = margin(0, 0, 0, 10),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt"),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom") +
  guides(pattern = guide_legend(override.aes = list(fill = "white", 
                                                    color = "black", 
                                                    pattern_fill = "black",
                                                    pattern_color = "black", 
                                                    pattern_density = 0.01, 
                                                    pattern_spacing = 0.08)), 
         color = "none", fill = "none")


# ---- 7. Assemble and save the figure ---------------------------------------

ggarrange(richness.plot, 
          ggarrange(diagram.plot, 
                    effects.plot + theme(legend.position = "bottom"), 
                    labels = c("B", "C"), 
                    font.label = list(size = 17),
                    widths = c(1, 0.9),
                    common.legend = TRUE), 
          labels = c("A", "", "", ""), 
          font.label = list(size = 17), 
          heights = c(0.6, 1),
          nrow = 2) + 
  theme(panel.background = element_rect(fill = "white", color = "white"))

ggsave("figures/fig5.png", height = 9, width = 9)
ggsave("figures/fig5.pdf", height = 9, width = 9)
