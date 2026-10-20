library(dplyr)
library(tibble)
library(V.PhyloMaker2)
library(vegan)
library(ape)
library(picante)

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

# ---- 2. Trait definitions (shared) ------------------------------------

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

# Ordinal traits are stored in the raw data as ordered text categories; these maps convert each category to its rank (1 = lowest) so the trait can be used in the Gower dissimilarity below as ordered/numeric rather than nominal
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


################################################################################
# ---- 3. LEAF analysis -------------------------------------------------------
################################################################################

# --- 3a. Restrict species-level data to the species profiled in leaf metabolomics ---
leaf.species.traits.raw <- species.traits.raw %>% filter(species %in% leaf.meta$Plants)
leaf.species.mat        <- species.mat %>% select(ugs, leaf.meta$Plants)
leaf.species.taxo       <- species.taxo %>% filter(species %in% leaf.meta$Plants)

# --- 3b. Restrict to sites with chemistry, floristic, and covariate data, so the resulting distance matrices only cover sites usable in Figures 2/3 ---
leaf.keep.sites  <- Reduce(intersect, list(leaf.chem.mat$ugs, leaf.species.mat$ugs, covariates.df$ugs))
leaf.species.mat <- leaf.species.mat %>% filter(ugs %in% leaf.keep.sites) %>% arrange(match(ugs, leaf.keep.sites))

# --- 3c. Clean the species trait table: convert ordinal traits to ranks, and collapse a few sparse/rare categories into a broader neighbouring category ---
leaf.traits.df <- leaf.species.traits.raw %>%
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  column_to_rownames("species")

for (nm in intersect(names(ordinal_maps), colnames(leaf.traits.df))) {
  leaf.traits.df[[nm]] <- ordinal_maps[[nm]][as.character(leaf.traits.df[[nm]])]
  leaf.traits.df[[nm]] <- as.ordered(leaf.traits.df[[nm]])
}

leaf.traits.df <- leaf.traits.df %>%
  dplyr::mutate(Flower.color = case_when(Flower.color %in% ("orange") ~ "red",
                                         Flower.color %in% c("pink/purple/violet") ~ "blue",
                                         Flower.color %in% c("brown/beige/gray") ~ "white",
                                         TRUE ~ Flower.color),
                Life.span = case_when(Life.span %in% c("annual", "biennial") ~ 1, TRUE ~ 0),
                Flower.bilateral = ifelse(Flower.bilateral == "bilateral", 1, 0),
                Leaf.greenness = ifelse(Leaf.greenness == "green", 1, 0),
                Leaf.elliptic = ifelse(Leaf.elliptic == "elliptic", 1, 0),
                Leaf.compoundness = ifelse(Leaf.compoundness == "compound", 1, 0),
                Leaf.deciduous = ifelse(Leaf.deciduous == "deciduous", 1, 0),
                Fall.Conspicuous = ifelse(Fall.Conspicuous == "Fall Conspicuous", 1, 0)) %>%
  dplyr::mutate(Flower.color = as.factor(Flower.color))

# --- 3d. Build the site x species community matrix, restricted to species present in both the floristic survey and the cleaned trait table ---
leaf.candidate_species <- intersect(leaf.species.taxo$species[leaf.species.taxo$class == "Magnoliopsida"],
                                    intersect(colnames(leaf.species.mat), rownames(leaf.traits.df)))
leaf.traits.df <- leaf.traits.df[leaf.candidate_species, , drop = FALSE]

leaf.species.comm <- leaf.species.mat %>%
  dplyr::select(all_of(leaf.candidate_species)) %>%
  as.matrix()
rownames(leaf.species.comm) <- leaf.species.mat$ugs

# --- 3e. Build a phylogeny for the candidate species (V.PhyloMaker2) and trim the community matrix / trait table down to species actually placed on the tree ---
leaf.tree.scenarios <- V.PhyloMaker2::phylo.maker(sp.list = leaf.species.taxo %>% dplyr::select(species, genus, family),
                                                  scenarios = c("S1", "S2", "S3"))

# scenario 3: species missing from the mega-tree are bound in as polytomies
# at the genus (or family) level -- the standard, most inclusive choice
leaf.tree <- leaf.tree.scenarios$scenario.3
leaf.tree$tip.label <- gsub("_", " ", leaf.tree$tip.label)
leaf.tree$tip.label <- gsub("  ", " _", leaf.tree$tip.label)

leaf.tree.species <- intersect(colnames(leaf.species.comm), leaf.tree$tip.label)
leaf.tree <- ape::keep.tip(leaf.tree, leaf.tree.species)
leaf.species.comm <- leaf.species.comm[, leaf.tree$tip.label, drop = FALSE]
leaf.traits.df <- leaf.traits.df[leaf.tree$tip.label, , drop = FALSE]

leaf.coph.mat <- ape::cophenetic.phylo(leaf.tree) # pairwise phylogenetic (branch-length) distance between species
leaf.coph.mat <- leaf.coph.mat[colnames(leaf.species.comm), colnames(leaf.species.comm)]

# --- 3f. Phylogenetic community dissimilarity between sites ---
leaf.phylo.dist.comdist <- picante::comdist(leaf.species.comm, leaf.coph.mat, abundance.weighted = FALSE)   # mean pairwise distance (MPD-based)
leaf.phylo.dist.comdistnt <- picante::comdistnt(leaf.species.comm, leaf.coph.mat, abundance.weighted = FALSE) # mean nearest-taxon distance (MNTD-based)

# --- 3g. Functional (trait-based) community dissimilarity between sites,
# for three trait subsets: all traits, human-associated only, ecological
# only. Gower distance handles the mix of numeric, ordinal, and nominal
# trait types. ---
leaf.gower.all <- cluster::daisy(leaf.traits.df %>% dplyr::select(any_of(c(human_traits, eco_traits))), metric = "gower")
leaf.gower.human <- cluster::daisy(leaf.traits.df %>% dplyr::select(any_of(human_traits)), metric = "gower")
leaf.gower.eco <- cluster::daisy(leaf.traits.df %>% dplyr::select(any_of(eco_traits)), metric = "gower")

leaf.all.dist.comdist <- picante::comdist(leaf.species.comm, as.matrix(leaf.gower.all),   abundance.weighted = FALSE)
leaf.human.dist.comdist <- picante::comdist(leaf.species.comm, as.matrix(leaf.gower.human), abundance.weighted = FALSE)
leaf.eco.dist.comdist <- picante::comdist(leaf.species.comm, as.matrix(leaf.gower.eco),   abundance.weighted = FALSE)

leaf.all.dist.comdistnt <- picante::comdistnt(leaf.species.comm, as.matrix(leaf.gower.all),   abundance.weighted = FALSE)
leaf.human.dist.comdistnt <- picante::comdistnt(leaf.species.comm, as.matrix(leaf.gower.human), abundance.weighted = FALSE)
leaf.eco.dist.comdistnt <- picante::comdistnt(leaf.species.comm, as.matrix(leaf.gower.eco),   abundance.weighted = FALSE)

################################################################################
# ---- 4. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 4a. Restrict species-level data to the species profiled in floral VOC sampling ---
flower.species.traits.raw <- species.traits.raw %>% filter(species %in% flower.meta$Plants)
flower.species.mat <- species.mat %>% select(ugs, flower.meta$Plants)
flower.species.taxo <- species.taxo %>% filter(species %in% flower.meta$Plants)

# --- 4b. Restrict to sites with chemistry, floristic, and covariate data ---
flower.keep.sites  <- Reduce(intersect, list(flower.chem.mat$ugs, flower.species.mat$ugs, covariates.df$ugs))
flower.species.mat <- flower.species.mat %>% filter(ugs %in% flower.keep.sites) %>% arrange(match(ugs, flower.keep.sites))

# --- 4c. Clean the species trait table ---
flower.traits.df <- flower.species.traits.raw %>%
  dplyr::mutate(across(where(is.character), as.factor)) %>%
  column_to_rownames("species")

for (nm in intersect(names(ordinal_maps), colnames(flower.traits.df))) {
  flower.traits.df[[nm]] <- ordinal_maps[[nm]][as.character(flower.traits.df[[nm]])]
  flower.traits.df[[nm]] <- as.ordered(flower.traits.df[[nm]])
}

flower.traits.df <- flower.traits.df %>%
  dplyr::mutate(Flower.color = case_when(Flower.color %in% ("orange") ~ "red",
                                         Flower.color %in% c("pink/purple/violet") ~ "blue",
                                         Flower.color %in% c("brown/beige/gray") ~ "white",
                                         TRUE ~ Flower.color),
                Life.span = case_when(Life.span %in% c("annual", "biennial") ~ 1, TRUE ~ 0),
                Flower.bilateral = ifelse(Flower.bilateral == "bilateral", 1, 0),
                Leaf.greenness = ifelse(Leaf.greenness == "green", 1, 0),
                Leaf.elliptic = ifelse(Leaf.elliptic == "elliptic", 1, 0),
                Leaf.compoundness = ifelse(Leaf.compoundness == "compound", 1, 0),
                Leaf.deciduous = ifelse(Leaf.deciduous == "deciduous", 1, 0),
                Fall.Conspicuous = ifelse(Fall.Conspicuous == "Fall Conspicuous", 1, 0)) %>%
  dplyr::mutate(Flower.color = as.factor(Flower.color))

# --- 4d. Build the site x species community matrix ---
flower.candidate_species <- intersect(flower.species.taxo$species[flower.species.taxo$class == "Magnoliopsida"],
                                      intersect(colnames(flower.species.mat), rownames(flower.traits.df)))
flower.traits.df <- flower.traits.df[flower.candidate_species, , drop = FALSE]

flower.species.comm <- flower.species.mat %>%
  dplyr::select(all_of(flower.candidate_species)) %>%
  as.matrix()
rownames(flower.species.comm) <- flower.species.mat$ugs

# --- 4e. Build a phylogeny and trim community/trait data to the tree's tips ---
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

# --- 4f. Phylogenetic community dissimilarity between sites ---
flower.phylo.dist.comdist <- picante::comdist(flower.species.comm, flower.coph.mat, abundance.weighted = FALSE)
flower.phylo.dist.comdistnt <- picante::comdistnt(flower.species.comm, flower.coph.mat, abundance.weighted = FALSE)

# --- 4g. Functional (trait-based) community dissimilarity between sites ---
flower.gower.all <- cluster::daisy(flower.traits.df %>% dplyr::select(any_of(c(human_traits, eco_traits))), metric = "gower")
flower.gower.human <- cluster::daisy(flower.traits.df %>% dplyr::select(any_of(human_traits)), metric = "gower")
flower.gower.eco <- cluster::daisy(flower.traits.df %>% dplyr::select(any_of(eco_traits)), metric = "gower")

flower.all.dist.comdist <- picante::comdist(flower.species.comm, as.matrix(flower.gower.all),   abundance.weighted = FALSE)
flower.human.dist.comdist <- picante::comdist(flower.species.comm, as.matrix(flower.gower.human), abundance.weighted = FALSE)
flower.eco.dist.comdist <- picante::comdist(flower.species.comm, as.matrix(flower.gower.eco),   abundance.weighted = FALSE)

flower.all.dist.comdistnt <- picante::comdistnt(flower.species.comm, as.matrix(flower.gower.all),   abundance.weighted = FALSE)
flower.human.dist.comdistnt <- picante::comdistnt(flower.species.comm, as.matrix(flower.gower.human), abundance.weighted = FALSE)
flower.eco.dist.comdistnt <- picante::comdistnt(flower.species.comm, as.matrix(flower.gower.eco),   abundance.weighted = FALSE)

# ---- 5. Save all distance matrices -----------------------------------------
save(leaf.phylo.dist.comdist, leaf.phylo.dist.comdistnt,
     flower.phylo.dist.comdist, flower.phylo.dist.comdistnt,
     leaf.all.dist.comdist, leaf.all.dist.comdistnt, leaf.human.dist.comdist, leaf.eco.dist.comdist, leaf.human.dist.comdistnt, leaf.eco.dist.comdistnt,
     flower.all.dist.comdist, flower.all.dist.comdistnt, flower.human.dist.comdist, flower.eco.dist.comdist, flower.human.dist.comdistnt, flower.eco.dist.comdistnt,
     file = "02_processed/phyto_funct_dist.RData")
