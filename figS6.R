library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ape)
library(V.PhyloMaker2)
library(ggtext)
library(vegan)

set.seed(123) # for reproducibility of the GAM smooths fitted in the plot below

# Color palette distinguishing the two tissues
tissue.colors <- c("Leaf"   = "#007C91",
                   "Flower" = "#8F2D56")


# ---- 1. Load input data --------------------------------------------------

# species x compound matrices (raw peak intensities, one row per species, with a "species" ID column)
# To recreate these datasets, use the prepare_metabolomics_df.R script. 
leaf.compound.mat   <- read.csv("02_processed/leaf_mtb_mat.csv")
flower.compound.mat <- read.csv("02_processed/flower_mtb_mat.csv")

# binarized (presence/absence) version of the leaf/flower compound matrices -- used below to compute species-pairwise chemical distance via Jaccard
leaf.compound.pa <- leaf.compound.mat %>%
  column_to_rownames("species") %>%
  dplyr::mutate(across(everything(), ~ as.numeric(. > 0))) %>%
  as.matrix()

flower.compound.pa <- flower.compound.mat %>%
  column_to_rownames("species") %>%
  dplyr::mutate(across(everything(), ~ as.numeric(. > 0))) %>%
  as.matrix()

# metadata listing which plant species were actually profiled for each tissue
leaf.meta <- read.csv("01_raw/leaf_metadata.csv")
flower.meta <- read.csv("01_raw/flower_metadata.csv")

species.taxo <- read.csv("01_raw/flora_taxo.csv", check.names = FALSE)  # species-level taxonomy (used to restrict to dicots and to build the phylogeny)
species.traits.raw <- read.csv("01_raw/flora_traits.csv", check.names = FALSE) # raw species-level trait table

leaf.species.taxo <- species.taxo %>% filter(species %in% leaf.meta$Plants)
leaf.species.traits.raw <- species.traits.raw %>% filter(species %in% leaf.meta$Plants)

flower.species.taxo <- species.taxo %>% filter(species %in% flower.meta$Plants)
flower.species.traits.raw <- species.traits.raw %>% filter(species %in% flower.meta$Plants)

# ---- 2. Trait definitions and preprocessing (shared) -----------------------

# Ordinal traits are stored in the raw data as ordered text categories; these maps convert each category to its rank (1 = lowest) so the trait is usable in the Gower dissimilarity below as ordered/numeric rather than nominal
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

# Cleans a raw species x trait table: converts character columns to factors, recodes ordinal traits to ranks, and collapses a few sparse/rare categories into a broader neighbouring category. 
# Unlike Figures 4/5, this uses every available trait column directly (no human/ecological subsetting).
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

# --- 3a. Clean species traits ---
leaf.traits.df <- process_traits(leaf.species.traits.raw)

# --- 3b. Restrict to species present in both the taxonomy table and the cleaned trait table, then 
# filter the compound presence/absence matrix down to those species and drop compounds with no variation among them ---
leaf.candidate_species <- intersect(leaf.species.taxo$species[leaf.species.taxo$class == "Magnoliopsida"],
                                    rownames(leaf.traits.df))
leaf.traits.df <- leaf.traits.df[leaf.candidate_species, , drop = FALSE]

leaf.compound.pa <- leaf.compound.pa[leaf.candidate_species, , drop = FALSE]
leaf.compound.pa <- leaf.compound.pa[, apply(leaf.compound.pa, 2, sd, na.rm = TRUE) > 0]

# --- 3c. Species-pairwise chemical distance (Jaccard, on presence/absence) ---
leaf.compound.dist <- as.matrix(vegdist(leaf.compound.pa, method = "jaccard", binary = TRUE))
leaf.lower.tri <- lower.tri(leaf.compound.dist)

# --- 3d. Build a phylogeny for the candidate species and compute species-pairwise phylogenetic (branch-length) distance ---
leaf.tree.scenarios <- phylo.maker(sp.list = leaf.species.taxo %>% dplyr::select(species, genus, family) %>% filter(species %in% leaf.candidate_species),
                                   output.tree = TRUE,
                                   scenarios = c("S1", "S2", "S3"))

# scenario 3: species missing from the mega-tree are bound in as polytomies at the genus (or family) level -- the standard, most inclusive choice
leaf.tree <- leaf.tree.scenarios$scenario.3
leaf.tree$tip.label <- gsub("_", " ", leaf.tree$tip.label)
leaf.tree$tip.label <- gsub("  ", " _", leaf.tree$tip.label)

leaf.phylo.dist <- as.matrix(cophenetic(leaf.tree))
leaf.phylo.dist <- leaf.phylo.dist[leaf.candidate_species, leaf.candidate_species]

# --- 3e. Species-pairwise functional-trait distance (Gower) ---
leaf.traits.dist <- as.matrix(cluster::daisy(leaf.traits.df, metric = "gower"))
leaf.traits.dist <- leaf.traits.dist[leaf.candidate_species, leaf.candidate_species]

################################################################################
# ---- 4. FLOWER analysis ------------------------------------------------------
################################################################################

# --- 4a. Clean species traits ---
flower.traits.df <- process_traits(flower.species.traits.raw)

# --- 4b. Candidate species: every species present in the trait table---
flower.candidate_species <- intersect(flower.species.taxo$species[flower.species.taxo$class == "Magnoliopsida"],
                                      rownames(flower.traits.df))
flower.traits.df <- flower.traits.df[flower.candidate_species, , drop = FALSE]

flower.compound.pa <- flower.compound.pa[flower.candidate_species, , drop = FALSE]
flower.compound.pa <- flower.compound.pa[, apply(flower.compound.pa, 2, sd, na.rm = TRUE) > 0]

# --- 4c. Species-pairwise chemical distance ---
flower.compound.dist <- as.matrix(vegdist(flower.compound.pa, method = "jaccard", binary = TRUE))

colnames(flower.compound.dist) <- flower.candidate_species
rownames(flower.compound.dist) <- flower.candidate_species
flower.lower.tri <- lower.tri(flower.compound.dist)

# --- 4d. Build a phylogeny and compute species-pairwise phylogenetic distance ---
flower.tree.scenarios <- phylo.maker(sp.list = flower.species.taxo %>% dplyr::select(species, genus, family) %>% filter(species %in% flower.candidate_species),
                                     output.tree = TRUE,
                                     scenarios = c("S1", "S2", "S3"))
flower.tree <- flower.tree.scenarios$scenario.3
flower.tree$tip.label <- gsub("_", " ", flower.tree$tip.label)
flower.tree$tip.label <- gsub("  ", " _", flower.tree$tip.label)

flower.phylo.dist <- as.matrix(cophenetic(flower.tree))
flower.phylo.dist <- flower.phylo.dist[flower.candidate_species, flower.candidate_species]

# --- 4e. Species-pairwise functional-trait distance (Gower) ---
flower.traits.dist <- as.matrix(cluster::daisy(flower.traits.df, metric = "gower"))
flower.traits.dist <- flower.traits.dist[flower.candidate_species, flower.candidate_species]

# ---- 5. Combine into one long-format pairwise-distance table ---------------

pairwise_distances.df <- rbind(data.frame(species_1 = rownames(leaf.compound.dist)[row(leaf.compound.dist)[leaf.lower.tri]],
                                          species_2 = colnames(leaf.compound.dist)[col(leaf.compound.dist)[leaf.lower.tri]],
                                          metabolite_distance = leaf.compound.dist[leaf.lower.tri],
                                          distance = leaf.phylo.dist[leaf.lower.tri],
                                          distance_type = "Phylogenetic distance",
                                          tissue = "Leaf"),
                               data.frame(species_1 = rownames(flower.compound.dist)[row(flower.compound.dist)[flower.lower.tri]],
                                          species_2 = colnames(flower.compound.dist)[col(flower.compound.dist)[flower.lower.tri]],
                                          metabolite_distance = flower.compound.dist[flower.lower.tri],
                                          distance = flower.phylo.dist[flower.lower.tri],
                                          distance_type = "Phylogenetic distance",
                                          tissue = "Flower"),
                               data.frame(species_1 = rownames(leaf.compound.dist)[row(leaf.compound.dist)[leaf.lower.tri]],
                                          species_2 = colnames(leaf.compound.dist)[col(leaf.compound.dist)[leaf.lower.tri]],
                                          metabolite_distance = leaf.compound.dist[leaf.lower.tri],
                                          distance = leaf.traits.dist[leaf.lower.tri],
                                          distance_type = "Functional distance",
                                          tissue = "Leaf"),
                               data.frame(species_1 = rownames(flower.compound.dist)[row(flower.compound.dist)[flower.lower.tri]],
                                          species_2 = colnames(flower.compound.dist)[col(flower.compound.dist)[flower.lower.tri]],
                                          metabolite_distance = flower.compound.dist[flower.lower.tri],
                                          distance = flower.traits.dist[flower.lower.tri],
                                          distance_type = "Functional distance",
                                          tissue = "Flower")) %>%  
  mutate(tissue = factor(tissue, levels = c("Leaf", "Flower")))

# ---- 6. Mantel tests and result labels --------------------------------------

leaf.phylo.mantel <- mantel(leaf.compound.dist,   as.dist(leaf.phylo.dist),   permutations = 9999)
flower.phylo.mantel <- mantel(flower.compound.dist, as.dist(flower.phylo.dist), permutations = 9999)
leaf.traits.mantel <- mantel(leaf.compound.dist,   as.dist(leaf.traits.dist),  permutations = 9999)
flower.traits.mantel <- mantel(flower.compound.dist, as.dist(flower.traits.dist), permutations = 9999)

# rich-text labels (leaf's Mantel r in leaf's colour, flower's in flower's color), one per facet (distance_type), positioned in the plot's corner
mantel_labels.df <- data.frame(distance = Inf,
                               metabolite_distance = -Inf,
                               distance_type = c("Phylogenetic distance", "Functional distance"),
                               label = c(sprintf("Mantel r = <span style='color:%s;'>%.3f</span> / <span style='color:%s;'>%.3f</span><br> p = %.4f",
                                                 tissue.colors["Leaf"], leaf.phylo.mantel$statistic,
                                                 tissue.colors["Flower"], flower.phylo.mantel$statistic,
                                                 leaf.phylo.mantel$signif),
                                         sprintf("Mantel r = <span style='color:%s;'>%.3f</span> / <span style='color:%s;'>%.3f</span><br> p = %.4f",
                                                 tissue.colors["Leaf"], leaf.traits.mantel$statistic,
                                                 tissue.colors["Flower"], flower.traits.mantel$statistic,
                                                 leaf.traits.mantel$signif)))

# ---- 7. Build and save the figure -----------------------------------------

ggplot(pairwise_distances.df, aes(x = distance, y = metabolite_distance)) +
  geom_hex(bins = 65, color = NA, alpha = 0.5) +
  geom_smooth(aes(color = tissue),
              method = "gam",
              formula = y ~ s(x, k = 3),
              linewidth = 1.5,
              se = TRUE) +
  ggtext::geom_richtext(data = mantel_labels.df,
                        aes(label = label),
                        hjust = 1,
                        vjust = -0.4,
                        fill = NA,
                        label.color = NA,
                        size = 4.5) +
  facet_wrap(~ distance_type, scales = "free") +
  scale_color_manual(values = tissue.colors) +
  scale_fill_viridis_c(name = "# pairwise\ncomparisons", trans = "log", labels = scales::label_number(accuracy = 1), breaks = c(1, 10, 100, 500, 2500)) +
  labs(x = "Pairwise distance",
       y = "Pairwise feature distance",
       color = NULL) +
  theme_pubr(base_size = 14) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  guides(fill = guide_colorbar(label.theme = element_text(angle = 45, hjust = 1, vjust = 1)))

ggsave("figures/figS6.png", height = 7, width = 10)
