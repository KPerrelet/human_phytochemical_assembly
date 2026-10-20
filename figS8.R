library(tidyverse)
library(vegan)
library(ggpubr)
library(ggrepel)

set.seed(123) # for reproducibility of the PCoA (cmdscale add-correction) and PERMANOVA below

# ---- 1. Load and harmonize input data --------------------------------------

species.mat <- read.csv("01_raw/flora_mat.csv", check.names = FALSE) # site x plant-species matrix, used only to check which species names are valid

# per-sample metadata: links each collected VOC sample (PlantID) to the plant species it was taken from
flora.meta <- read.csv("01_raw/flower_metadata.csv")

voc.df <- read.csv("01_raw/flower_vocs.csv")
voc.df <- voc.df %>% 
  filter(row.ID != 1421) %>% # Remove naphthalene
  dplyr::select(1:3, paste0("X", sprintf("%02d", flora.meta$PlantID), ".peak.area"))

# relabel sample columns from PlantID to species name
colnames(voc.df)[-c(1:3)] <- flora.meta$Plants
colnames(voc.df) <- make.names(colnames(voc.df), unique = TRUE)

# restrict to species that have more than one replicate sample
keep_cols <- names(voc.df)[-(1:3)] %>%
  .[sub("\\.[0-9]+$", "", .) %in%
      names(table(sub("\\.[0-9]+$", "", .)))[table(sub("\\.[0-9]+$", "", .)) > 1]]

col_order <- tibble(col = keep_cols) %>%
  mutate(species = str_remove(col, "\\.[0-9]+$"),
         rep = as.integer(str_extract(col, "[0-9]+$")),
         rep = replace_na(rep, 0)) %>%
  arrange(species, rep) %>%
  pull(col)

voc.replicated.mat <- voc.df %>% select(all_of(col_order)) # feature x sample matrix, replicated species only

# ---- 2. PCoA ordination of VOC samples, coloured by species (panel A) ------

voc.dist <- vegdist(t(voc.replicated.mat), method = "jaccard", binary = TRUE)

sample.info.df <- tibble(sample = rownames(as.matrix(voc.dist)),
                         species = str_remove(sample, "\\.[0-9]+$")) %>%
  mutate(species = str_replace(species, "\\.", " "))

# classical (metric) multidimensional scaling == PCoA; add = TRUE applies a
# correction (Cailliez) so all eigenvalues are usable even if some are negative
voc.pcoa <- cmdscale(voc.dist, k = 2, eig = TRUE, add = TRUE)

pcoa.scores.df <- as.data.frame(voc.pcoa$points) %>%
  rownames_to_column("sample") %>%
  rename(PCoA1 = V1, PCoA2 = V2) %>%
  mutate(species = str_remove(sample, "\\.[0-9]+$"),
         species = str_replace(species, "\\.", " "))

pcoa.var.pct <- round(voc.pcoa$eig[1:2] / sum(voc.pcoa$eig[voc.pcoa$eig > 0]) * 100, 1)

hulls.df <- pcoa.scores.df %>%
  group_by(species) %>%
  slice(chull(PCoA1, PCoA2)) # convex hull points per species, for shaded polygons

centroids.df <- pcoa.scores.df %>%
  group_by(species) %>%
  summarise(across(c(PCoA1, PCoA2), mean))

pcoa.plot <- ggplot(pcoa.scores.df, aes(PCoA1, PCoA2, color = species)) +
  geom_point(size = 2, alpha = 0.3) +
  geom_polygon(data = hulls.df, aes(fill = species, group = species), alpha = 0.15, linewidth = 0.25) +
  geom_text_repel(data = centroids.df, aes(label = species), size = 4,
                  min.segment.length = 0.1, box.padding = 0.7, point.padding = 0.6) +
  coord_equal() +
  labs(x = paste0("PCoA1 (", pcoa.var.pct[1], "%)"),
       y = paste0("PCoA2 (", pcoa.var.pct[2], "%)")) +
  theme_pubr(base_size = 14) +
  theme(legend.position = "none")


# ---- 3. Within- vs. between-species VOC-profile distances (panel B) --------

# distance of each sample to its own species' group centroid (multivariate dispersion)
betadisper.fit <- betadisper(voc.dist, group = sample.info.df$species)

within.dist.df <- tibble(sample = names(betadisper.fit$distances),
                         distance = betadisper.fit$distances,
                         comparison = "Within\nspecies")

# for each sample, the mean distance to every OTHER species (averaged first within each other species, then across those species means)
sample.species <- sample.info.df$species
names(sample.species) <- sample.info.df$sample

between.dist.df <- lapply(names(sample.species), function(samp) {
  focal_species <- sample.species[samp]
  other_samples <- names(sample.species)[sample.species != focal_species]
  other_distances <- tapply(other_samples, sample.species[other_samples],
                            function(x) mean(as.matrix(voc.dist)[samp, x]))
  tibble(sample = samp, distance = mean(other_distances), comparison = "Between\nspecies")
}) %>%
  bind_rows()

distance.variation.df <- bind_rows(within.dist.df, between.dist.df)

beta.variation.plot <- ggplot(distance.variation.df, aes(x = comparison, y = distance, fill = comparison)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.5) +
  labs(x = NULL, y = "Jaccard distance") +
  theme_pubclean(base_size = 14) +
  theme(legend.position = "none")

# PERMANOVA: tests whether species identity explains a significant share of VOC compositional variance, complementing panel B
voc.permanova <- adonis2(voc.dist ~ species, data = sample.info.df, permutations = 999)


# ---- 4. Within- vs. between-species VOC richness (panels C & D) ------------

richness.df <- voc.replicated.mat %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(richness = vegan::specnumber(select(., -sample)),
         species = str_remove(sample, "\\.[0-9]+$"),
         species = str_replace(species, "\\.", " ")) %>%
  select(sample, species, richness)

richness.plot <- ggplot(richness.df, aes(x = species, y = richness, group = sample, fill = species)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  labs(x = NULL, y = "VOC richness") +
  theme_pubclean(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(l = 20),
        legend.position = "none")

richness.by.species.df <- richness.df %>%
  group_by(species) %>%
  summarise(mean_richness = mean(richness), .groups = "drop")

# within-species: each sample's absolute deviation from its own species' mean richness
richness.within.df <- richness.df %>%
  left_join(richness.by.species.df, by = "species") %>%
  mutate(diff = abs(richness - mean_richness), comparison = "Within\nspecies") %>%
  select(sample, species, diff, comparison)

# between-species: each sample's mean absolute difference from every OTHER species' mean richness
richness.between.df <- richness.df %>%
  rowwise() %>%
  mutate(diff = mean(abs(richness - richness.by.species.df$mean_richness[richness.by.species.df$species != species])),
         comparison = "Between\nspecies") %>%
  ungroup() %>%
  select(sample, species, diff, comparison)

richness.variation.df <- bind_rows(richness.within.df, richness.between.df)

richness.variation.plot <- ggplot(richness.variation.df, aes(x = comparison, y = diff, fill = comparison)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.5) +
  labs(x = NULL, y = "Difference in VOC richness") +
  theme_pubclean(base_size = 14) +
  theme(legend.position = "none")


# ---- 5. Assemble and save the figure ---------------------------------------

ggarrange(pcoa.plot, beta.variation.plot, richness.plot, richness.variation.plot,
          labels = c("A", "B", "C", "D"),
          widths = c(1, 0.5, 1, 0.5),
          align = "v")

ggsave("figures/figS8.png", height = 10, width = 10)
