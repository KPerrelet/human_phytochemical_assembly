library(tidyverse)

# ---- 1. Load the shared plant-community matrix -----------------------------

# site x plant-species matrix
species.mat <- read.csv("01_raw/flora_mat.csv", check.names = FALSE) %>%
  column_to_rownames("ugs")

################################################################################
# ---- 2. LEAF metabolomics processing ---------------------------------------
################################################################################

# --- 2a. Load and filter sample metadata (keep only real samples, not blanks) ---
leaf.meta <- read.csv("01_raw/leaf_metadata.csv")
leaf.meta <- leaf.meta[leaf.meta$SampleType == "Sample", ]

# --- 2b. Load the raw peak table and keep only columns for retained samples ---
# rows = detected compounds (with 3 metadata columns: row ID, m/z, retention time, etc.); columns 4+ = one peak-area column per sample
leaf.peaks.raw <- read.csv("01_raw/leaf_metabolom.csv")
leaf.peaks.raw <- leaf.peaks.raw %>%
  dplyr::select(1:3, any_of(paste0("X", leaf.meta$FileName, ".Peak.area")))

# --- 2c. Rename sample columns to species names, then average replicate
# samples of the same species down to one column per species ---
leaf.meta_cols <- leaf.peaks.raw[, 1:3] # compound metadata, set aside to reattach after aggregation

colnames(leaf.peaks.raw)[4:ncol(leaf.peaks.raw)] <- leaf.meta$Plants
colnames(leaf.peaks.raw) <- make.names(colnames(leaf.peaks.raw), unique = TRUE)

leaf.compound.df <- leaf.peaks.raw[, -c(1:3)] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  mutate(species = gsub("\\.", " ", species),
         species = gsub(" [0-9]+$", "", species)) %>% # strip the ".1", ".2"... replicate suffixes make.names() added
  group_by(species) %>%
  summarise(across(everything(), ~ mean(as.numeric(.), na.rm = TRUE))) %>% # average across replicate samples
  column_to_rownames("species") %>%
  t() %>%
  as.data.frame()

leaf.compound.df <- bind_cols(leaf.meta_cols, leaf.compound.df)

# fix species-name artefacts introduced by make.names() (dots replacing spaces/punctuation that needs to be restored to match flora.mat's species names)
colnames(leaf.compound.df) <- gsub("sp ", "sp.", colnames(leaf.compound.df))
colnames(leaf.compound.df)[colnames(leaf.compound.df) == "Dryopteris filix mas"] <- "Dryopteris filix-mas"
colnames(leaf.compound.df)[colnames(leaf.compound.df) == "Malus domestica  Idared "] <- "Malus domestica \"Idared\""
colnames(leaf.compound.df)[colnames(leaf.compound.df) == "Ribes uva crispa"] <- "Ribes uva-crispa"
colnames(leaf.compound.df)[colnames(leaf.compound.df) == "Symphyotrichum novae angliae"] <- "Symphyotrichum novae-angliae"
colnames(leaf.compound.df)[colnames(leaf.compound.df) == "Symphyotrichum novi belgii"] <- "Symphyotrichum novi-belgii"
colnames(leaf.compound.df)[colnames(leaf.compound.df) == "Vitex agnus castus"] <- "Vitex agnus-castus"

# --- 2d. Species x compound matrix: transpose, coerce to numeric, drop compounds with zero signal across all species, and label compound columns by their original row ID ---
leaf.species.mat <- leaf.compound.df[, -c(1:3)] %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  replace(is.na(.), 0) %>%
  select(where(~ sum(.x) > 0)) %>%
  mutate(species = rownames(.)) %>%
  relocate(species)
colnames(leaf.species.mat)[-1] <- leaf.compound.df$row.ID[rowSums(leaf.compound.df[, -c(1:3)]) > 0]

# --- 2e. Site x compound matrix, by combining plant-community matrix with species-level chemistry: for each site, sum each species' chemical profile weighted by that species' abundance at the site ---
leaf.shared_species <- intersect(colnames(species.mat), rownames(leaf.species.mat))
leaf.site.mat <- as.data.frame(as.matrix(species.mat[, leaf.shared_species]) %*% as.matrix(leaf.species.mat[leaf.shared_species, -1]))
leaf.site.mat$ugs <- rownames(leaf.site.mat)
leaf.site.mat <- leaf.site.mat %>% relocate(ugs)

# --- 2f. Save leaf outputs ---
write.csv(leaf.compound.df, "02_processed/leaf_mtb_df.csv", row.names = FALSE)
write.csv(leaf.species.mat, "02_processed/leaf_mtb_mat.csv", row.names = FALSE)
write.csv(leaf.site.mat,    "02_processed/leaf_site_mtb_mat.csv", row.names = FALSE)

################################################################################
# ---- 3. FLOWER VOC processing -----------------------------------------------
################################################################################

# --- 3a. Load sample metadata (no Blank/Sample split needed for this dataset) ---
flower.meta <- read.csv("01_raw/flower_metadata.csv")

# --- 3b. Load the raw peak table and keep only columns for retained samples ---
flower.peaks.raw <- read.csv("01_raw/flower_vocs.csv")
flower.peaks.raw <- flower.peaks.raw %>%
  filter(row.ID != 1421) %>% # remove naphthalene (a known contaminant/artefact peak)
  dplyr::select(1:3, paste0("X", sprintf("%02d", flower.meta$PlantID), ".peak.area"))

# --- 3c. Rename sample columns to species names, then average replicate samples of the same species down to one column per species ---
flower.meta_cols <- flower.peaks.raw[, 1:3]

colnames(flower.peaks.raw)[-c(1:3)] <- flower.meta$Plants
colnames(flower.peaks.raw) <- make.names(colnames(flower.peaks.raw), unique = TRUE)

flower.compound.df <- flower.peaks.raw[, -c(1:3)] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  mutate(species = gsub("\\.", " ", species),
         species = gsub(" [0-9]+$", "", species)) %>%
  group_by(species) %>%
  summarise(across(everything(), ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  column_to_rownames("species") %>%
  t() %>%
  as.data.frame()

flower.compound.df <- bind_cols(flower.meta_cols, flower.compound.df)

# fix species-name artefacts introduced by make.names()
colnames(flower.compound.df)[colnames(flower.compound.df) == "Campanula glomerata subsp  glomerata"] <- "Campanula glomerata subsp. glomerata"
colnames(flower.compound.df)[colnames(flower.compound.df) == "Rosa  Little Mozart "]                  <- "Rosa \"Little Mozart\""
colnames(flower.compound.df)[colnames(flower.compound.df) == "Symphyotrichum novae angliae"]          <- "Symphyotrichum novae-angliae"

# --- 3e. Species x compound matrix ---
flower.species.mat <- flower.compound.df[, -c(1:3)] %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  replace(is.na(.), 0) %>%
  select(where(~ sum(.x) > 0)) %>%
  mutate(species = rownames(.)) %>%
  relocate(species)
colnames(flower.species.mat)[-1] <- flower.compound.df$row.ID[rowSums(flower.compound.df[, -c(1:3)]) > 0]

# --- 3f. Site x compound matrix ---
flower.shared_species <- intersect(colnames(species.mat), rownames(flower.species.mat))
flower.site.mat <- as.data.frame(as.matrix(species.mat[, flower.shared_species]) %*% as.matrix(flower.species.mat[flower.shared_species, -1]))
flower.site.mat$ugs <- rownames(flower.site.mat)
flower.site.mat <- flower.site.mat %>% relocate(ugs)

# --- 3g. Save flower outputs ---
write.csv(flower.compound.df, "02_processed/flower_mtb_df.csv", row.names = FALSE)
write.csv(flower.species.mat, "02_processed/flower_mtb_mat.csv", row.names = FALSE)
write.csv(flower.site.mat,    "02_processed/flower_site_mtb_mat.csv", row.names = FALSE)
