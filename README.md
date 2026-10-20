# Data and code for the article titled "Human-mediated community assembly reshapes urban phytochemical diversity"
Authors: Kilian Perrelet, Gregory Röder, Sebastian Ruile, Sergio Rasmann, Marco Moretti
Corresponding author: Kilian Perrelet, Zuercherstrasse, Birmensdorf, Zürich, Switzerland, kilian.perrelet@gmail.com

This repository contains the raw data, processed data, and R scripts used to characterize the plant communities of 60 urban sites and the leaf metabolome and flower volatile profiles of the plant species found there.



### Workflow summary

1. Run "prepare_metabolomics_df.R".
2. Run "phylo_funct_dist.R".
3. Run "fig2.R"–"fig5.R" and "figS2.R"–"figS10.R".



### Raw data (data/01_raw/)

covariates.csv — Site-level covariates, one row per site. See the data dictionary below.

flora_mat.csv — Site x species plant occurrence matrix. 
Rows = sites, columns = plant species. Values are binary: 1 = species recorded at the site, 0 = not recorded.

flora_taxo.csv — Taxonomy table for all plant species in the study.
Pne row per species, with columns "class", "order", "family", "genus", "species".

flora_traits.csv — Functional traits table. 
See the data dictionary below.

leaf_metabolom.csv — Leaf metabolomics feature table. 
Rows = detected features, columns = feature descriptors followed by one peak-area column per injection.
- "row.ID", "row.m.z", "row.retention.time" — feature ID, mass-to-charge ratio, and retention time (min).
- "<FileName>.mzML.Peak.area" — integrated peak area of each feature in each injection. Column names correspond to the raw-data file names (prefixed with "X" by R). Injections include plant samples ("Sample_<n>"), quality-control injections ("QC"), blanks ("Blank", "B1"–"B10"), and extraction replicates ("ExR"). Only the plant samples are linked to a species via "leaf_metadata.csv".

leaf_metadata.csv — Sample metadata for the leaf metabolomics injections.
One row per plant sample.
- "FileName" — raw-data file name, which matches the column names in "leaf_metabolom.csv" (without the "X" prefix and ".Peak.area" suffix).
- "Sample_name_researcher" — sample number assigned during sample preparation.
- "SampleType" — type of injection (all "Sample").
- "Extraction.date", "Measurment.date" — dates of extraction and LC-MS measurement (dd.mm.yyyy; 5–12 May 2025).
- "Plants" — plant species sampled.
- "PlantID" — unique sample ID.

flower_vocs.csv — Floral volatile organic compound (VOC) feature table. 
Rows = detected features, columns = "row ID", "row m/z", "row retention time", followed by one "<PlantID> peak area" column per sample. 

flower_metadata.csv — Sample metadata for the flower VOC samples. 
One row per sample.
- "PlantID" — sample ID, which matches the column prefixes in "flower_vocs.csv".
- "Plants" — plant species sampled.
- "flower_vocs.csv" contains 10 more sample columns than listed here, which are not linked to a species.

# Data dictionary — "covariates.csv"

Variable | Description 
ugs | Site ID (matches the rows of "flora_mat.csv") 
ugs_type | Urban green space type: "Allotment", "Garden", "Park", "Real Estate", or "Ruderal" (12 sites each) 
human_index | Index of human influence at the site, ranging from 0 to 1 

# Data dictionary — "flora_traits.csv"

Variable | Description
"species" | Plant species (matches "flora_taxo.csv" and the columns of "flora_mat.csv")
"Growth.form" | "herb", "shrub", "tree", or "liana"
"Plant.height.vegetative..m." | Vegetative plant height (m)
"Life.span" | "annual", "biennial", or "perennial"
"Flower.color" | Dominant flower colour ("yellow", "white", "pink/purple/violet", "blue", "green", "red", "orange", "brown/beige/gray")
"Flower.bilateral" | Flower symmetry: "radial" or "bilateral"
"Leaf.greenness" | "green" or "non-green" leaves
"Leaf.elliptic" | Leaf shape: "elliptic", "oblong", or "coniferous"
"Leaf.compoundness" | "simple" or "compound" leaves
"Leaf.deciduous" | "deciduous" or "evergreen"
"Fall.Conspicuous" | Whether autumn foliage color is conspicuous
"is_alien" | Non-native status (1 = non-native, 0 = native)
"SLA..mm2.mg." | Specific leaf area (mm2 mg-1)
"Seed.dry.mass..mg." | Seed dry mass (mg)
"Leaf.water.content.per.leaf.dry.mass..g.W..g.DM." | Leaf water content per unit leaf dry mass (g g-1)
"LDMC..g.g." | Leaf dry matter content (g g-1)
"Winter.hardiness" | Winter hardiness (minimum tolerated temperature, ≤ 0)
"Temperature", "Light", "Moisture", "Nutrients", "Humus", "Aeration" | Ecological indicator values describing the species' habitat preferences
"CRS.life.strategies..competitive.", "..stress_tolerant.", "..ruderal." | Proportion (0–1) of Grime's competitor, stress-tolerator, and ruderal life strategies
"Flowering.duration" | Flowering duration (months, 1–12)



### Processed data ("data/02_processed/")

Leaf files are built from "leaf_metabolom.csv" and "leaf_metadata.csv". 
Flower files are built from "flower_vocs.csv" and "flower_metadata.csv". 
In all files, metabolite features are identified by "row.ID", the feature ID from the raw feature tables.

leaf_mtb_df.csv — Leaf metabolomic features by species. 
Rows = features, columns = "row.ID", "row.m.z" (mass-to-charge ratio), "row.retention.time" (min), followed by one peak-area column per plant species. 

leaf_mtb_mat.csv — The same leaf data as a species x feature matrix (wide format). 
Rows = species, columns = features.

flower_mtb_df.csv — Floral VOC features by species. 
Rows = features, columns = "row.ID", "row.m.z", "row.retention.time", followed by one peak-area column per plant species.

flower_mtb_mat.csv — The same flower data as a species x feature matrix (wide format). 
Rows = species, columns = features.

leaf_site_mtb_mat.csv — Site x leaf-feature matrix. 
Rows = site, columns = features.
Each value is the sum, across the species recorded at the site ("flora_mat.csv"), of that species' peak area for the feature. Only species with leaf metabolomic data contribute.

flower_site_mtb_mat.csv — Site x flower-feature matrix.
Rows = site, columns = features.

phyto_funct_dist.RData — An R workspace with 16 "dist" objects. Each holds the pairwise distances among the 60 sites. Object names follow the pattern "<organ>.<trait set>.dist.<metric>":

Component | Values
organ | "leaf", "flower" (the species with metabolomic data for that organ)
trait set | "phylo" (phylogenetic distance), "all" (functional distance, all traits), "human" (functional distance, human-associated traits), "eco" (functional distance, broader ecological traits)
metric | "comdist" (mean pairwise distance between all species of two communities), "comdistnt" (mean distance between each species and its nearest neighbor in the other community)



### Analysis scripts

prepare_metabolomics_df.R — Builds all processed chemistry tables from the raw leaf and flower data.
Builds the species x feature matrix and site x feature matrix by summing the species profiles over the species present at each site.

phylo_funct_dist.R — Builds the pairwise site dissimilarity objects saved in "phyto_funct_dist.RData", separately for leaf and flower species sets.
Computes site-by-site MPD- and MNTD-type distances from the phylogenetic distances and from Gower functional distances (all, human-associated and ecological traits).

fig2.R — Procrustes comparison of ordination spaces. 

fig3.R — Relationships between site-pair chemical dissimilarity (Jaccard) and functional or phylogenetic MPD and MNTD distances, for leaf and flower. 

fig4.R — Partial least squares regression of site chemistry on community-weighted traits. 

fig5.R — Chemical feature richness along the human-influence gradient. 

figS2.R — Procrustes correlation tests between all pairs of the four ordination spaces: human traits, ecological traits, phylogeny and chemistry.

figS3.R — Feature richness vs plant richness per site type for leaf and flower.

figS4.R — Feature richness by site type.

figS5.R — Variation partitioning of site chemical dissimilarity among human traits, ecological traits and phylogeny, using the first five PCoA axes of each. 

figS6.R — Species-level relationships. Plots pairwise chemical (Jaccard) distance between species against phylogenetic and against Gower functional distance (all traits).

figS7.R — The same diversity metrics as "fig5.R", plus plant richness, plotted against the human-influence index: feature, plant, trait and phylogenetic richness. 

figS8.R — Within- vs between-species variation in flower VOC profiles, using species with more than one replicate.

figS9.R — PCoA of plant-community composition (Jaccard).

figS10.R — Cross-validation of the PLS models from "fig4.R". Fits five components and plots total Q² per component (10-fold, 100 repeats), to justify using two components.
