
# Dissertation code workflow for Ecological and Environmental Science
# Eve Boswood
# Finished on 03.04.2026


# Libraries 
library(tidyverse)
library(readxl)
library(boot)
library(vegan)
library(betapart)
library(tibble)
library(patchwork)
library(FD)
library(lme4)
library(lmerTest)



# APPENDIX A1.1 data prep
# From raw subcell records to aspect-level community objects

set.seed(123)

# Import raw species data from excel
data_species <- read_excel("dissertation_data_species.xlsx")

# Clean
data_species <- data_species %>%
  mutate(
    species_present = na_if(species_present, "N/A"),
    species_present = str_squish(species_present)
  )

# Site metadata
site_meta <- data_species %>%
  distinct(site_id) %>%
  mutate(
    cove = str_extract(site_id, "C\\d+"),
    side = case_when(
      str_detect(site_id, "NF") ~ "north-facing",
      str_detect(site_id, "SF") ~ "south-facing",
      TRUE ~ NA_character_
    )
  )

# long format
species_data_long <- data_species %>%
  filter(!is.na(species_present), species_present != "") %>%
  separate_rows(species_present, sep = ",") %>%
  mutate(species_present = str_squish(species_present)) %>%
  filter(species_present != "") %>%
  left_join(site_meta, by = "site_id")

# Presence/absence at subcell level
pa_subcell <- species_data_long %>%
  distinct(site_id, cove, side, quadrat_cell, species_present) %>%
  mutate(present = 1L)

# Sampled subcells per site 
site_denoms <- data_species %>%
  distinct(site_id, quadrat_cell) %>%
  count(site_id, name = "n_subcells") %>%
  left_join(site_meta, by = "site_id")

# Species-level cover per site
species_cover_site <- pa_subcell %>%
  count(site_id, cove, side, species_present, name = "n_cells_present") %>%
  left_join(site_denoms, by = c("site_id", "cove", "side")) %>%
  mutate(
    cover_prop = n_cells_present / n_subcells,
    cover_pct  = cover_prop * 100
  )

# Species richness per site
site_richness <- species_data_long %>%
  distinct(site_id, cove, side, species_present) %>%
  count(site_id, cove, side, name = "richness")

# Macroalgal cover per site
site_total_cover <- species_cover_site %>%
  group_by(site_id, cove, side) %>%
  summarise(
    total_cover = sum(cover_pct, na.rm = TRUE),
    .groups = "drop"
  )

# Wide community matrix for proportional cover
community_mat <- species_cover_site %>%
  select(site_id, species_present, cover_prop) %>%
  pivot_wider(
    names_from  = species_present,
    values_from = cover_prop,
    values_fill = 0
  ) %>%
  left_join(site_meta, by = "site_id")

# Community matrix for presence-absence
community_pa <- community_mat %>%
  mutate(
    across(
      .cols = -c(site_id, cove, side),
      .fns  = ~ as.integer(.x > 0)
    )
  )

# summary table
site_summary <- site_meta %>%
  left_join(site_richness, by = c("site_id", "cove", "side")) %>%
  left_join(site_total_cover, by = c("site_id", "cove", "side")) %>%
  left_join(site_denoms, by = c("site_id", "cove", "side"))













# APPENDIX A1.2. Environmental differences between aspects (H1)

# environmental summary table
env_data <- read_excel("env_data.xlsx")

# paired cove-level structure
env <- env_data %>%
  mutate(
    cove = str_extract(site_id, "^C[0-9]+"),
    side = str_extract(site_id, "(SF|NF)$"),
    side = factor(side, levels = c("NF", "SF"))
  )

env_pairs <- env %>%
  select(cove, side, air_temp_mean, lux_mean, relative_wave_exposure) %>%
  pivot_wider(
    names_from = side,
    values_from = c(air_temp_mean, lux_mean, relative_wave_exposure),
    names_sep = "__"
  ) %>%
  mutate(
    d_temp = air_temp_mean__SF - air_temp_mean__NF,
    d_lux  = lux_mean__SF - lux_mean__NF,
    d_wave = relative_wave_exposure__SF - relative_wave_exposure__NF
  )

# Exact paired sign-flip test (because overdispersed)
exact_signflip_p <- function(diffs) {
  diffs <- diffs[!is.na(diffs)]
  n <- length(diffs)
  
  all_signs <- expand.grid(rep(list(c(-1, 1)), n))
  perm_means <- apply(all_signs, 1, function(s) mean(as.numeric(s) * diffs))
  obs <- mean(diffs)
  
  mean(abs(perm_means) >= abs(obs))
}

p_temp <- exact_signflip_p(env_pairs$d_temp)
p_lux  <- exact_signflip_p(env_pairs$d_lux)
p_wave <- exact_signflip_p(env_pairs$d_wave)

# Bootstrap 95% confidence intervals
boot_mean_diff <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

boot_ci_mean <- function(diffs, R = 10000) {
  diffs <- diffs[!is.na(diffs)]
  
  b <- boot(diffs, statistic = boot_mean_diff, R = R)
  ci <- boot.ci(b, type = "perc")
  
  tibble(
    mean_diff = mean(diffs),
    ci_low = ci$percent[4],
    ci_high = ci$percent[5]
  )
}

temp_summary <- boot_ci_mean(env_pairs$d_temp)
lux_summary  <- boot_ci_mean(env_pairs$d_lux)
wave_summary <- boot_ci_mean(env_pairs$d_wave)

# Summary 
h1_summary <- bind_rows(
  temp_summary %>% mutate(metric = "Mean air temperature (°C)", p_value = p_temp),
  lux_summary  %>% mutate(metric = "Mean light intensity (lux)", p_value = p_lux),
  wave_summary %>% mutate(metric = "Relative wave exposure", p_value = p_wave)
) %>%
  select(metric, mean_diff, ci_low, ci_high, p_value)

h1_summary














# APPENDIX A1.3. Univariate community responses (H2)
# Species richness and total macroalgal cover


# paired cove-level dataset
h2_pairs <- site_summary %>%
  mutate(
    side_code = str_extract(site_id, "(NF|SF)$")
  ) %>%
  select(cove, side_code, richness, total_cover) %>%
  pivot_wider(
    id_cols = cove,
    names_from = side_code,
    values_from = c(richness, total_cover),
    names_sep = "__"
  ) %>%
  mutate(
    d_richness = richness__SF - richness__NF,
    d_cover    = total_cover__SF - total_cover__NF
  )

# Paired t-tests (not overdispersed)
t_richness <- t.test(
  h2_pairs$richness__SF,
  h2_pairs$richness__NF,
  paired = TRUE
)

t_cover <- t.test(
  h2_pairs$total_cover__SF,
  h2_pairs$total_cover__NF,
  paired = TRUE
)

# Bootstrap 95% confidence intervals
boot_mean_diff <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

boot_ci_mean <- function(diffs, R = 10000) {
  diffs <- diffs[!is.na(diffs)]
  
  b <- boot(diffs, statistic = boot_mean_diff, R = R)
  ci <- boot.ci(b, type = "perc")
  
  tibble(
    mean_diff = mean(diffs),
    ci_low = ci$percent[4],
    ci_high = ci$percent[5]
  )
}

cover_summary <- boot_ci_mean(h2_pairs$d_cover)
richness_summary <- boot_ci_mean(h2_pairs$d_richness)

# summary table
h2_summary <- bind_rows(
  cover_summary %>%
    mutate(
      metric = "Total macroalgal cover (summed % cover)",
      test_statistic = unname(t_cover$statistic),
      df = unname(t_cover$parameter),
      p_value = t_cover$p.value
    ),
  richness_summary %>%
    mutate(
      metric = "Species richness",
      test_statistic = unname(t_richness$statistic),
      df = unname(t_richness$parameter),
      p_value = t_richness$p.value
    )
) %>%
  select(metric, mean_diff, ci_low, ci_high, test_statistic, df, p_value)

h2_summary













# APPENDIX A.1.4. Multivariate community structure (H3)
# PERMANOVA, dispersion, NMDS, beta-diversity partitioning

# community matrix and metadata
comm <- community_mat %>%
  column_to_rownames("site_id")

meta_h3 <- comm %>%
  rownames_to_column("site_id") %>%
  select(site_id, cove, side)

comm <- comm %>%
  select(-cove, -side)

# 4th-root transformation, Bray-Curtis dissimilarity
comm_4rt <- comm^(1/4)
dist_h3 <- vegdist(comm_4rt, method = "bray")

# PERMANOVA
adon_h3 <- adonis2(
  dist_h3 ~ side,
  data = meta_h3,
  permutations = 999,
  strata = meta_h3$cove
)

adon_h3

# Homogeneity of dispersion
bd_h3 <- betadisper(dist_h3, group = meta_h3$side)

anova(bd_h3)
permutest(bd_h3, permutations = 999)

# NMDS on Bray-Curtis
set.seed(123)

nmds_h3 <- metaMDS(
  comm_4rt,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  wascores = FALSE
)

nmds_h3$stress

scores_h3 <- as.data.frame(scores(nmds_h3, display = "sites")) %>%
  rownames_to_column("site_id") %>%
  left_join(meta_h3, by = "site_id")

#  NMDS plot data with C2_SF omitted for plotting 
scores_h3_plot <- scores_h3 %>%
  filter(site_id != "C2_SF")

#Presence–absence matrix for beta-diversity 
comm_pa <- decostand(comm, method = "pa")

# north-facing vs south-facing comparisons within coves
pair_beta_by_cove <- map_dfr(unique(meta_h3$cove), function(cv) {
  ids <- meta_h3 %>%
    filter(cove == cv) %>%
    arrange(side) %>%
    pull(site_id)
  
  mat_cv <- comm_pa[ids, , drop = FALSE]
  
  beta_pair <- betapart::beta.pair(mat_cv, index.family = "jaccard")
  
  tibble(
    cove = cv,
    turnover = as.matrix(beta_pair$beta.jtu)[1, 2],
    nestedness = as.matrix(beta_pair$beta.jne)[1, 2],
    total_beta = as.matrix(beta_pair$beta.jac)[1, 2]
  )
})

pair_beta_by_cove

beta_summary <- pair_beta_by_cove %>%
  summarise(
    mean_turnover = mean(turnover),
    mean_nestedness = mean(nestedness),
    mean_total_beta = mean(total_beta),
    turnover_prop = mean(turnover / total_beta),
    nestedness_prop = mean(nestedness / total_beta)
  )

beta_summary

# Species-pool asymmetry
nf_species <- meta_h3 %>%
  filter(side == "north-facing") %>%
  pull(site_id) %>%
  {\(ids) colnames(comm_pa)[colSums(comm_pa[ids, , drop = FALSE]) > 0]}()

sf_species <- meta_h3 %>%
  filter(side == "south-facing") %>%
  pull(site_id) %>%
  {\(ids) colnames(comm_pa)[colSums(comm_pa[ids, , drop = FALSE]) > 0]}()

species_pool_summary <- tibble(
  shared_species = length(intersect(nf_species, sf_species)),
  nf_unique_species = length(setdiff(nf_species, sf_species)),
  sf_unique_species = length(setdiff(sf_species, nf_species))
)

species_pool_summary










# APPENDIX A1.5. Environmental drivers of community structure (H4)
# dbRDA, PERMANOVA, richness GLM


# community matrix and metadata
comm_h4 <- community_mat %>%
  column_to_rownames("site_id")

meta_h4 <- comm_h4 %>%
  rownames_to_column("site_id") %>%
  select(site_id, cove, side)

comm_h4 <- comm_h4 %>%
  select(-cove, -side)

# fourth-root transform community data
comm_4rt <- comm_h4^(1/4)
dist_h4 <- vegdist(comm_4rt, method = "bray")

# Environmental predictors
env_h4 <- env_data %>%
  mutate(
    cove = str_extract(site_id, "^C[0-9]+"),
    side = str_extract(site_id, "(SF|NF)$"),
    side = factor(side, levels = c("NF", "SF")),
    cove = factor(cove),
    temp_z = as.numeric(scale(air_temp_mean)),
    lux_z  = as.numeric(scale(lux_mean)),
    wave_z = as.numeric(scale(relative_wave_exposure))
  ) %>%
  filter(site_id %in% rownames(comm_4rt)) %>%
  arrange(match(site_id, rownames(comm_4rt)))

stopifnot(all(env_h4$site_id == rownames(comm_4rt)))

# check collinearity 
cor(env_h4 %>% select(temp_z, lux_z, wave_z), use = "pairwise.complete.obs")

# dbRDA
dbrda_mean <- capscale(
  comm_4rt ~ temp_z + lux_z + wave_z,
  data = env_h4,
  distance = "bray"
)

anova(dbrda_mean, permutations = 999)
anova(dbrda_mean, by = "margin", permutations = 999)

# PERMANOVA side-only 
adon_side <- adonis2(
  dist_h4 ~ side,
  data = env_h4,
  permutations = 999,
  strata = env_h4$cove
)

adon_side

# PERMANOVA full  model
adon_full <- adonis2(
  dist_h4 ~ side + temp_z + lux_z + wave_z,
  data = env_h4,
  permutations = 999,
  strata = env_h4$cove,
  by = "margin"
)

adon_full

# Species richness GLM
rich_h4 <- site_summary %>%
  left_join(
    env_h4 %>%
      select(site_id, temp_z, lux_z, wave_z),
    by = "site_id"
  ) %>%
  mutate(
    side = factor(side, levels = c("north-facing", "south-facing")),
    cove = factor(cove)
  )

m_rich_env <- glm(
  richness ~ cove + side + temp_z + lux_z + wave_z,
  data = rich_h4,
  family = poisson(link = "log")
)

summary(m_rich_env)
drop1(m_rich_env, test = "Chisq")










# APPENDIX A1.6. Grazer effects (H5)



# quadrat identity from subcell 
# assumption... each site contains 10 quadrats x 100 subcells
cell_lookup <- data_species %>%
  distinct(site_id, quadrat_cell) %>%
  arrange(site_id, quadrat_cell) %>%
  group_by(site_id) %>%
  mutate(
    cell_rank = row_number(),
    quadrat_id = ((cell_rank - 1) %/% 100) + 1
  ) %>%
  ungroup()

# Limpet counts per subcell
# Replace NA's with 0 
data_limp <- data_species %>%
  mutate(
    limpet_common = replace_na(limpet_common, 0),
    limpet_china  = replace_na(limpet_china, 0),
    limpet_total  = limpet_common + limpet_china
  ) %>%
  left_join(cell_lookup, by = c("site_id", "quadrat_cell"))

# Quadrat-level limpet totals
limp_quadrat <- data_limp %>%
  group_by(site_id, quadrat_id) %>%
  summarise(
    limpet_total_quadrat = sum(limpet_total, na.rm = TRUE),
    .groups = "drop"
  )

# site-level: mean limpet abundance per quadrat
limp_site <- limp_quadrat %>%
  group_by(site_id) %>%
  summarise(
    mean_limp = mean(limpet_total_quadrat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cove = str_extract(site_id, "^C\\d+"),
    side = case_when(
      str_detect(site_id, "NF$") ~ "north-facing",
      str_detect(site_id, "SF$") ~ "south-facing",
      TRUE ~ NA_character_
    )
  )

# cove-level 
limp_pairs <- limp_site %>%
  mutate(
    side_code = str_extract(site_id, "(NF|SF)$")
  ) %>%
  select(cove, side_code, mean_limp) %>%
  pivot_wider(
    id_cols = cove,
    names_from = side_code,
    values_from = mean_limp,
    names_glue = "mean_limp__{side_code}"
  ) %>%
  mutate(
    d_limp = mean_limp__SF - mean_limp__NF
  )

# Paired t-test
t_limp <- t.test(
  limp_pairs$mean_limp__SF,
  limp_pairs$mean_limp__NF,
  paired = TRUE
)

# Bootstrap 95% CI
boot_mean_diff <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

boot_ci_mean <- function(diffs, R = 10000) {
  diffs <- diffs[!is.na(diffs)]
  
  b <- boot(diffs, statistic = boot_mean_diff, R = R)
  ci <- boot.ci(b, type = "perc")
  
  tibble(
    mean_diff = mean(diffs),
    ci_low = ci$percent[4],
    ci_high = ci$percent[5]
  )
}

limp_summary <- boot_ci_mean(limp_pairs$d_limp)

h5_pair_summary <- limp_summary %>%
  mutate(
    metric = "Mean limpet abundance per quadrat",
    test_statistic = unname(t_limp$statistic),
    df = unname(t_limp$parameter),
    p_value = t_limp$p.value
  ) %>%
  select(metric, mean_diff, ci_low, ci_high, test_statistic, df, p_value)

h5_pair_summary

# Site-level dataframe
limp_macro <- limp_site %>%
  left_join(
    site_summary %>% select(site_id, richness, total_cover, cove),
    by = c("site_id", "cove")
  ) %>%
  mutate(
    cove = factor(cove)
  )

# Richness ~ limpets + cove (quasi-Poisson)
m_limp_rich <- glm(
  richness ~ mean_limp + cove,
  data = limp_macro,
  family = quasipoisson(link = "log")
)

summary(m_limp_rich)
drop1(m_limp_rich, test = "F")

# Cover ~ limpets + cove (Gaussian)
m_limp_cover <- lm(
  total_cover ~ mean_limp + cove,
  data = limp_macro
)

summary(m_limp_cover)
drop1(m_limp_cover, test = "F")


















# APPENDIX A1.7. Functional traits and diversity (H6)
# Functional richness, functional dispersion, and CWM PCA


set.seed(123)

# import matrix
trait_matrix <- read_excel("quantitative_trait_matrix.xlsx") %>%
  mutate(
    species = str_squish(as.character(species)),
    species = na_if(species, "")
  ) %>%
  filter(!is.na(species)) %>%
  distinct(species, .keep_all = TRUE) %>%
  column_to_rownames("species")

cat_cols <- c(
  "pneumatocysts", "MFGs", "FGs_simple", "FGs",
  "canopy_turf", "canopy/subcanopy/turf"
)

drop_cols <- c(
  "group", "EGs_k_medoids", "FGs_simple", "canopy_turf"
)

exclude_taxa <- c("Encrusting Brown", "Encrusting Red", "Encrusting Red 2")

# categorical traits -> factors
trait_matrix <- trait_matrix %>%
  mutate(
    across(
      any_of(cat_cols),
      ~ as.factor(trimws(as.character(.x)))
    )
  )

# Remove SD/SE columns not needed for analysis
trait_fd <- trait_matrix %>%
  select(
    -matches("(_SD|_SE)$"),
    -any_of(drop_cols)
  )

# aspect-level community matrix
comm_h6 <- community_mat %>%
  column_to_rownames("site_id") %>%
  select(-cove, -side) %>%
  as.data.frame()

comm_h6[] <- lapply(comm_h6, as.numeric)
colnames(comm_h6) <- str_squish(colnames(comm_h6))

# Match taxa with trait matrices
taxa_keep <- intersect(colnames(comm_h6), rownames(trait_fd))
taxa_keep <- setdiff(taxa_keep, exclude_taxa)

comm_fd <- comm_h6[, taxa_keep, drop = FALSE]
trait_fd <- trait_fd[taxa_keep, , drop = FALSE]

# Drop trait columns if all NA or no variation
trait_fd <- trait_fd %>%
  mutate(across(where(is.factor), droplevels))

keep_cols <- names(trait_fd)[vapply(trait_fd, function(x) {
  vals <- x[!is.na(x)]
  length(vals) > 0 && dplyr::n_distinct(vals) > 1
}, logical(1))]

trait_fd <- trait_fd[, keep_cols, drop = FALSE]

# checkk it
stopifnot(length(taxa_keep) > 2)
stopifnot(ncol(trait_fd) > 1)
stopifnot(identical(colnames(comm_fd), rownames(trait_fd)))
stopifnot(all(rowSums(comm_fd, na.rm = TRUE) > 0))

# Functional diversity 
fd_out <- dbFD(
  x = trait_fd,
  a = as.matrix(comm_fd),
  stand.x = TRUE,
  w.abun = TRUE,
  corr = "cailliez",
  calc.FRic = TRUE,
  messages = FALSE
)

# site-level indices and no.trait-bearing taxa
fd_idx <- tibble(
  site_id = rownames(comm_fd),
  FRic = fd_out$FRic,
  FDis = fd_out$FDis,
  n_trait_taxa = rowSums(comm_fd > 0)
) %>%
  left_join(
    community_mat %>% select(site_id, cove, side),
    by = "site_id"
  ) %>%
  mutate(
    side_code = case_when(
      side == "north-facing" ~ "NF",
      side == "south-facing" ~ "SF",
      TRUE ~ NA_character_
    ),
    # FDis is not informative for single-species communities
    FDis = if_else(n_trait_taxa < 2, NA_real_, FDis)
  )

fd_idx

# Paired differences
fric_diff <- fd_idx %>%
  select(cove, side_code, FRic) %>%
  pivot_wider(names_from = side_code, values_from = FRic) %>%
  filter(!is.na(NF) & !is.na(SF)) %>%
  mutate(d_FRic = SF - NF)

fdis_diff <- fd_idx %>%
  select(cove, side_code, FDis, n_trait_taxa) %>%
  pivot_wider(
    names_from = side_code,
    values_from = c(FDis, n_trait_taxa)
  ) %>%
  filter(
    !is.na(FDis_NF) & !is.na(FDis_SF),
    n_trait_taxa_NF >= 2,
    n_trait_taxa_SF >= 2
  ) %>%
  mutate(d_FDis = FDis_SF - FDis_NF)

# Normality?
shapiro.test(fric_diff$d_FRic)
shapiro.test(fdis_diff$d_FDis)

# t-tests
t_fric <- t.test(fric_diff$d_FRic)
t_fdis <- t.test(fdis_diff$d_FDis)

# table
h6_summary <- tibble(
  metric = c("Functional richness (FRic)", "Functional dispersion (FDis)"),
  n_pairs = c(length(fric_diff$d_FRic), length(fdis_diff$d_FDis)),
  mean_diff = c(mean(fric_diff$d_FRic), mean(fdis_diff$d_FDis)),
  ci_low = c(t_fric$conf.int[1], t_fdis$conf.int[1]),
  ci_high = c(t_fric$conf.int[2], t_fdis$conf.int[2]),
  test_statistic = c(unname(t_fric$statistic), unname(t_fdis$statistic)),
  df = c(unname(t_fric$parameter), unname(t_fdis$parameter)),
  p_value = c(t_fric$p.value, t_fdis$p.value)
)

h6_summary


# Community-weighted means 
cwm_df <- as.data.frame(fd_out$CWM) %>%
  rownames_to_column("site_id") %>%
  left_join(
    fd_idx %>% select(site_id, cove, side),
    by = "site_id"
  )

cwm_df

# PCA of CWM 
cwm_num <- cwm_df %>%
  select(where(is.numeric)) %>%
  as.data.frame()

rownames(cwm_num) <- cwm_df$site_id

pca_cwm <- prcomp(cwm_num, scale. = TRUE)

cwm_scores <- as.data.frame(pca_cwm$x[, 1:2]) %>%
  rownames_to_column("site_id") %>%
  left_join(
    cwm_df %>% select(site_id, cove, side),
    by = "site_id"
  )

cwm_scores

# Convex hulls 
cwm_hulls <- cwm_scores %>%
  group_by(side) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

cwm_hulls














# APPENDIX A1.8. Quadrat-scale control analyses
# Does rugosity have an effect on proportional cover and species richness?


# rugosity data
rugosity <- read_excel("rugosity_quadrat.xlsx")

# quadrat identity from subcells
# Again assumes 10 quadrats per site and 100 subcells per quadrat
cell_lookup <- data_species %>%
  distinct(site_id, quadrat_cell) %>%
  arrange(site_id, quadrat_cell) %>%
  group_by(site_id) %>%
  mutate(
    cell_rank = row_number(),
    quadrat = ((cell_rank - 1) %/% 100) + 1
  ) %>%
  ungroup()

# Full list of quadrats 
all_quadrats <- cell_lookup %>%
  distinct(site_id, quadrat)

# Quadrat proportional cover
cover_quadrat <- data_species %>%
  left_join(cell_lookup, by = c("site_id", "quadrat_cell")) %>%
  mutate(
    algae_present = if_else(
      is.na(species_present) | species_present == "" | species_present == "N/A",
      0L, 1L
    )
  ) %>%
  group_by(site_id, quadrat) %>%
  summarise(
    total_cover = mean(algae_present),
    .groups = "drop"
  ) %>%
  right_join(all_quadrats, by = c("site_id", "quadrat")) %>%
  mutate(total_cover = replace_na(total_cover, 0))

# Quadrat-level species richness
richness_quadrat <- pa_subcell %>%
  left_join(cell_lookup, by = c("site_id", "quadrat_cell")) %>%
  distinct(site_id, quadrat, species_present) %>%
  count(site_id, quadrat, name = "richness") %>%
  right_join(all_quadrats, by = c("site_id", "quadrat")) %>%
  mutate(richness = replace_na(richness, 0L))

# joining rugosity onto quadrat responses
rug_total <- cover_quadrat %>%
  left_join(rugosity, by = c("site_id", "quadrat"))

rug_rich <- richness_quadrat %>%
  left_join(rugosity, by = c("site_id", "quadrat"))

# Check check
sum(is.na(rug_total$rugosity_index))
sum(is.na(rug_rich$rugosity_index))
nrow(rug_total)
nrow(rug_rich)

# Gaussian mixed-effects model: proportional cover ~ rugosity
m_rug_cover <- lmer(
  total_cover ~ rugosity_index + (1 | site_id),
  data = rug_total
)

summary(m_rug_cover)

# Poisson mixed-effects model: richness ~ rugosity
m_rug_rich <- glmer(
  richness ~ rugosity_index + (1 | site_id),
  data = rug_rich,
  family = poisson(link = "log")
)

summary(m_rug_rich)

# overdispersion?
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  chisq <- sum(rp^2)
  ratio <- chisq / rdf
  p <- pchisq(chisq, df = rdf, lower.tail = FALSE)
  c(chisq = chisq, ratio = ratio, rdf = rdf, p_value = p)
}

overdisp_fun(m_rug_rich)










# APPENDIX A1.9. Validation of percent-cover estimates
# calibration of Fucus vesiculosus cover and holdfast abundance



# raw data with holdfast column as numeric
data_species <- read_excel(
  "dissertation_data_species.xlsx",
  col_types = c(
    "text",    # site_id
    "numeric", # quadrat_cell
    "numeric", # limpet_common
    "numeric", # limpet_china
    "text",    # species_present
    "numeric", # holdfasts_f_vesiculosus
    "numeric", # holdfasts_m_stellatus
    "numeric", # holdfasts_p_umbilicalis
    "numeric"  # holdfasts_u_intestinalis
  )
)

# Site-level proportional cover and holdfast totals for Fucus vesiculosus
fv_calibration <- data_species %>%
  group_by(site_id) %>%
  summarise(
    cover_prop = mean(grepl("Fucus vesiculosus", species_present)),
    holdfast_total = if (all(is.na(holdfasts_f_vesiculosus))) NA_real_
    else sum(holdfasts_f_vesiculosus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(holdfast_total))

fv_calibration

# Poisson model and overdispersion 
m_fv_pois <- glm(
  holdfast_total ~ cover_prop,
  family = poisson(link = "log"),
  data = fv_calibration
)

summary(m_fv_pois)

dispersion <- sum(residuals(m_fv_pois, type = "pearson")^2) / df.residual(m_fv_pois)
dispersion

# Therefore, negative binomial regression
m_fv_nb <- MASS::glm.nb(
  holdfast_total ~ cover_prop,
  data = fv_calibration
)

summary(m_fv_nb)

# Predicted values -> observed-predicted correlation
fv_calibration <- fv_calibration %>%
  mutate(
    predicted_holdfasts = predict(m_fv_nb, type = "response")
  )

cor(
  fv_calibration$holdfast_total,
  fv_calibration$predicted_holdfasts,
  use = "complete.obs"
)













# APPENDIX A1.10a. Extreme temperature metric
# Cove-specific 95th percentile threshold exceedance


# data from all the temperatures 
temp_all <- read_excel("temp_data.xlsx", col_types = "text") %>%
  mutate(
    temperature_C = as.numeric(temperature_C),
    light_lux = as.numeric(light_lux)
  )

# Make 95th percentile metric
logger_extreme <- temp_all %>%
  mutate(
    cove = if_else(
      str_detect(as.character(cove), "^C"),
      as.character(cove),
      paste0("C", cove)
    ),
    side_code = case_when(
      side %in% c("NF", "SF") ~ as.character(side),
      str_detect(tolower(as.character(side)), "north") ~ "NF",
      str_detect(tolower(as.character(side)), "south") ~ "SF",
      TRUE ~ as.character(side)
    ),
    site_id = paste0(cove, "_", side_code)
  ) %>%
  group_by(cove) %>%
  mutate(
    cove_temp95_thresh = quantile(temperature_C, probs = 0.95, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(site_id, cove, side_code) %>%
  summarise(
    prop_time_above_temp95 = mean(temperature_C > cove_temp95_thresh, na.rm = TRUE),
    pct_time_above_temp95 = 100 * mean(temperature_C > cove_temp95_thresh, na.rm = TRUE),
    .groups = "drop"
  )

logger_extreme

# Join it onto other data
env_data_95 <- env_data %>%
  left_join(
    logger_extreme %>%
      select(site_id, prop_time_above_temp95, pct_time_above_temp95),
    by = "site_id"
  )

env_data_95










# APPENDIX A1.10b. Mean temperature versus extreme temperature
# H1 comparison using cove-specific 95th percentile threshold exceedance

# Prepare data
env_95 <- env_data_95 %>%
  mutate(
    cove = str_extract(site_id, "^C[0-9]+"),
    side = str_extract(site_id, "(SF|NF)$"),
    side = factor(side, levels = c("NF", "SF"))
  )

# Make paired format
env_pairs_95 <- env_95 %>%
  select(cove, side, air_temp_mean, pct_time_above_temp95) %>%
  pivot_wider(
    names_from = side,
    values_from = c(air_temp_mean, pct_time_above_temp95),
    names_sep = "__"
  ) %>%
  mutate(
    d_temp_mean = air_temp_mean__SF - air_temp_mean__NF,
    d_temp95    = pct_time_above_temp95__SF - pct_time_above_temp95__NF
  )

# Define the function
exact_signflip_p <- function(diffs) {
  diffs <- diffs[!is.na(diffs)]
  n <- length(diffs)
  
  all_signs <- expand.grid(rep(list(c(-1, 1)), n))
  perm_means <- apply(all_signs, 1, function(s) mean(as.numeric(s) * diffs))
  obs <- mean(diffs)
  
  mean(abs(perm_means) >= abs(obs))
}

# Exact sign-flip test to mean-temp and extreme-temp
p_temp_mean <- exact_signflip_p(env_pairs_95$d_temp_mean)
p_temp95    <- exact_signflip_p(env_pairs_95$d_temp95)

#define for bootstrapping
boot_mean_diff <- function(data, indices) {
  d <- data[indices]
  mean(d, na.rm = TRUE)
}

boot_ci_mean <- function(diffs, R = 10000) {
  diffs <- diffs[!is.na(diffs)]
  
  b <- boot(diffs, statistic = boot_mean_diff, R = R)
  ci <- boot.ci(b, type = "perc")
  
  tibble(
    mean_diff = mean(diffs),
    ci_low = ci$percent[4],
    ci_high = ci$percent[5]
  )
}

#Summarise mean temp
temp_mean_summary <- boot_ci_mean(env_pairs_95$d_temp_mean) %>%
  mutate(
    metric = "Mean air temperature (°C)",
    p_value = p_temp_mean
  )

#Summarise 95th percentile temp
temp95_summary <- boot_ci_mean(env_pairs_95$d_temp95) %>%
  mutate(
    metric = "% time above cove-specific 95th percentile threshold",
    p_value = p_temp95
  )

# COmpare
h1_temp_compare <- bind_rows(temp_mean_summary, temp95_summary) %>%
  select(metric, mean_diff, ci_low, ci_high, p_value)

h1_temp_compare














# APPENDIX A1.11. Environmental drivers of community structure (H4)
# Replace mean temperature with extreme temp

# Comm matrix for H4 sesitivity
comm_h4_95 <- community_mat %>%
  column_to_rownames("site_id")

# Extract metadata and remove some columns
meta_h4_95 <- comm_h4_95 %>%
  rownames_to_column("site_id") %>%
  select(site_id, cove, side)

comm_h4_95 <- comm_h4_95 %>%
  select(-cove, -side)

# 4th-root transform and calculate Bray-Curtis
comm_4rt_95 <- comm_h4_95^(1/4)
dist_h4_95 <- vegdist(comm_4rt_95, method = "bray")

#env dataset using 95th percentile
env_h4_95 <- env_data_95 %>%
  mutate(
    cove = str_extract(site_id, "^C[0-9]+"),
    side = str_extract(site_id, "(SF|NF)$"),
    side = factor(side, levels = c("NF", "SF")),
    cove = factor(cove),
    temp95_z = as.numeric(scale(pct_time_above_temp95)),
    lux_z    = as.numeric(scale(lux_mean)),
    wave_z   = as.numeric(scale(relative_wave_exposure))
  ) %>%
  filter(site_id %in% rownames(comm_4rt_95)) %>%
  arrange(match(site_id, rownames(comm_4rt_95)))

stopifnot(all(env_h4_95$site_id == rownames(comm_4rt_95)))

# Inspect
cor(
  env_h4_95 %>% select(temp95_z, lux_z, wave_z),
  use = "pairwise.complete.obs"
)

# dbRDA
dbrda_95_appendix <- capscale(
  comm_4rt_95 ~ temp95_z + lux_z + wave_z,
  data = env_h4_95,
  distance = "bray"
)

# Test dbRDA and marginal contribution of each predictor
anova(dbrda_95_appendix, permutations = 999)
anova(dbrda_95_appendix, by = "margin", permutations = 999)

# PERMANOVA - side only
adon_side_95_appendix <- adonis2(
  dist_h4_95 ~ side,
  data = env_h4_95,
  permutations = 999,
  strata = env_h4_95$cove
)

adon_side_95_appendix

# Full PERMANOVA including env predictors
adon_full_95_appendix <- adonis2(
  dist_h4_95 ~ side + temp95_z + lux_z + wave_z,
  data = env_h4_95,
  permutations = 999,
  strata = env_h4_95$cove,
  by = "margin"
)

adon_full_95_appendix

#join 95th percentile env predictors -> species richness dataset
rich_h4_95 <- site_summary %>%
  left_join(
    env_h4_95 %>%
      select(site_id, temp95_z, lux_z, wave_z),
    by = "site_id"
  ) %>%
  mutate(
    side = factor(side, levels = c("north-facing", "south-facing")),
    cove = factor(cove)
  )

# Species richness model
m_rich_env_temp95 <- glm(
  richness ~ cove + side + temp95_z + lux_z + wave_z,
  data = rich_h4_95,
  family = poisson(link = "log")
)

summary(m_rich_env_temp95)
drop1(m_rich_env_temp95, test = "Chisq")









# APPENDIX A.1.12a. Sensitivity analyses excluding the outlier cove
# sensitivity analyses repeated after excluding cove C2 (extreme site was C2_SF)


# make the sensitivity datasets
cove_exclude <- "C2"

community_mat_sens <- community_mat %>%
  filter(cove != cove_exclude)

site_summary_sens <- site_summary %>%
  filter(cove != cove_exclude)

env_data_sens <- env_data %>%
  filter(!str_detect(site_id, paste0("^", cove_exclude, "_")))

limp_site_sens <- limp_site %>%
  filter(cove != cove_exclude)







# A1.12b. H1 sensitivity: environmental differences

# Prepare dataset
env_sens <- env_data_sens %>%
  mutate(
    cove = str_extract(site_id, "^C[0-9]+"),
    side = str_extract(site_id, "(SF|NF)$"),
    side = factor(side, levels = c("NF", "SF"))
  )

# NF vs SF and calculate within-cove differences
env_pairs_sens <- env_sens %>%
  select(cove, side, air_temp_mean, lux_mean, relative_wave_exposure) %>%
  pivot_wider(
    names_from = side,
    values_from = c(air_temp_mean, lux_mean, relative_wave_exposure),
    names_sep = "__"
  ) %>%
  mutate(
    d_temp = air_temp_mean__SF - air_temp_mean__NF,
    d_lux  = lux_mean__SF - lux_mean__NF,
    d_wave = relative_wave_exposure__SF - relative_wave_exposure__NF
  )

# EXact paired sign-flip p-values 
p_temp_sens <- exact_signflip_p(env_pairs_sens$d_temp)
p_lux_sens  <- exact_signflip_p(env_pairs_sens$d_lux)
p_wave_sens <- exact_signflip_p(env_pairs_sens$d_wave)

#Summarise
temp_summary_sens <- boot_ci_mean(env_pairs_sens$d_temp) %>%
  mutate(metric = "Mean air temperature (°C)", p_value = p_temp_sens)

lux_summary_sens <- boot_ci_mean(env_pairs_sens$d_lux) %>%
  mutate(metric = "Mean light intensity (lux)", p_value = p_lux_sens)

wave_summary_sens <- boot_ci_mean(env_pairs_sens$d_wave) %>%
  mutate(metric = "Relative wave exposure", p_value = p_wave_sens)

h1_summary_sens <- bind_rows(
  temp_summary_sens,
  lux_summary_sens,
  wave_summary_sens
) %>%
  select(metric, mean_diff, ci_low, ci_high, p_value)

h1_summary_sens






# APPENDIX A1.12c. H2 sensitivity: univariate community responses

# Paried cover and richness data for senstivity 
h2_pairs_sens <- site_summary_sens %>%
  mutate(
    side_code = str_extract(site_id, "(NF|SF)$")
  ) %>%
  select(cove, side_code, richness, total_cover) %>%
  pivot_wider(
    id_cols = cove,
    names_from = side_code,
    values_from = c(richness, total_cover),
    names_sep = "__"
  ) %>%
  mutate(
    d_richness = richness__SF - richness__NF,
    d_cover    = total_cover__SF - total_cover__NF
  )

# Paired t-tests
t_richness_sens <- t.test(
  h2_pairs_sens$richness__SF,
  h2_pairs_sens$richness__NF,
  paired = TRUE
)

t_cover_sens <- t.test(
  h2_pairs_sens$total_cover__SF,
  h2_pairs_sens$total_cover__NF,
  paired = TRUE
)

# Summarising results
cover_summary_sens <- boot_ci_mean(h2_pairs_sens$d_cover) %>%
  mutate(
    metric = "Total macroalgal cover (summed % cover)",
    test_statistic = unname(t_cover_sens$statistic),
    df = unname(t_cover_sens$parameter),
    p_value = t_cover_sens$p.value
  )

richness_summary_sens <- boot_ci_mean(h2_pairs_sens$d_richness) %>%
  mutate(
    metric = "Species richness",
    test_statistic = unname(t_richness_sens$statistic),
    df = unname(t_richness_sens$parameter),
    p_value = t_richness_sens$p.value
  )

h2_summary_sens <- bind_rows(cover_summary_sens, richness_summary_sens) %>%
  select(metric, mean_diff, ci_low, ci_high, test_statistic, df, p_value)

h2_summary_sens













# APPENDIX A.1.12d. H3 sensitivity: multivariate community structure

# Prepare sensitivity community matrix
comm_sens <- community_mat_sens %>%
  column_to_rownames("site_id")

# extract site metadata
meta_h3_sens <- comm_sens %>%
  rownames_to_column("site_id") %>%
  select(site_id, cove, side)

# Remove metadata so only species remains
comm_sens <- comm_sens %>%
  select(-cove, -side)

# 4th-root transformation and Bray-Curtis
comm_4rt_sens <- comm_sens^(1/4)
dist_h3_sens <- vegdist(comm_4rt_sens, method = "bray")

#Paired PERMANOVA
adon_h3_sens <- adonis2(
  dist_h3_sens ~ side,
  data = meta_h3_sens,
  permutations = 999,
  strata = meta_h3_sens$cove
)

adon_h3_sens

#multivariate dispersion
bd_h3_sens <- betadisper(dist_h3_sens, group = meta_h3_sens$side)

anova(bd_h3_sens)
permutest(bd_h3_sens, permutations = 999)

# NMDS 
set.seed(123)
nmds_h3_sens <- metaMDS(
  comm_4rt_sens,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  wascores = FALSE
)

nmds_h3_sens$stress

# Extract site scores and join to metadata
scores_h3_sens <- as.data.frame(scores(nmds_h3_sens, display = "sites")) %>%
  rownames_to_column("site_id") %>%
  left_join(meta_h3_sens, by = "site_id")

scores_h3_sens

# Convert into presene-absence data
comm_pa_sens <- decostand(comm_sens, method = "pa")

# Pairwise beta diveristy -> turnover vs nestedness
pair_beta_by_cove_sens <- map_dfr(unique(meta_h3_sens$cove), function(cv) {
  ids <- meta_h3_sens %>%
    filter(cove == cv) %>%
    arrange(side) %>%
    pull(site_id)
  
  mat_cv <- comm_pa_sens[ids, , drop = FALSE]
  
  beta_pair <- betapart::beta.pair(mat_cv, index.family = "jaccard")
  
  tibble(
    cove = cv,
    turnover = as.matrix(beta_pair$beta.jtu)[1, 2],
    nestedness = as.matrix(beta_pair$beta.jne)[1, 2],
    total_beta = as.matrix(beta_pair$beta.jac)[1, 2]
  )
})

pair_beta_by_cove_sens

#Summarise
beta_summary_sens <- pair_beta_by_cove_sens %>%
  summarise(
    mean_turnover = mean(turnover),
    mean_nestedness = mean(nestedness),
    mean_total_beta = mean(total_beta),
    turnover_prop = mean(turnover / total_beta),
    nestedness_prop = mean(nestedness / total_beta)
  )

beta_summary_sens

# Identify taxa present on NF and SF
nf_species_sens <- meta_h3_sens %>%
  filter(side == "north-facing") %>%
  pull(site_id) %>%
  {\(ids) colnames(comm_pa_sens)[colSums(comm_pa_sens[ids, , drop = FALSE]) > 0]}()

sf_species_sens <- meta_h3_sens %>%
  filter(side == "south-facing") %>%
  pull(site_id) %>%
  {\(ids) colnames(comm_pa_sens)[colSums(comm_pa_sens[ids, , drop = FALSE]) > 0]}()

# Shared and aspect-unique species pools
species_pool_summary_sens <- tibble(
  shared_species = length(intersect(nf_species_sens, sf_species_sens)),
  nf_unique_species = length(setdiff(nf_species_sens, sf_species_sens)),
  sf_unique_species = length(setdiff(sf_species_sens, nf_species_sens))
)

species_pool_summary_sens










# APPENDIX A1.12e. H4 sensitivity: environmental drivers

# Prepare sensitivity matrix
comm_h4_sens <- community_mat_sens %>%
  column_to_rownames("site_id")

# Extract and remove metadata
meta_h4_sens <- comm_h4_sens %>%
  rownames_to_column("site_id") %>%
  select(site_id, cove, side)

comm_h4_sens <- comm_h4_sens %>%
  select(-cove, -side)

# 4th root and Bray-Curtis
comm_4rt_h4_sens <- comm_h4_sens^(1/4)
dist_h4_sens <- vegdist(comm_4rt_h4_sens, method = "bray")

# Prepare dataset
env_h4_sens <- env_data_sens %>%
  mutate(
    cove = str_extract(site_id, "^C[0-9]+"),
    side = str_extract(site_id, "(SF|NF)$"),
    side = factor(side, levels = c("NF", "SF")),
    cove = factor(cove),
    temp_z = as.numeric(scale(air_temp_mean)),
    lux_z  = as.numeric(scale(lux_mean)),
    wave_z = as.numeric(scale(relative_wave_exposure))
  ) %>%
  filter(site_id %in% rownames(comm_4rt_h4_sens)) %>%
  arrange(match(site_id, rownames(comm_4rt_h4_sens)))

# CHeck allignment
stopifnot(all(env_h4_sens$site_id == rownames(comm_4rt_h4_sens)))

# dbRDA
dbrda_mean_sens <- capscale(
  comm_4rt_h4_sens ~ temp_z + lux_z + wave_z,
  data = env_h4_sens,
  distance = "bray"
)

anova(dbrda_mean_sens, permutations = 999)
anova(dbrda_mean_sens, by = "margin", permutations = 999)

# PERMANOVA - shore aspect only
adon_side_sens <- adonis2(
  dist_h4_sens ~ side,
  data = env_h4_sens,
  permutations = 999,
  strata = env_h4_sens$cove
)

adon_side_sens

#Full PERMANOVA
adon_full_sens <- adonis2(
  dist_h4_sens ~ side + temp_z + lux_z + wave_z,
  data = env_h4_sens,
  permutations = 999,
  strata = env_h4_sens$cove,
  by = "margin"
)

adon_full_sens

# standardised env predictors -> species-richness dataset
rich_h4_sens <- site_summary_sens %>%
  left_join(
    env_h4_sens %>% select(site_id, temp_z, lux_z, wave_z),
    by = "site_id"
  ) %>%
  mutate(
    side = factor(side, levels = c("north-facing", "south-facing")),
    cove = factor(cove)
  )

# sensitivity
m_rich_env_sens <- glm(
  richness ~ cove + side + temp_z + lux_z + wave_z,
  data = rich_h4_sens,
  family = poisson(link = "log")
)

summary(m_rich_env_sens)

#Test indepednent contribution of each predictor
drop1(m_rich_env_sens, test = "Chisq")












# APPENDIX A.1.12f. H5 sensitivity: grazer effects


# Prepare paried limpet-abudnance data
limp_pairs_sens <- limp_site_sens %>%
  mutate(
    side_code = str_extract(site_id, "(NF|SF)$")
  ) %>%
  select(cove, side_code, mean_limp) %>%
  pivot_wider(
    id_cols = cove,
    names_from = side_code,
    values_from = mean_limp,
    names_glue = "mean_limp__{side_code}"
  ) %>%
  mutate(
    d_limp = mean_limp__SF - mean_limp__NF
  )

# Paried t-test 
t_limp_sens <- t.test(
  limp_pairs_sens$mean_limp__SF,
  limp_pairs_sens$mean_limp__NF,
  paired = TRUE
)

#Summarise sensitivity
limp_summary_sens <- boot_ci_mean(limp_pairs_sens$d_limp) %>%
  mutate(
    metric = "Mean limpet abundance per quadrat",
    test_statistic = unname(t_limp_sens$statistic),
    df = unname(t_limp_sens$parameter),
    p_value = t_limp_sens$p.value
  ) %>%
  select(metric, mean_diff, ci_low, ci_high, test_statistic, df, p_value)

limp_summary_sens

# Join limpets to macroalgal richness
limp_macro_sens <- limp_site_sens %>%
  left_join(
    site_summary_sens %>% select(site_id, richness, total_cover, cove),
    by = c("site_id", "cove")
  ) %>%
  mutate(
    cove = factor(cove)
  )

#Sensitivity richness with limpet abundance
m_limp_rich_sens <- glm(
  richness ~ mean_limp + cove,
  data = limp_macro_sens,
  family = quasipoisson(link = "log")
)

summary(m_limp_rich_sens)
drop1(m_limp_rich_sens, test = "F")

#Sensitivity cover with limpet abudnance 
m_limp_cover_sens <- lm(
  total_cover ~ mean_limp + cove,
  data = limp_macro_sens
)

summary(m_limp_cover_sens)
drop1(m_limp_cover_sens, test = "F")












# APPENDIX A.1.12f H6 sensitivity: functional diversity and CWM PCA


# Site-level indices 
fd_idx_sens <- fd_idx %>%
  filter(cove != cove_exclude)

# Prepare FRic sensitivity
fric_diff_sens <- fd_idx_sens %>%
  select(cove, side_code, FRic) %>%
  pivot_wider(names_from = side_code, values_from = FRic) %>%
  filter(!is.na(NF) & !is.na(SF)) %>%
  mutate(d_FRic = SF - NF)

# t-test FRic
t_fric_sens <- t.test(fric_diff_sens$d_FRic)

# Summarise
fric_summary_sens <- tibble(
  metric = "Functional richness (FRic)",
  n_pairs = length(fric_diff_sens$d_FRic),
  mean_diff = mean(fric_diff_sens$d_FRic),
  ci_low = t_fric_sens$conf.int[1],
  ci_high = t_fric_sens$conf.int[2],
  test_statistic = unname(t_fric_sens$statistic),
  df = unname(t_fric_sens$parameter),
  p_value = t_fric_sens$p.value
)

fric_summary_sens

# FDis sensitivity
# Note: C2 was already excluded for FDis because C2_SF had <2 trait-bearing taxa
# So sensitivity will  make no difference 
fdis_diff_sens <- fd_idx_sens %>%
  select(cove, side_code, FDis, n_trait_taxa) %>%
  pivot_wider(
    names_from = side_code,
    values_from = c(FDis, n_trait_taxa)
  ) %>%
  filter(
    !is.na(FDis_NF) & !is.na(FDis_SF),
    n_trait_taxa_NF >= 2,
    n_trait_taxa_SF >= 2
  ) %>%
  mutate(d_FDis = FDis_SF - FDis_NF)

# Paried t-test for FDis
t_fdis_sens <- t.test(fdis_diff_sens$d_FDis)

# Summarise
fdis_summary_sens <- tibble(
  metric = "Functional dispersion (FDis)",
  n_pairs = length(fdis_diff_sens$d_FDis),
  mean_diff = mean(fdis_diff_sens$d_FDis),
  ci_low = t_fdis_sens$conf.int[1],
  ci_high = t_fdis_sens$conf.int[2],
  test_statistic = unname(t_fdis_sens$statistic),
  df = unname(t_fdis_sens$parameter),
  p_value = t_fdis_sens$p.value
)

fdis_summary_sens

h6_summary_sens <- bind_rows(fric_summary_sens, fdis_summary_sens)
h6_summary_sens

# Recalculate CWM PCA without cove C2
cwm_df_sens <- cwm_df %>%
  filter(cove != cove_exclude)

cwm_num_sens <- cwm_df_sens %>%
  select(where(is.numeric)) %>%
  as.data.frame()

rownames(cwm_num_sens) <- cwm_df_sens$site_id

# Run sensitivity PCA
pca_cwm_sens <- prcomp(cwm_num_sens, scale. = TRUE)

#Extract site scores
cwm_scores_sens <- as.data.frame(pca_cwm_sens$x[, 1:2]) %>%
  rownames_to_column("site_id") %>%
  left_join(
    cwm_df_sens %>% select(site_id, cove, side),
    by = "site_id"
  )

# Calculate convex hulls
cwm_hulls_sens <- cwm_scores_sens %>%
  group_by(side) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

# View scores, hull coordinates and variance
cwm_scores_sens
cwm_hulls_sens
summary(pca_cwm_sens)

