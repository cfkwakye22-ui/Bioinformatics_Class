####################################################################
## Single-experiment alpha-lattice mixed model
## Maize accessions | 2 REPs | incomplete BLOCKs | BLUPs of true
## genotypic performance (one environment/trial only)
####################################################################

## ---- 0. Packages -------------------------------------------------
# install.packages(c("lme4","lmerTest","tidyverse","sommer","pbkrtest","emmeans", "dplyr"))

install.packages("pbkrtest")
install.packages("sommer")

#Load installed packages

library(tidyverse) #
library(lme4)      # For fitting linear mixed-effects models
library(lmerTest)  # For p-values and ANOVA tables
library(sommer)    #
library(dplyr)     # For data manipulation
library(clipr)     # For moving data copyied in excel to Import Data into R
library(emmeans)   # print means
library(pbkrtest)  # 
library(dplyr)     # Called to hep use the 'select' function sitting in dplyr only

# ==============================================================================
#Import Data
#===============================================================================
## ---- 1. Expected data structure -----------------------------------
# One row per plot, single trial only:
#   REP      - factor, 2 levels
#   BLOCK    - incomplete BLOCK factor, nested in REP
#   GENOTYPE - factor, one level per accession
# Trait1..TraitN - numeric trait columns
# If BLOCK numbers REPeat across REPs (e.g., both REPs have BLOCKs 1-9),
# make the BLOCK ID unique within REP first:
# SingleBLUPSdata$BLOCK <- factor(paste(SingleBLUPSdata$REP, SingleBLUPSdata$BLOCK, sep = "_"))
#------------------------------------------------------------------------

##=================================================================
##Import Data by copying from clipboard
#-------------------------------------
SingleBLUPSdata = clipr::read_clip_tbl()
# REPlace 'SingleBLUPSdata' with your dataset name if different
# Ensure factors are correctly set as categorical variables:
str(SingleBLUPSdata)
SingleBLUPSdata$REP = as.factor(SingleBLUPSdata$REP)
SingleBLUPSdata$BLOCK = as.factor(SingleBLUPSdata$BLOCK)
SingleBLUPSdata$GENOTYPE = as.factor(SingleBLUPSdata$GENOTYPE)

# i. Rename factors first (if not already done) as it's in the Dataset
SingleBLUPSdata <- SingleBLUPSdata %>%
  rename(Genotype = GENOTYPE, Block = BLOCK, Rep = REP) %>%
  mutate(
    Rep      = factor(Rep),
    Block    = factor(paste(Rep, Block, sep = "_")),
    Genotype = factor(Genotype)
  )
  
trait_names <- c("POLLEN", "SILK", "ASI", "PLHT", "EHT", "STYG", "PASP", "HUSKC",
                  "RLPERC", "SLPERC", "EPP", "EASP", "EROT", "YIELD")
#trait_names <- stores that list of your real trait-column names/counts in an object called trait_names.
# — this code doesn't care how many traits are in the vector 
# and it's just a list the rest of the script loops over later.

##========================================================================
# OR import Data from csv file
#-------------------------------------
SingleBLUPSdata <- read_csv("filename for the_single_trial_data.csv") %>%
  mutate(
    REP      = factor(REP),
    BLOCK    = factor(paste(REP, BLOCK, sep = "_")),
    GENOTYPE = factor(GENOTYPE)
  )
trait_names <- c("POLLEN", "SILK", "ASI", "PLHT", "EHT", "STYG", "PASP", "HUSKC",
                 "RLPERC", "SLPERC", "EPP", "EASP", "EROT", "YIELD")
#trait_names <- stores that list of your real trait-column names/counts in an object called trait_names.
# — this code doesn't care how many traits are in the vector 
# and it's just a list the rest of the script loops over later (map(trait_names, ...)) 
# so the same model gets fit once per trait automatically, instead of you writing out 14 separate BLOCKs of code.


###################################################################
## 2. Fit the alpha-lattice model per trait
####################################################################
## Standard single-trial alpha-lattice model for BLUPs:
##   y = mu + REP + (1|BLOCK:REP) + (1|GENOTYPE) + e
## REP is typically fixed (for this analysis it's only 2 levels, and REPs REPresent as fixed BLOCKing factor of the design rather than a random sample of REPs).
## BLOCK within REP is random (recovers inter-BLOCK information).
## GENOTYPE is random -> BLUPs (shrunken toward the overall mean).

# ii. Re-run the function definition (unchanged, but make sure it's the current version)
fit_one_trait <- function(trait, data) {
  f <- as.formula(paste0(trait, " ~ Rep + (1|Block) + (1|Genotype)"))
  lmer(f, data = data, REML = TRUE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

# iii. REFIT all 14 models using the renamed data
lmer_fits <- map(trait_names, ~ fit_one_trait(.x, SingleBLUPSdata)) %>%
  set_names(trait_names)

# Summary of variance components
summary(fit_one_trait)
summary(lmer_fits)

# Check convergence / singular fits before trusting any results
convergence_check <- tibble(
  Trait = trait_names,
  singular = map_lgl(lmer_fits, isSingular),
  messages = map_chr(lmer_fits, ~ paste(.x@optinfo$conv$lme4$messages, collapse = "; "))
)
print(convergence_check)
# A trait with singular = TRUE usually means BLOCK (or GENOTYPE) variance
# were estimated at ~0 — check whether blocking helped that trait at all
# (compare against a simpler RCBD-style model without BLOCK via AIC).

####################################################################
## 3. Extract GENOTYPE BLUPs
####################################################################

get_genotype_blups <- function(model, trait) {
  intercept <- fixef(model)[["(Intercept)"]]
  rep_fx <- fixef(model)[grepl("^Rep", names(fixef(model)))]
  rep_mean_adj <- mean(c(0, rep_fx))
  
  re <- ranef(model, condVar = TRUE)$Genotype
  pv <- attr(re, "postVar")
  se <- sqrt(pv[1, 1, ])
  
  tibble(
    Genotype = rownames(re),
    Trait = trait,
    BLUP_deviation = re[, 1],
    BLUP = intercept + rep_mean_adj + re[, 1],
    SE = se
  )
}

genotype_blups_all <- map2_dfr(lmer_fits, trait_names, get_genotype_blups)

genotype_blups_wide <- genotype_blups_all %>%
dplyr::select(Genotype, Trait, BLUP) %>%
pivot_wider(names_from = Trait, values_from = BLUP)

### Quick check to confirm it looks right before moving on:
# Call column 'Genotype' in 'genotype_blups_all'.
colnames(genotype_blups_all)
dim(genotype_blups_wide)   # should be ~72 rows (genotypes) x 15 columns (Genotype + 14 traits)
head(genotype_blups_wide)  # peek at the first few rows

#Save to CSV:
write_csv(genotype_blups_wide, "genotype_BLUPs_single_trial_wide.csv")
write_csv(genotype_blups_all,  "genotype_BLUPs_single_trial_long.csv")

# Locate file on laptop
getwd()

####################################################################
## 4. Variance components + heritability + Repeatability
####################################################################
## Two heritability estimates are given:
##  (i)  "Standard" broad-sense heritability on an entry-mean basis:
##         H2 = Vg / (Vg + Ve/r)
##       Simple, widely reorted, but assumes a balanced design.

#==Script to use for the standard-only version
n_rep <- nlevels(SingleBLUPSdata$Rep)

get_h2_standard_only <- function(model, trait) {
  vc <- as.data.frame(VarCorr(model))
  Vg <- vc$vcov[vc$grp == "Genotype"]
  Vb <- vc$vcov[vc$grp == "Block"]
  Ve <- attr(VarCorr(model), "sc")^2
  H2_standard <- Vg / (Vg + Ve / n_rep)
  
  tibble(Trait = trait, Vg = Vg, Vblock = Vb, Ve = Ve, H2_standard = H2_standard)
}

##  (ii) Cullis et al. (2006) generalized heritability — more appropriate
##       for alpha-lattice / unbalanced designs where PEV varies by
##       GENOTYPE (some Genotypes are better connected across Blocks
##       than others):
##         H2_Cullis = 1 - (mean_vdBLUP / (2 * Vg))
##       where mean_vdBLUP is the mean pairwise prediction error variance
##       of the GENOTYPE BLUP differences, approximated here as
##       2 * mean(PEV_i) under the standard independence approximation.

#==Script to use for the Cullis-only version:
n_rep <- nlevels(SingleBLUPSdata$Rep)

get_h2_cullis_only <- function(model, trait, blups) {
  vc <- as.data.frame(VarCorr(model))
  Vg <- vc$vcov[vc$grp == "Genotype"]
  Vb <- vc$vcov[vc$grp == "Block"]
  Ve <- attr(VarCorr(model), "sc")^2
  
  pev <- blups %>% dplyr::filter(Trait == trait) %>% dplyr::pull(SE) %>% {.^2}
  mean_vdBLUP <- 2 * mean(pev)
  H2_cullis <- 1 - (mean_vdBLUP / (2 * Vg))
  
  tibble(Trait = trait, Vg = Vg, Vblock = Vb, Ve = Ve, H2_Cullis = H2_cullis)
}

#### Running the scripts for both H2 methods in one pass gives 
#### a table with both H2_standard and H2_Cullis side by side, per trait.

n_rep <- nlevels(SingleBLUPSdata$Rep)

get_varcomp_h2 <- function(model, trait, blups) {
  vc <- as.data.frame(VarCorr(model))
  Vg <- vc$vcov[vc$grp == "Genotype"]
  Vb <- vc$vcov[vc$grp == "Block"]
  Ve <- attr(VarCorr(model), "sc")^2
  H2_standard <- Vg / (Vg + Ve / n_rep)
  
  pev <- blups %>% dplyr::filter(Trait == trait) %>% dplyr::pull(SE) %>% {.^2}
  mean_vdBLUP <- 2 * mean(pev)
  H2_cullis <- 1 - (mean_vdBLUP / (2 * Vg))
  
  tibble(Trait = trait, Vg = Vg, Vblock = Vb, Ve = Ve,
         H2_standard = H2_standard, H2_Cullis = H2_cullis)
}

# Run the function for H2 methods above as-is, then just use/report whichever 
# column you prefer from the output:This gives you a table 
# with both H2_standard and H2_Cullis side by side, per trait.
varcomp_summary <- map2_dfr(lmer_fits, trait_names,
                            ~ get_varcomp_h2(.x, .y, genotype_blups_all))
write_csv(varcomp_summary, "variance_components_H2_single_trial.csv")
print(varcomp_summary)


####################################################################
############ Repeatability estimates per trait (single-trial alpha lattice)
## Run AFTER the H2 section (needs: lmer_fits, trait_names, n_rep,
## varcomp_summary already built)
####################################################################
## Repeatability vs heritability — the distinction:
##  - Heritability (H2) asks: what proportion of phenotypic variance
##    among genotypes is genetic?
##  - Repeatability (R) asks: how consistent is the SAME genotype's
##    measured performance across its repeated observations (here,
##    across the 2 reps)? In a single-trial context like this, R is
##    calculated the same way but is more naturally interpreted as
##    "how reliably would I get the same ranking if I repeated this
##    trial layout again" — useful for traits measured with any
##    subjectivity/noise (e.g., visual scores) or before deciding
##    how many reps you'd need in a future trial.
##
## Two versions are given, matching common usage:
##  (a) Plot-basis (single-observation) repeatability — treats a
##      single plot record as "one observation" of a genotype:
##        R_plot = Vg / (Vg + Vblock + Ve)
##  (b) Entry-mean-basis repeatability — treats a genotype's overall
##      trial mean (averaged over its reps and the blocks it appeared
##      in) as "one observation":
##        R_mean = Vg / (Vg + Vblock/n_block + Ve/n_rep)
##      This is directly comparable to H2_standard, but additionally
##      accounts for block variance rather than ignoring it.
#   Any trait with Vg = 0, expect R = 0; That's a correct result, not something to debug further

## NB: ---- Make sure required objects exist --------------------------
if (!exists("varcomp_summary")) {
  stop("Run the H2 section first (get_varcomp_h2 / varcomp_summary) before this script.")
}

### Before running, quick check as usual:
n_rep   <- nlevels(SingleBLUPSdata$Rep)
n_block <- SingleBLUPSdata %>%
  distinct(Rep, Block) %>%
  nrow() %>%
  {. / n_rep}   # blocks per rep (assumes equal block count per rep)

####################################################################
## 1. Calculate repeatability from the variance components already
##    stored in varcomp_summary (Vg, Vblock, Ve)
####################################################################

repeatability_summary <- varcomp_summary %>%
  mutate(
    R_plot_basis = Vg / (Vg + Vblock + Ve),
    R_entry_mean = Vg / (Vg + Vblock / n_block + Ve / n_rep)
  ) %>%
  dplyr::select(Trait, Vg, Vblock, Ve, R_plot_basis, R_entry_mean)

print(repeatability_summary)
write_csv(repeatability_summary, "repeatability_single_trial.csv")

####################################################################
## 2. Optional: combine repeatability with the H2 table into one
##    single summary file, since they're calculated from the same
##    variance components and usually get reported together
####################################################################

full_summary <- varcomp_summary %>%
  left_join(
    repeatability_summary %>% select(Trait, R_plot_basis, R_entry_mean),
    by = "Trait"
  )

write_csv(full_summary, "variance_components_H2_repeatability_single_trial.csv")
print(full_summary)

####################################################################
## Notes
####################################################################
# 1. Vg = 0 traits (e.g. EPP from your earlier output) will also show
#    R = 0 here, for the same reason H2 came out 0/NaN — no detectable
#    genotype variance means no repeatability either. This is a
#    legitimate result to report, not a coding issue.
#
# 2. R_plot_basis will generally be LOWER than R_entry_mean, since a
#    single plot observation is noisier than an averaged genotype
#    mean. If you're deciding how many reps to use in a FUTURE trial
#    for a given trait, R_plot_basis is the more relevant number to
#    reason from (it tells you how noisy one observation is before
#    any averaging benefit is applied).
#
# 3. If a trait had a singular Block fit (Vblock = 0, as with POLLEN
#    earlier), R_entry_mean collapses toward the same value as
#    H2_standard, since the block term contributes nothing either way.


####################################################################
## Base Index (BI) for selection under drought stress
## Badu-Apraku et al. multiple-trait-based selection index:
##   BI = (2 x YIELD) + EPP - ASI - PASP - EASP - STYG
##
## Provided for BOTH:
##   (a) a single trial's genotype BLUPs
##   (b) combined across-environment (multi-trial) genotype BLUPs
##
## Run AFTER you have a genotype BLUP wide table already built —
## either `genotype_blups_wide` from the single-trial script, or
## the combined-BLUP wide table from the multi-environment script.
####################################################################


## ---- Base Index (BI)--reminder (Badu-Apraku index) ------------------
## YIELD - grain yield: HIGHER is better  -> positive weight, doubled
##         (yield is given double weight, reflecting its primary
##         importance as a selection target)
## EPP   - ears per plant: HIGHER is better -> positive weight
## ASI   - anthesis-silking interval: LOWER is better (shorter ASI =
##         less drought-induced flowering asynchrony) -> negative weight
## PASP  - plant aspect score: LOWER is better (visual score, lower
##         score = healthier/more desirable plant type) -> negative weight
## EASP  - ear aspect score: LOWER is better (same logic as PASP)
##         -> negative weight
## STYG  - stay-green / senescence score: LOWER is better (lower score
##         = more stress tolerance, e.g. greener at maturity)
##         -> negative weight
##
## IMPORTANT: double-check these directions and the exact scoring
## scale (e.g., 1-9 vs 1-5) used in YOUR dataset before trusting the
## sign convention below — visual scores are sometimes coded in the
## opposite direction depending on the scoring sheet used.

####################################################################
## 1. Reusable function: standardize traits + compute BI + rank
####################################################################
## Traits are standardized (converted to Z-scores: mean 0, SD 1)
## before combining, since YIELD (t/ha), EPP (count), and the
## aspect/ASI scores (small integer scales) are on completely
## different measurement scales — combining raw values directly
## would let whichever trait has the largest numeric range dominate
## the index regardless of its actual importance.

calculate_BI <- function(blup_data, genotype_col = "Genotype") {
  
  required_traits <- c("YIELD", "EPP", "ASI", "PASP", "EASP", "STYG")
  missing_cols <- setdiff(required_traits, colnames(blup_data))
  if (length(missing_cols) > 0) {
    stop("Missing required trait columns: ", paste(missing_cols, collapse = ", "),
         "\nCheck your column names match exactly (case-sensitive).")
  }
 
  z_data <- blup_data %>%
    mutate(across(all_of(required_traits),
                  ~ as.numeric(scale(.x)),
                  .names = "Z_{.col}"))
  
  z_data %>%
    mutate(
      BI = (2 * Z_YIELD) + Z_EPP - Z_ASI - Z_PASP - Z_EASP - Z_STYG,
      BI_Rank = rank(-BI, ties.method = "min")   # rank 1 = best (highest BI)
    ) %>%
    arrange(BI_Rank) %>%
    dplyr::select(all_of(genotype_col), all_of(required_traits),
                  starts_with("Z_"), BI, BI_Rank)
}

####################################################################
## 2a. SINGLE TRIAL — apply to genotype_blups_wide
####################################################################
## Assumes `genotype_blups_wide` already exists from your single-trial
## script, with columns: Genotype, YIELD, EPP, ASI, PASP, EASP, STYG,
## (plus any other traits you fitted — extra columns are fine, only
## the 6 required ones are used).

BI_single_trial <- calculate_BI(genotype_blups_wide)

print(BI_single_trial)
write_csv(BI_single_trial, "BaseIndex_BI_single_trial.csv")

####################################################################
## 2b. MULTIPLE TRIALS (combined across-environment BLUPs)
####################################################################
## Assumes you have a combined-BLUP wide table from the multi-
## environment script — e.g. `genotype_blups_wide` from that script
## (rename here to avoid clashing with the single-trial object above
## if both are in the same R session).

# Example: if your combined multi-environment wide table is called
# `combined_genotype_blups_wide`, rename it below to match:
# combined_genotype_blups_wide <- genotype_blups_wide   # (from the MET script)

BI_multi_trial <- calculate_BI(combined_genotype_blups_wide)

print(BI_multi_trial)
write_csv(BI_multi_trial, "BaseIndex_BI_combined_multitrial.csv")

####################################################################
## 3. Optional: compare single-trial vs combined BI rankings
####################################################################
## Useful to see whether top genotypes selected from one trial alone
## agree with their ranking once multi-environment data is pooled —
## large disagreements flag genotypes whose apparent performance in
## a single trial may not hold up as true, stable drought tolerance.

BI_comparison <- BI_single_trial %>%
  dplyr::select(Genotype, BI_SingleTrial = BI, Rank_SingleTrial = BI_Rank) %>%
  full_join(
    BI_multi_trial %>%
      dplyr::select(Genotype, BI_Combined = BI, Rank_Combined = BI_Rank),
    by = "Genotype"
  ) %>%
  arrange(Rank_Combined)

write_csv(BI_comparison, "BaseIndex_BI_comparison_single_vs_combined.csv")
print(BI_comparison)

####################################################################
## Notes
####################################################################
# 1. Z-scores are calculated WITHIN each dataset separately (single
#    trial standardized against itself; combined BLUPs standardized
#    against themselves) — this is standard practice, since the point
#    is relative ranking of genotypes within a given evaluation, not
#    comparing raw index values across different trials/datasets.
#
# 2. If any of the 6 required trait columns is missing or misspelled
#    in your BLUP table, calculate_BI() will stop with a clear message
#    listing exactly which column(s) it couldn't find — check spelling/
#    capitalization against your actual trait_names first (same class
#    of issue we've hit repeatedly with Genotype/GENOTYPE case).
#
# 3. If ASI, PASP, EASP, or STYG in your dataset are scored on a scale
#    where a HIGHER number is actually better (opposite of what's
#    assumed here), flip the sign for that term in the BI formula
#    inside calculate_BI() accordingly — the formula's signs only make
#    sense given the "lower is better" scoring convention noted above.


####################################################################
## 5. Ranked GENOTYPE table (single trial)
####################################################################
## EDIT to match your actual breeding objective per trait (see note
## in the multi-environment ranking script — same principle applies).

trait_direction <- tibble::tribble(
  ~Trait,              ~higher_is_better,
  "GrainYIELD",         TRUE,
  "PlantHeight",        TRUE,
  "EarHeight",          TRUE,
  "DaysToAnthesis",     FALSE,
  "DaysToSilking",      FALSE,
  "EarLength",          TRUE,
  "EarDiameter",        TRUE,
  "KernelRows",         TRUE,
  "KernelsPerRow",      TRUE,
  "HundredKernelWt",    TRUE,
  "ShellingPct",        TRUE,
  "GrainMoisture",      FALSE,
  "LodgingPct",         FALSE,
  "ChlorophyllIndex",   TRUE
)

ranked_blups <- GENOTYPE_blups_all %>%
  left_join(trait_direction, by = "Trait") %>%
  mutate(higher_is_better = REPlace_na(higher_is_better, TRUE)) %>%
  group_by(Trait) %>%
  mutate(Rank = if_else(higher_is_better,
                         rank(-BLUP, ties.method = "min"),
                         rank(BLUP,  ties.method = "min"))) %>%
  ungroup() %>%
  select(-higher_is_better)

write_csv(ranked_blups, "GENOTYPE_BLUPs_single_trial_ranked_long.csv")

ranked_wide <- ranked_blups %>%
  select(GENOTYPE, Trait, BLUP, Rank) %>%
  pivot_wider(names_from = Trait, values_from = c(BLUP, Rank),
              names_glue = "{Trait}_{.value}")
write_csv(ranked_wide, "GENOTYPE_BLUPs_single_trial_ranked_wide.csv")

####################################################################
## 6. Alternative engine: sommer (often more stable for small,
##    unbalanced alpha-lattice trials; gives heritability directly)
####################################################################

fit_sommer_one_trait <- function(trait, data) {
  data$y <- data[[trait]]
  mmer(
    y ~ REP,
    random = ~ BLOCK + GENOTYPE,
    rcov   = ~ units,
    data   = data,
    verbose = FALSE
  )
}

sommer_fits <- map(trait_names, ~ fit_sommer_one_trait(.x, SingleBLUPSdata)) %>%
  set_names(trait_names)

get_sommer_blups_h2 <- function(model, trait, n_REP) {
  u <- model$U$Genotype$y
  vc <- summary(model)$varcomp
  Vg <- vc["Genotype", "VarComp"]
  Ve <- vc["units", "VarComp"]
  H2 <- Vg / (Vg + Ve / n_REP)
  list(
    blups = tibble(Genotype = names(u), Trait = trait, BLUP_deviation = as.numeric(u)),
    h2    = tibble(Trait = trait, Vg = Vg, Ve = Ve, H2_standard = H2)
  )
}

sommer_results <- map2(sommer_fits, trait_names, ~ get_sommer_blups_h2(.x, .y, n_REP))
sommer_blups   <- map_SingleBLUPSdatar(sommer_results, "blups")
sommer_h2      <- map_SingleBLUPSdatar(sommer_results, "h2")

write_csv(sommer_blups, "sommer_GENOTYPE_BLUPs_single_trial.csv")
write_csv(sommer_h2,    "sommer_H2_single_trial.csv")

####################################################################
## Notes
####################################################################
# 1. REP as fixed vs random: with only 2 REPs, REP is conventionally
#    treated as fixed in single-trial alpha-lattice analysis (as coded
#    here) since 2 levels give an unreliable variance-component estimate
#    if treated as random. This does not affect the GENOTYPE BLUPs.
#
# 2. If a trait shows isSingular(model) == TRUE, the BLOCK variance is
#    likely ~0 for that trait — inspect whether a simpler model without
#    BLOCK fits equally well (anova/AIC). BLUPs are still usable either way.
#
# 3. Cullis heritability (H2_Cullis) is generally preferred over the
#    standard formula for alpha-lattice designs because GENOTYPEs are
#    not equally well-connected across incomplete BLOCKs, so their
#    individual prediction error variances differ — the standard formula
#    implicitly assumes they don't.
