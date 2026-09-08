####################################################################
## Single-experiment alpha-lattice mixed model
## Maize accessions | 2 reps | incomplete blocks | BLUPs of true
## genotypic performance (one environment/trial only)
####################################################################

## ---- 0. Packages -------------------------------------------------
# install.packages(c("lme4","lmerTest","tidyverse","sommer"))

library(tidyverse)
library(lme4)
library(lmerTest)
library(sommer)

## ---- 1. Expected data structure -----------------------------------
# Your raw data (SingleBLUPSdata) has columns GENOTYPE, BLOCK, REP
# (uppercase) — this script needs Genotype/Block/Rep internally, so
# rather than overwriting SingleBLUPSdata itself, it builds a SEPARATE
# renamed copy here. This keeps SingleBLUPSdata untouched and usable
# by the ANOVA script (which expects the original uppercase names),
# while this BLUP script works from its own copy without conflict.

stopifnot(all(c("GENOTYPE","BLOCK","REP") %in% colnames(SingleBLUPSdata)))
# If this stopifnot fails, SingleBLUPSdata currently has lowercase-first
# names already (Genotype/Block/Rep) rather than the original uppercase —
# check colnames(SingleBLUPSdata) and adjust the rename() below to match
# whatever it currently contains.

SingleBLUPSdata_renamed <- SingleBLUPSdata %>%
  rename(Genotype = GENOTYPE, Block = BLOCK, Rep = REP) %>%
  mutate(
    Rep      = factor(Rep),
    Block    = factor(paste(Rep, Block, sep = "_")),
    Genotype = factor(Genotype)
  )

df <- SingleBLUPSdata_renamed
# Everything below still refers to `df`, unchanged — only this loading
# step changed, so the rest of the script needs no further edits.

trait_names <- c("GrainYield","PlantHeight","EarHeight","DaysToAnthesis",
                  "DaysToSilking","EarLength","EarDiameter","KernelRows",
                  "KernelsPerRow","HundredKernelWt","ShellingPct",
                  "GrainMoisture","LodgingPct","ChlorophyllIndex")
# Replace with your real trait column names / count — the code below
# doesn't care how many traits are in the vector.

####################################################################
## 2. Fit the alpha-lattice model per trait
####################################################################
## Standard single-trial alpha-lattice model for BLUPs:
##   y = mu + Rep + (1|Block:Rep) + (1|Genotype) + e
## Rep is typically fixed (only 2 levels, and reps represent a fixed
## blocking factor of the design rather than a random sample of reps).
## Block within Rep is random (recovers inter-block information).
## Genotype is random -> BLUPs (shrunken toward the overall mean).

fit_one_trait <- function(trait, data) {
  f <- as.formula(paste0(trait, " ~ Rep + (1|Block) + (1|Genotype)"))
  lmer(f, data = data, REML = TRUE,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}

lmer_fits <- map(trait_names, ~ fit_one_trait(.x, df)) %>%
  set_names(trait_names)

# Check convergence / singular fits before trusting any results
convergence_check <- tibble(
  Trait = trait_names,
  singular = map_lgl(lmer_fits, isSingular),
  messages = map_chr(lmer_fits, ~ paste(.x@optinfo$conv$lme4$messages, collapse = "; "))
)
print(convergence_check)
# A trait with singular = TRUE usually means Block (or Genotype) variance
# estimated at ~0 — check whether blocking helped that trait at all
# (compare against a simpler RCBD-style model without Block via AIC).

####################################################################
## 3. Extract genotype BLUPs
####################################################################

get_genotype_blups <- function(model, trait) {
  intercept <- fixef(model)[["(Intercept)"]]
  # mean over Rep fixed-effect levels so BLUPs represent the trial-mean
  # genotype performance rather than the Rep-1 reference level
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

write_csv(genotype_blups_wide, "genotype_BLUPs_single_trial_wide.csv")
write_csv(genotype_blups_all,  "genotype_BLUPs_single_trial_long.csv")

####################################################################
## 4. Variance components + heritability
####################################################################
## Two heritability estimates are given:
##  (i)  "Standard" broad-sense heritability on an entry-mean basis:
##         H2 = Vg / (Vg + Ve/r)
##       Simple, widely reported, but assumes a balanced design.
##  (ii) Cullis et al. (2006) generalized heritability — more appropriate
##       for alpha-lattice / unbalanced designs where PEV varies by
##       genotype (some genotypes are better connected across blocks
##       than others):
##         H2_Cullis = 1 - (mean_vdBLUP / (2 * Vg))
##       where mean_vdBLUP is the mean pairwise prediction error variance
##       of the genotype BLUP differences, approximated here as
##       2 * mean(PEV_i) under the standard independence approximation.

n_rep <- nlevels(df$Rep)

get_varcomp_h2 <- function(model, trait, blups) {
  vc <- as.data.frame(VarCorr(model))
  Vg <- vc$vcov[vc$grp == "Genotype"]
  Vb <- vc$vcov[vc$grp == "Block"]
  Ve <- attr(VarCorr(model), "sc")^2

  H2_standard <- Vg / (Vg + Ve / n_rep)

  # Cullis heritability using PEV from ranef condVar (SE^2 = PEV)
  pev <- blups %>% filter(Trait == trait) %>% pull(SE) %>% {.^2}
  mean_vdBLUP <- 2 * mean(pev)   # standard approximation
  H2_cullis <- 1 - (mean_vdBLUP / (2 * Vg))

  tibble(Trait = trait, Vg = Vg, Vblock = Vb, Ve = Ve,
         H2_standard = H2_standard, H2_Cullis = H2_cullis)
}

varcomp_summary <- map2_dfr(lmer_fits, trait_names,
                             ~ get_varcomp_h2(.x, .y, genotype_blups_all))
write_csv(varcomp_summary, "variance_components_H2_single_trial.csv")
print(varcomp_summary)

####################################################################
## 5. Ranked genotype table (single trial)
####################################################################
## EDIT to match your actual breeding objective per trait (see note
## in the multi-environment ranking script — same principle applies).

trait_direction <- tibble::tribble(
  ~Trait,              ~higher_is_better,
  "GrainYield",         TRUE,
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

ranked_blups <- genotype_blups_all %>%
  left_join(trait_direction, by = "Trait") %>%
  mutate(higher_is_better = replace_na(higher_is_better, TRUE)) %>%
  group_by(Trait) %>%
  mutate(Rank = if_else(higher_is_better,
                         rank(-BLUP, ties.method = "min"),
                         rank(BLUP,  ties.method = "min"))) %>%
  ungroup() %>%
  select(-higher_is_better)

write_csv(ranked_blups, "genotype_BLUPs_single_trial_ranked_long.csv")

ranked_wide <- ranked_blups %>%
  select(Genotype, Trait, BLUP, Rank) %>%
  pivot_wider(names_from = Trait, values_from = c(BLUP, Rank),
              names_glue = "{Trait}_{.value}")
write_csv(ranked_wide, "genotype_BLUPs_single_trial_ranked_wide.csv")

####################################################################
## 6. Alternative engine: sommer (often more stable for small,
##    unbalanced alpha-lattice trials; gives heritability directly)
####################################################################

fit_sommer_one_trait <- function(trait, data) {
  data$y <- data[[trait]]
  mmer(
    y ~ Rep,
    random = ~ Block + Genotype,
    rcov   = ~ units,
    data   = data,
    verbose = FALSE
  )
}

sommer_fits <- map(trait_names, ~ fit_sommer_one_trait(.x, df)) %>%
  set_names(trait_names)

get_sommer_blups_h2 <- function(model, trait, n_rep) {
  u <- model$U$Genotype$y
  vc <- summary(model)$varcomp
  Vg <- vc["Genotype", "VarComp"]
  Ve <- vc["units", "VarComp"]
  H2 <- Vg / (Vg + Ve / n_rep)
  list(
    blups = tibble(Genotype = names(u), Trait = trait, BLUP_deviation = as.numeric(u)),
    h2    = tibble(Trait = trait, Vg = Vg, Ve = Ve, H2_standard = H2)
  )
}

sommer_results <- map2(sommer_fits, trait_names, ~ get_sommer_blups_h2(.x, .y, n_rep))
sommer_blups   <- map_dfr(sommer_results, "blups")
sommer_h2      <- map_dfr(sommer_results, "h2")

write_csv(sommer_blups, "sommer_genotype_BLUPs_single_trial.csv")
write_csv(sommer_h2,    "sommer_H2_single_trial.csv")

####################################################################
## Notes
####################################################################
# 1. Rep as fixed vs random: with only 2 reps, Rep is conventionally
#    treated as fixed in single-trial alpha-lattice analysis (as coded
#    here) since 2 levels give an unreliable variance-component estimate
#    if treated as random. This does not affect the Genotype BLUPs.
#
# 2. If a trait shows isSingular(model) == TRUE, the Block variance is
#    likely ~0 for that trait — inspect whether a simpler model without
#    Block fits equally well (anova/AIC). BLUPs are still usable either way.
#
# 3. Cullis heritability (H2_Cullis) is generally preferred over the
#    standard formula for alpha-lattice designs because genotypes are
#    not equally well-connected across incomplete blocks, so their
#    individual prediction error variances differ — the standard formula
#    implicitly assumes they don't.
