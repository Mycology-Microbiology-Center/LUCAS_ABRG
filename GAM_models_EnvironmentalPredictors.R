#!/usr/bin/env Rscript

## Environmental predictors importance

## Usage:
# ./GAM_models_EnvironmentalPredictors.R --geneid "1"

## Input data = list with data.tables


############################################## Parse input parameters

## Check time
start_time <- Sys.time()


cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse))

## Parse arguments
option_list <- list(
  make_option(c("-g", "--geneid"), action="store", default=1L, type="integer", help="Slot number in the input list")
  )
opt <- parse_args(OptionParser(option_list=option_list))


## Assign variables
GENE <- opt$geneid

## Log assigned variables
cat(paste("Slot ID: ", GENE, "\n", sep=""))

cat("\n")


############################################## Load packages and data

cat("Loading R packages...\n")

load_pckg <- function(pkg = "data.table"){
  suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
  cat(paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("data.table")
load_pckg("mgcv")
load_pckg("ggplot2")
load_pckg("marginaleffects")
load_pckg("effectsize")
load_pckg("plyr")

## Number of CPU threads
setDTthreads(threads = 10)  # for data.table

## Activate parallel computation for marginaleffects
load_pckg("future.apply")
plan(strategy = multicore, workers = 10)

cat("\n")

############################################## Main pipeline

## Load the data
cat("Loading the input data..\n")
DATT <- readRDS("ARG_and_AB_READCOUNTS_for_Regression.RData")


## Main effects
main_predz <- c("pH_H2O", "H2O_content_volumetric", "CN", "Organic_carbon",
  "Annual_Mean_Temperature", "Annual_Precipitation",
  "Phosphorus_content", "Carbonate_content", "Clay_content", "Electrical_conductivity",
  "MeanTemperature_Avg31day", "PrecipitationSum_Avg31day")

## Predictor interactions
interacts <- c("pH_H2O,CN", "pH_H2O,Annual_Mean_Temperature")

## Function to fit models
fit_model <- function(datt, main_predz, interacts, spatial.modtype = "smoother"){

  cat("..Preparing formulae\n")

  ## Prepare the formula for main effects (add interactions and the dimension of the basis)
  prep_maineffs <- function(main_predz, interact = NULL, k = 3, spatial.modtype = "none"){
    if(is.null(interact)){
      mm <- paste0("s(", main_predz, ", k = ", k, ")")
    } else {
      mm <- paste0("s(", main_predz, ", k = ", k, ", by = ", interact, ")")
    }
    mm <- paste(mm, collapse = " + ")

    ## Interaction of latitude/longitude in form of a smoother term
    if(spatial.modtype %in% "smoother"){
      mm <- paste0(mm, " + s(GPS_LAT, GPS_LONG)")
    }

    return(mm)
  }

  ## Tensor product interaction (ti) - terms should be present in main effects as well!
  interacts <- paste0("ti(", interacts, ")")
  interacts <- paste(interacts, collapse = " + ")
  #     Anisotropic smooth for interactions (t2) - already includes main effs

  ## Formulae with interactions
  FRMS <- list(
    LC_simplest = paste0("Resids_Abundance ~ LC_simplest + ", prep_maineffs(main_predz, interact = "LC_simplest", k = 3, spatial.modtype = spatial.modtype)),
    Bfn_simplest_orig = paste0("Resids_Abundance ~ Bfn_simplest_orig + ", prep_maineffs(main_predz, interact = "Bfn_simplest_orig", k = 3, spatial.modtype = spatial.modtype))
    )

  ## Formulae without interactions
  FRMSNI <- list(
    LC_simplest = paste0("Resids_Abundance ~ LC_simplest + ", prep_maineffs(main_predz, interact = NULL, k = 3, spatial.modtype = spatial.modtype)),
    Bfn_simplest_orig = paste0("Resids_Abundance ~ Bfn_simplest_orig + ", prep_maineffs(main_predz, interact = NULL, k = 3, spatial.modtype = spatial.modtype))
    )

  ## Model wrapper
  run_the_model <- function(frm, x, verbose = TRUE){

    cat("..Fitting a model\n")
    mod <- try( gam(
      formula = as.formula(frm),
      family = gaussian(),
      data = x,
      select = TRUE, gamma = 1  # add an extra penalty to each term so that it can be penalized to zero
      ), silent = ! verbose )

    return(mod)
  }

  ## Fit the models with ecoregions and interactions
  cat("Fitting models with interactions..\n")
  MODS_ECOREG <- llply(.data = FRMS, .fun = function(m){ run_the_model(frm = m, x = datt) })

  ## Fit the models without ecoregion interactions (for overall variable importance)
  cat("Fitting models without interactions..\n")
  MODS_ECOREGNI <- llply(.data = FRMSNI, .fun = function(m){ run_the_model(frm = m, x = datt) })

  ## Extract model summaries and effects
  extract_effects <- function(mod){

    if(!"try-error" %in% class(mod)){

    ## Variable importance
    ## how much variance in the response variables is accounted for by the explanatory variables?
    ## Epsilon squared (unbiased estimator of eta squared) is analogous to R2
    cat("..Estimating variable importance\n")
    imp <- effectsize::epsilon_squared(mod, ci = 0.95, alternative = "two.sided")

    ## Land-use term name
    lu <- names(attr(x = mod$pterms, which = "dataClasses"))[2]

    ## Estimate marginal effects ("AME" of marginaleffects)
    cat("..Estimating marginal effects\n")
    ME <- marginaleffects(model = mod, newdata = NULL, variables = main_predz,
      type = "response", conf_level = 0.95,
      eps = 1e-04  # “step” size to use when calculating numerical derivatives
      )

    ## Overtall average marginal effects
    MA <- summary(ME)

    ## Median marginal effects, by land type
    cat("..Marginal effects summary\n")
    MM <- marginaleffects(model = mod, newdata = NULL, variables = main_predz,
      type = "response", conf_level = 0.95,
      eps = 1e-04,
      by = lu)
    setDT(MM)
    MM[ , WithinVar := lu ]

    eff <- list()
    eff$imp <- imp    # variable importance (full model)
    eff$MM <- MM      # median slopes by ecoregion
    eff$MA <- MA      # overall average marginal effects
    eff$ME <- ME      # observation-level effects

    } else {
      eff <- NULL
    } 
    return(eff)
  } # end of `extract_effects`


  ## Effects for all models
  cat("Effects for all models..\n")
  EFFS <- llply(.data = MODS_ECOREG, .fun = extract_effects)

  ## Variable importance for non-interaction models
  cat("Variable importance for non-interaction models..\n")
  VI <- llply(.data = MODS_ECOREGNI, .fun = function(mod){
    imp <- effectsize::epsilon_squared(mod, ci = 0.95, alternative = "two.sided")
    return(imp)
    })

  res <- list()
  res$mod1 <- MODS_ECOREG     # model with interactions
  res$mod2 <- MODS_ECOREGNI   # model without interactions
  res$eff <- EFFS             # Variable importance + median slopes (models with interactions)
  res$vi <- VI                # Variable importance for non-interaction models

  ## Add dataset ID as metadata
  attr(res, "Dataset") <- datt$GroupName[1]

  return(res)
}


## Fit the models to the selected gene
cat("Running analysis\n")
RES <- fit_model(
  datt = DATT[[ GENE ]],
  main_predz = main_predz, interacts = interacts,
  spatial.modtype = "smoother")

## Export results
cat("Exporting data\n")
saveRDS(
  object = RES,
  file = paste0("Env_", GENE, ".RData"),
  compress = "xz")

cat("ALL DONE\n")

#####################

## Check time
end_time <- Sys.time()

tmm <- as.numeric(difftime(end_time, start_time, units = "min"))
cat("\nElapsed time: ", tmm, " minutes\n")

cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
