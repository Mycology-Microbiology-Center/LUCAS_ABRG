#!/usr/bin/env Rscript

## Land-use differences within ecoregion

## Usage:
# ./GAM_models_Landuse.R --geneid "1"

## Input data = list with data.tables, `geneid` parameter defines the element in this list.


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
load_pckg("emmeans")
load_pckg("ggeffects")
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


## List of the models for testing
model_list <- list(
  ifPasturesCLCbased = "Resids_Abundance ~ LC_simplest_1 * ifPasturesCLCbased + s(MeanTemperature_Avg31day, k = 3) + s(PrecipitationSum_Avg31day, k = 3)",
  Bfn_simplest_orig  = "Resids_Abundance ~ LC_simplest_1 * Bfn_simplest_orig  + s(MeanTemperature_Avg31day, k = 3) + s(PrecipitationSum_Avg31day, k = 3)"
  )


## Function to fit a single model and estimate contrasts of interest
fit_model <- function(x, frm, spatial.modtype = "smoother"){

  ## To account for spatial autocorrelation, include latitude and longitude in a GAM model
  if(spatial.modtype %in% "smoother"){
    ## Interaction of latitude/longitude in form of a smoother term
    frm <- paste0(frm, " + s(GPS_LAT, GPS_LONG)")
  }

  cat(".Fitting a model\n")

  ## Fit the model
  mod <- try( gam(
    formula = as.formula(frm),
    family = gaussian(),
    data = x,
    select = FALSE)
    )

  if(!"try-error" %in% class(mod)){
    
    ## Check model performance
    # performance::performance(mod)
    # parameters::model_parameters(mod)
    # effectsize::eta_squared(mod)
  
    ## Terms of interest
    land_cover <- attr(x = mod$pterms, which = "term.labels")[1]
    ecoregion  <- attr(x = mod$pterms, which = "term.labels")[2]
    
    vrr <- list(
      v1 = unique(x[, ..land_cover])[[1]],
      v2 = unique(x[, ..ecoregion])[[1]])
    names(vrr) <- c(land_cover, ecoregion)

    ## Marginal effects
    lvl_land_cover <- as.character( levels(droplevels(mod$model[, land_cover])) )
    lvl_ecoregion <- as.character( levels(droplevels(mod$model[, ecoregion])) )
    avgvars <- expand.grid(lvl_land_cover, lvl_ecoregion)
    attr(avgvars, "out.attrs") <- NULL
    colnames(avgvars) <- c(land_cover, ecoregion)
    avgvars$MeanTemperature_Avg31day  <- 0   # Z-scores
    avgvars$PrecipitationSum_Avg31day <- 0
    avgvars$GPS_LAT <- mean(x$GPS_LAT, na.rm = T)
    avgvars$GPS_LONG <- mean(x$GPS_LONG, na.rm = T)

    ## Marginal means for combinations of categories
    cat("..Marginal means\n")
    eff <- predictions(mod, newdata = avgvars)
    setDT(eff)
    eff[, c("rowid", "type", "Resids_Abundance", "MeanTemperature_Avg31day", "PrecipitationSum_Avg31day") := NULL ]
    if(spatial.modtype %in% "smoother"){ eff[, c("GPS_LAT", "GPS_LONG") := NULL ] }
    setcolorder(eff, c(land_cover, ecoregion))

    ### Comparisons with average contrasts
    ## Land covers within Ecoregions
    cat("..Comparisons within Ecoregions\n")
    cnt_LC <- comparisons(model = mod, variables = land_cover,
      type = "response", conf_level = 0.95, contrast_factor = "pairwise",
      by = ecoregion)

    setDT(cnt_LC)
    setcolorder(cnt_LC, c("type", "term", ecoregion))
    cnt_LC[, contrast := gsub(pattern = "mean\\(", replacement = "", x = contrast) ]
    cnt_LC[, contrast := gsub(pattern = "\\)", replacement = "", x = contrast) ]

    ## Ecoregions within Land covers
    cat("..Comparisons within Land covers\n")
    cnt_ER <- comparisons(model = mod, variables = ecoregion,
      type = "response", conf_level = 0.95, contrast_factor = "pairwise",
      by = land_cover)

    setDT(cnt_ER)
    setcolorder(cnt_ER, c("type", "term", land_cover))
    cnt_ER[, contrast := gsub(pattern = "mean\\(", replacement = "", x = contrast) ]
    cnt_ER[, contrast := gsub(pattern = "\\)", replacement = "", x = contrast) ]

    res <- list()
    res$mod <- mod         # Fitted model
    res$eff <- eff         # Expected read counts
    res$cnt_LC <- cnt_LC   # Land covers within Ecoregions
    res$cnt_ER <- cnt_ER   # Ecoregions within Land covers

    ## Add dataset ID as metadata
    attr(res, "Dataset") <- x$GroupName[1]
  
  } else {
    res <- NULL
  }

  return(res)
}

## Function to fit all models
fit_all_models <- function(x){
  llply(.data = model_list, .fun = function(z){
    fit_model(x = x, frm = z)
    })
}

## Fit the models to the selected gene
cat("Running analysis\n")
RES <- fit_all_models( x = DATT[[ GENE ]] )

## Export results
cat("Exporting data\n")
saveRDS(
  object = RES,
  file = paste0("Int_", GENE, ".RData"),
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
