# Copyright 2025 David Novak
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## iidx internal module: 99_Aux.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines constants and auxiliary functions used in the workflow.


## Resolve fixed paths ----

## Paths for saving results
fpath_res00 <- 'Results_00_Preprocessing'
fpath_res01 <- 'Results_01_BatchEffectCorrection'
fpath_res02 <- 'Results_02_FlowSOMFeatureExtraction'
fpath_res03 <- 'Results_03_OutlierAndNoiseDetection'
fpath_res04 <- 'Results_04_StatisticalModelling'
fpath_res05 <- 'Results_05_Profiling'
fpath_res06 <- 'Results_06_Dashboard'

## Paths to pre-processing module outputs
fpath_fcs_comptrans <- 
  file.path(fpath_res00, 'PreprocessedInputs')
fpath_fcs_unstained_comptrans <- 
  file.path(fpath_res00, 'PreprocessedUnstainedSamples')
fname_agg_unstained <- 
  file.path(fpath_res00, 'UnstainedSamplesAggregate.fcs')

## Paths to normalisation module outputs
fpath_fcs_norm          <- file.path(fpath_res01, 'NormalisedInputs')
fpath_norm_model        <- file.path(fpath_res01, 'NormModelFiles')
fname_norm_model        <- file.path(fpath_res01, 'NormModel.RDS')
fpath_agg_normtrain_batchwise <- file.path(
  fpath_res01, 'BatchwiseAggregates_NormTrain'
)
fpath_agg_batchwise <- file.path(
  fpath_res01, 'BatchwiseAggregates'
)
fpath_norm_tmp           <- file.path(fpath_res01, 'Tmp')
fname_viz_emd_prenorm    <- file.path(fpath_res01, 'EMD_PreNorm.pdf')
fname_viz_emd_postnorm   <- file.path(fpath_res01, 'EMD_PostNorm.pdf')
fname_viz_emd_compressed <- file.path(fpath_res01, 'EMD_Compressed.pdf')
fname_viz_splines        <- file.path(fpath_res01, 'NormSplines.pdf')
fname_viz_signal_distributions <- file.path(
  fpath_res01, 'SignalDistributions.pdf'
)

## Paths to feature extraction module outputs
fname_agg            <- file.path(fpath_res02, 'Aggregate.fcs')
fname_fsom           <- file.path(fpath_res02, 'FlowSOMModel.RDS')
fname_metaclustering <- file.path(fpath_res02, 'FlowSOMMetaclustering.RDS')
fname_fsommary       <- file.path(fpath_res02, 'FlowSOMmary.pdf')
fname_viz_norm_lineage     <- file.path(fpath_res02, 'NormCheck_MCPerc.pdf')
fname_viz_norm_state_mfi   <- file.path(fpath_res02, 'NormCheck_MFI.pdf')
fname_viz_norm_state_pheno <- file.path(fpath_res02, 'NormCheck_Pheno.pdf')
fpath_features                 <- file.path(fpath_res02, 'FeatureMatrices')
fname_features_prefix_prenorm  <- file.path(fpath_res02, 'Features_PreNorm_')
fname_features_prefix_postnorm <- file.path(fpath_res02, 'Features_PostNorm_')
fname_features_prenorm  <- file.path(fpath_res02, 'Features_PreNorm.RDS')
fname_features_postnorm <- file.path(fpath_res02, 'Features_PostNorm.RDS')
fname_noise <- file.path(fpath_res02, 'NoiseLevels.RDS')

## Paths to outlier and noise detection module outputs
fname_sample_outliers     <- file.path(fpath_res03, 'SampleOutliers.RDS')
fname_viz_sample_outliers <- file.path(fpath_res03, 'SampleOutliers.pdf')
fname_robustness          <- file.path(fpath_res03, 'SignalRobustness.RDS')

## Paths to statistical modelling module outputs
fname_fsom_da          <- file.path(fpath_res04, 'FlowSOM_DA_Results.RDS')
fname_fsom_ds_mfi      <- file.path(fpath_res04, 'FlowSOM_DS-MFI_Results.RDS')
fname_fsom_ds_pheno    <- file.path(fpath_res04, 'FlowSOM_DS-Pheno_Results.RDS')
fname_cellcnn_da       <- file.path(fpath_res04, 'CellCnn_DA_Results.RDS')
fname_cellcnn_ds_mfi   <- file.path(fpath_res04, 'CellCnn_DS-MFI_Results.RDS')
fname_cellcnn_ds_pheno <- file.path(fpath_res04, 'CellCnn_DS-Pheno_Results.RDS')

## Paths to profiling module outputs
fname_mc_profiles  <- file.path(fpath_res05, 'MetaclusterProfiles.RDS')
fname_agg_gm       <- file.path(fpath_res05, 'AggregateGatingMatrix.RDS')
fname_agg_labels   <- file.path(fpath_res05, 'AggregateLabels.RDS')
fname_mc_quartiles <- file.path(fpath_res05, 'MetaclusterQuartiles.RDS')

## Resolve availability of unstained samples ----

unstained_samples_available <- (
  !is.null(fpath_fcs_unstained) && length(fpath_fcs_unstained)==1
)

## Define fixed hyperparameters ----

batch_ref_flag <- 999
cohorts        <- unique(annotation$Cohort)

## Resolve reference vs target batches ----

null_ref  <- is.null(batch_ref) || length(batch_ref)==0
null_norm <- is.null(batch_norm) || length(batch_norm)==0
if (null_ref || null_norm) {
  
  batches <- sort(unique(annotation$Batch))
  if (!is.null(batch_remove)) {
    batches <- batches[!batches%in%batch_remove]
  }
  if (null_ref) {
    batch_ref <- batches[1]
  }
  if (null_norm) {
    batch_norm <- batches[!batches%in%batch_ref]
  }
}

## Resolve whether input FCS files are needed ----

resolve_inputs <- TRUE
if (exists('input_files_needed')) {
  
  if (!input_files_needed) { # input files unnecessary
    resolve_inputs <- FALSE
  }
}
if (resolve_inputs) { # input files necessary
  
  source('99_CheckInputs.R')
  inputs <- resolve_input_fcs(
    annotation,
    fpath_fcs,
    fpath_fcs_comptrans,
    fpath_fcs_norm,
    normalise,
    fpath_fcs_unstained,
    cohorts,
    batch_remove,
    verbose = TRUE
  )
  fnames_fcs_raw       <- inputs[['fnames_fcs_raw']]
  fnames_fcs           <- fnames_fcs_raw
  fnames_fcs_comptrans <- inputs[['fnames_fcs_comptrans']]
  fnames_fcs_norm      <- inputs[['fnames_fcs_norm']]
  
  fnames_fcs_unstained_raw       <- inputs[['fnames_fcs_unstained_raw']]
  fnames_fcs_unstained_comptrans <- inputs[['fnames_fcs_unstained_comptrans']]
  
  annotation              <- inputs[['annotation']]
  
  fnames_fcs_preprocessed <- function() {
    if (normalise) fnames_fcs_norm else fnames_fcs_comptrans
  }
  fnames_fcs_unstained_preprocessed <- function() {
    fnames_fcs_unstained_comptrans
  }
} else {
  fnames_fcs_preprocessed <- function() {
    if (normalise) {
      dir(fpath_fcs_norm, pattern = '[.]fcs$', full.names = TRUE)
    } else {
      dir(fpath_fcs_comptrans, pattern = '[.]fcs$', full.names = TRUE)
    }
  }
  fnames_fcs_unstained_preprocessed <- function() {
    dir(fpath_fcs_unstained_comptrans, pattern = '[.]fcs$', full.names = TRUE)
  }
}

## Function: pre-processed stained samples directory ----

fpath_fcs_preprocessed <- function() {
  if (normalise) fpath_fcs_norm else fpath_fcs_comptrans
}

## Function: pre-processed unstained samples directory ----

fpath_fcs_unstained_preprocessed <- function() {
  fpath_fcs_unstained_comptrans
}

## Function: p-value to asterisks ----

sig_asterisks <- function(p) {
  return(
    ifelse(p<0.001, '***',
           ifelse(p<0.01, '**',
                  ifelse(p<0.05,
                         '*', '')))
  )
}

## Function: p-value to significance level ----

sig_level <- function(p, prefix = TRUE) {
  if (prefix) { s <- 'p' } else { s <- '' }
  return(
    ifelse(p<0.001, paste0(s, '<.001'),
           ifelse(p<0.01, paste0(s, '<.01'),
                  ifelse(p<0.05, paste0(s, '<.05'),
                         paste0(s, '=', round(p, 2)))))
  )
}

## Function: design linear modelling experiment ----

prep_experiment <- function(
  files,            # names of samples to include in experiment
  annotation,       # sample-level annotation
  fixed_effects,    # columns of `annotation` to model as fixed effects
  random_effects,   # columns of `annotation`to model as random effects
  y        = NULL,  # outcome values per sample (if only 1 compartment tested)
  olre     = FALSE, # whether to model observation-level random effects
  force_ls = FALSE  # whether to include LHS of formula even if `y` not given
) {
  
  ## Check whether all samples are included in the annotation data frame
  missing_files <- files[!files%in%annotation$FileName]
  n_missing <- length(missing_files)
  if (n_missing>1) {
    stop(
      n_missing, ifelse(n_missing>1, ' samples', ' sample'),
      ' not found in `annotation`')
  }
  
  ## Gather annotation of analysed samples
  ann <- annotation[
    match(files, annotation$FileName),
    c('FileName', unique(c(fixed_effects, random_effects))),
    drop = FALSE
  ]
  
  ## Include sample IDs as factor
  colnames(ann)[1] <- 'sample_id'
  ann[, 1] <- as.factor(ann[, 1])
  
  ## Convert character and logical columns to factors
  for (idx_col in colnames(ann)) {
    
    if (is.character(ann[, idx_col])) {
      
      ann[, idx_col] <- as.factor(ann[, idx_col])
    } else if (is.logical(ann[, idx_col])) {
      
      ann[, idx_col] <- factor(
        ann[, idx_col],
        levels = c(FALSE, TRUE)
      )
    }
  }
  
  ## Convert batch labels to factors
  if ('Batch' %in% colnames(ann)) {
    ann$Batch <- as.factor(ann$Batch)
  }
  
  ## Assume unpaired design
  experiment_info  <- data.frame(
    'group_id' = as.factor(1),
    ann
  )
  
  ## Include outcome values if provided
  if (!is.null(y)) {
    d <- cbind(
      'y' = y,
      experiment_info
    )
  } else {
    d <- experiment_info
  }
  
  ## Exclude samples with NAs (important: 1 NA value already causes exclusion!)
  na_rows  <- apply(X = d, MARGIN = 1, FUN = function(x) any(is.na(x)))
  na_files <- d$sample_id[na_rows]
  d        <- d[!na_rows, , drop = FALSE]
  
  ## Resolve RHS of model formula
  if (olre) { # observation-level random effects modelled
    
    random_effects <- c(random_effects, 'sample_id')
  }
  if (
    is.null(random_effects) ||
    length(random_effects) == 0
  ) { # no specific random effects given
    
    rhs <- paste(fixed_effects, collapse = ' + ')
  } else { # specific random effects given
    
    rhs <- paste0(
      paste0(fixed_effects, collapse = ' + '), ' + ',
      paste(paste0('(1|', random_effects, ')'), collapse = ' + ')
    )
  }
  
  ## Resolve LHS of model formula
  if (!is.null(y)) {
    lhs <- 'y'
  } else {
    lhs <- ''
  }
  
  ## Resolve full model formula
  formula <- stats::as.formula(paste0(lhs, '~', rhs))
  
  ## Create design matrix
  design <- stats::model.matrix(
    object = stats::as.formula(
      paste0('~', paste(fixed_effects, collapse = ' + '))
    ),
    data = d
  )
  
  ## Fill in LHS if needed
  if (is.null(y) && force_ls) {
    formula <- stats::as.formula(paste0('y ~', rhs))
  }
  
  list(
    'Data'          = d,
    'Formula'       = formula,
    'Design'        = design,
    'NA'            = na_files,
    'FixedEffects'  = fixed_effects,
    'RandomEffects' = random_effects
  )
}

## Function: gather inputs for DE modelling ----

prep_inputs <- function(
    fname_features,        # path to saved FlowSOM-derived features
    fname_sample_outliers, # path to saved sample-level outlier names
    annotation,            # sample-level annotation
    batch_remove           # batch(es) to remove from analysis
) {
  
  ## Load matrices of feature values
  f <- readRDS(fname_features)
  if (is.null(fname_sample_outliers)) {
    out_samples <- c()
  } else {
    out_samples <- readRDS(fname_sample_outliers)
  }
  mc_counts  <- f$metacluster_counts
  mc_perc    <- f$metacluster_percentages
  mc_mfi     <- f$metacluster_MFIs
  mc_percpos <- f$metacluster_percentages_positive
  
  ## Check that rows (samples) are aligned between the feature matrices
  if (
    nrow(mc_counts)!=nrow(mc_perc) ||
      nrow(mc_counts)!=nrow(mc_mfi) ||
      (!is.null(mc_percpos) && nrow(mc_counts)!=nrow(mc_percpos))
  ) {
    stop('Row (sample) counts are not the same between the feature matrices')
  }
  if (
    any(rownames(mc_counts)!=rownames(mc_perc)) ||
    any(rownames(mc_counts)!=rownames(mc_mfi)) ||
    (!is.null(mc_percpos) && any(rownames(mc_counts)!=rownames(mc_percpos)))
  ) {
    stop('Row (sample) names are not the same between the feature matrices')
  }
  
  ## Use file names without paths as row names
  samples <- basename(rownames(mc_counts))
  rownames(mc_counts) <- 
    rownames(mc_perc) <- 
    rownames(mc_mfi)  <- samples
  if (!is.null(mc_percpos)) {
    rownames(mc_percpos) <- samples
  }
  
  ## Initialise excluded samples
  excluded_samples <- list(
    'batch'    = c(),
    'outliers' = c()
  )
  
  ## Exclude unused batches
  if (
    !is.null(batch_remove) && length(batch_remove)>0
  ) {
    
    mask <- samples%in%annotation$FileName[annotation$Batch%in%batch_remove]
    
    if (any(mask)) {
      
      excluded_samples[['batch']] <- samples[mask]
      mc_counts <- mc_counts[!mask, , drop = FALSE]
      mc_perc   <- mc_perc[!mask, , drop = FALSE]
      mc_mfi    <- mc_mfi[!mask, , drop = FALSE]
      if (!is.null(mc_percpos)) {
        mc_percpos <- mc_percpos[!mask, , drop = FALSE]
      }
      samples <- samples[!mask]
    }
  }
  
  ## Exclude outliers
  if (
    length(out_samples)>0
  ) {
    
    mask <- samples%in%out_samples
    
    excluded_samples[['outliers']] <- samples[mask]
    mc_counts  <- mc_counts[!mask, , drop = FALSE]
    mc_perc    <- mc_perc[!mask, , drop = FALSE]
    mc_mfi     <- mc_mfi[!mask, , drop = FALSE]
    if (!is.null(mc_percpos)) {
      mc_percpos <- mc_percpos[!mask, , drop = FALSE]
    }
    samples <- samples[!mask]
  }
  
  ## Gather abundance and state feature matrices
  res <- list(
    'counts'      = mc_counts,
    'percentages' = mc_perc, 
    'MFIs'        = mc_mfi,
    'PercPos'     = mc_percpos
  )
  
  ## Note excluded samples
  attributes(res)[['excluded_samples']] <- excluded_samples
  
  res
}

## Function: match channel names using regex ----

match_channel_names <- function(
    ref, # flowCore::transformList channels to be matched
    real # flowCore::flowFrame channels to match against
) {
  
  ## Initialise 1-to-1 mapping
  matched_idcs <- rep(NA, length(ref))
  
  ## Iterate over queried channel names
  for (idx in seq_along(ref)) {
    
    one_ref <- ref[idx]
    
    ## Match by numeric part
    idcs <- which(
      stringr::str_extract(
        one_ref, '[[:digit:]]+'
      ) == stringr::str_extract(
        real, '[[:digit:]]+'
      )
    )
    
    if (length(idcs) != 1) {
      
      ## Match by alphanumeric part
      idcs <- which(
        stringr::str_extract(
          one_ref, '[A-Za-z][[0-9]]{3}'
        ) == stringr::str_extract(
          real,  '[A-Za-z][[0-9]]{3}'
        )
      )
    }
    
    if (length(idcs) != 1) {
      
      ## Match exactly
      idcs <- which(real == one_ref)
    }
    
    if (length(idcs) == 0) {
      
      idcs <- NA
    } else if (length(idcs) > 1) {
      
      stop(
        'transformList channel "', one_ref,
        '" was matched to multiple flowFrame channels:\n\t',
        paste(real[idcs], collapse = ',\n\r')
      )
    }
    
    matched_idcs[idx] <- idcs
  }
  
  matched <- real[matched_idcs]
  
  ## Return data frame of aligned channel names
  data.frame(
    'Reference' = ref,
    'Real'      = matched,
    'Index'     = matched_idcs
  )
  
}

## Function: extract confounders from DE results ----

get_confounders <- function(
  res,      # DE analysis results
  predictor # predictor of interest
) {
  
  if (attributes(res)$AnalysisType=='Differential Abundance') { # DA analysis
    
    if (!predictor%in%names(res$confounders_joint)) { # predictor not found
      
      return(c())
    } else  { # predictor found
      
      return(names(res$confounders_joint[[predictor]]))
    }
  } else { # DS-MFI or DS-Pheno analysis
    
    if (!predictor%in%names(res$confounders_main)) { # predictor not found
      
      return(c())
    } else { # predictor found
      
      return(names(res$confounders_main[[predictor]]))
    }
  }
}

## Function: extract predictors from DE results ----

get_predictors <- function(
  res # DE analysis results
) {
  
  if (attributes(res)$AnalysisType=='Differential Abundance') { # DA analysis
    
    return(names(res$joint))
  } else { # DS-MFI or DS-Pheno analysis
    
    return(names(res$main))
  }
}

## Function: parallelised flowFrame aggregation ----

ParallelAggregate <- function(
    fnames,                  # paths to FCS files to aggregate
    N,                       # target number of cells in total
    cols            = NULL,  # specific columns to include
    keep_file_order = FALSE, # whether to preserve files order
    keep_cell_order = FALSE, # whether to preserve cells order within files
    seed            = NULL,  # random seed for reproducibility
    cores           = parallel::detectCores(), # number of cores to use
    as_flowFrame    = TRUE,  # whether to return flowFrame (not data frame)
    fsom            = NULL,  # FlowSOM model to select metacluster to aggregate
    metacluster     = 1      # only metacluster to aggregate if `fsom` given
) {
  require(parallel)
  require(doParallel)
  require(foreach)
  require(snow)
  require(doSNOW)
  
  ## Make sure aggregate size is at least as large as the number of files
  n_files <- length(fnames)
  stopifnot(N>=n_files)
  
  ## Determine sample size per file
  if (is.null(fsom)) { # all cells to be sampled (no metacluster choice)
    
    ## Get cell counts per file
    s_full <- sapply(
      fnames,
      function(f) {
        as.integer(
          flowCore::read.FCSheader(f, keyword = '$TOT')[[1]][['$TOT']]
        )
      })
    if (is.null(N) || is.infinite(N)) { # no sub-sampling
      
      s_perfile <- s_full
      N <- sum(s_perfile)
    } else {  # sub-sampling
      
      ## Determine a priori average cell count taken per file
      s_generic  <- ceiling(N/n_files)
      
      ## Initialise real cell counts take per file
      s_perfile  <- rep(NA, times = n_files)
      
      ## Flag and count files with fewer cells than the a priori average count
      mask_small <- s_full<s_generic
      n_small    <- sum(mask_small)
      
      ## Flag and count files with more cells than the a priori average count
      mask_large <- !mask_small
      n_large    <- sum(mask_large)
      
      ## Take all cells from files with fewer cells than a priori
      s_perfile[mask_small] <- s_full[mask_small]
      
      ## Adjust number of cell taken from files with more cells than a priori
      s_adjusted            <- ceiling((N-sum(s_perfile, na.rm = TRUE))/n_large)
      s_perfile[mask_large] <- s_adjusted
    }
  } else { # only cells from chosen metacluster to be sampled
    
    ## Set (upper-bound) cell counts taken per file naively
    s_perfile <- rep(ceiling(N/n_files), times = n_files)
  }
  
  ## Set up parallel processing
  cl <- snow::makeCluster(cores-1)
  registerDoSNOW(cl)
  
  pb <- utils::txtProgressBar(
    min = 0, max = n_files, style = 3
  ) # progress bar
  opts <- list(
    'progress' = function(i) utils::setTxtProgressBar(pb, i)
  )
  
  ## Define sample ordering function
  if (keep_file_order) {
    f_sort <- function(res) res[order(res[, 'File']), , drop = FALSE]
  } else {
    f_sort <- NULL
  }
  
  ## Define cell ordering function
  if (keep_cell_order) {
    f_cell_idcs <- function(n, s) sort(sample.int(n, s))
  } else {
    f_cell_idcs <- function(n, s) sample.int(n, s)
  }
  
  ## Prepare inputs
  select_mc <- FALSE
  packages  <- c('flowCore')
  if (!is.null(fsom)) {
    
    select_mc <- TRUE
    packages  <- c(packages, 'FlowSOM')
    codes     <- fsom$map$codes
    meta      <- fsom$metaclustering
  }
  
  ## Generate data aggregate
  res <- foreach(
    i             = seq_along(fnames),
    .combine      = rbind,
    .inorder      = FALSE,
    .final        = f_sort,
    .packages     = packages,
    .options.snow = opts
  ) %dopar% {
    
    ## Load expression data
    d <- flowCore::read.FCS(fnames[i])[, cols, drop = FALSE]@exprs
    
    ## Subset by metacluster
    if (select_mc) {
      
      mc_idcs <- which(
        meta[FlowSOM:::MapDataToCodes(codes, d)[, 1]]==metacluster
      )
      d <- d[mc_idcs, , drop = FALSE]
    }
    
    ## No sub-sampling if 0-1 cells remain
    nd <- nrow(d)
    if (nd<2) {
      return(d)
    }
    
    ## Sub-sample cells if needed
    s <- min(s_perfile[i], nd)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    ## Resolve cell ordering
    idcs <- f_cell_idcs(nd, s)
    d    <- d[idcs, , drop = FALSE]
    
    ## Add file of origin and original cell index as columns
    cbind(
      d,
      'File' = i,
      'OriginalIndex' = idcs
    )
  }
  close(pb)
  snow::stopCluster(cl)
  
  ## Return data frame if requested
  if (!as_flowFrame) {
    return(res)
  }
  
  ## Load template flowFrame object
  ff <- flowCore::read.FCS(fnames[1], which.lines = 1)[, cols, drop = FALSE]
  
  ## Fill in aggregate expression data
  ff@exprs <- res[, cols, drop = FALSE]
  
  ## Add file of origin and original cell index as columns
  flowCore::fr_append_cols(
    fr   = ff,
    cols = res[, c('File', 'OriginalIndex')]
  )
}

## Function: resolve which samples were used per test ----

samples_per_test <- function(
  res,               # DE analysis results
  predictor,         # predictor of interest
  confounder = NULL, # confounder of interest
  comps      = NULL  # compartments of interest
) {
  
  ## Resolve analysis type
  type <- attributes(res)$AnalysisType
  
  ## Resolve whether confounder specified
  wconf <- !is.null(confounder) && confounder!='none'
  
  ## Resolve input samples
  samples   <- attributes(res)$InputSamples
  n_samples <- attributes(res)$NInputSamples
  
  ## Resolve compartments of interest if not specified
  if (is.null(comps)) {
    if (type=='Differential Abundance') {
      
      if (wconf) {
        comps <- unique(
          res$confounders_joint[[predictor]][[confounder]]$Compartment
        )
      } else {
        comps <- unique(
          res$joint[[predictor]]$Compartment
        )
      }
    } else { # DS
      
      if (wconf) {
        comps <- unique(
          res$confounders_main[[predictor]][[confounder]]$Compartment
        )
      } else {
        comps <- unique(
          res$main[[predictor]]$Compartment
        )
      }
    }
  } else {
    comps <- unique(comps)
  }
  n_comps <- length(comps)
  
  ## Identify samples w/ missing annotation (globally) or outcome (per test)
  if (type=='Differential Abundance') {
    
    na_annotation <- as.character(
      if (wconf) {
        attributes(
          res$confounders_joint[[predictor]][[confounder]]
        )$na_annotation
      } else {
        attributes(
          res$joint[[predictor]]
        )$na_annotation
      }
    )
    na_outcome <- c()
  } else { # DS
    
    na_annotation <- as.character(
      if (wconf) {
        res$confounders_na_annotation[[predictor]][[confounder]]
      } else {
        res$na_annotation[[predictor]]
      }
    )
    na_outcome <-
      if (wconf) {
        res$confounders_na_outcome[[predictor]][[confounder]]
      } else {
        res$na_outcome[[predictor]]
      }
  }
  
  ## Generate a Boolean matrix for inclusion of samples per text (compartment)
  if (type=='Differential Abundance') {
    
    out <- matrix(TRUE, nrow = n_samples, ncol = n_comps)
    rownames(out) <- samples
    colnames(out) <- comps
    out[na_annotation, ] <- FALSE
  } else {
    
    excl <- # excluded samples per test
      lapply(
        na_outcome[comps],
        function(x) {
          unique(c(x, na_annotation))
        }
      )
    
    out <- do.call(cbind, lapply(excl, function(x) !samples%in%x))
    rownames(out) <- samples
  }
  
  out
}

## Function: compute robustness rates for DS results ----

robustness_rate_per_test <- function(
  res,                # DS analysis results
  robustness,         # list of robustness matrices per signal filter strength
  samples_mat = NULL, # matrix created by `samples_per_test` if pre-computed
  predictor   = NULL, # biological predictor of interest
  confounder  = NULL, # biological confounder of interest
  comps       = NULL, # compartments of interest
  tidy        = FALSE # whether to return as list of data frames per cut-off
) {

  ## Check if predictor of samples-per-test matrix provided
  if (is.null(samples_mat)) {
    if (is.null(predictor)) {
      stop('If `samples_mat` not given, `predictor` must be specified')
    }
  }
  
  ## Resolve samples-per-test matrix
  if (is.null(samples_mat)) {
    samples_mat <- samples_per_test(
      res        = res,
      predictor  = predictor,
      confounder = confounder,
      comps      = comps
    )
  } else {
    if (!is.null(comps)) {
      samples_mat <- samples_mat[, comps, drop = FALSE]
    }
  }
  
  ## Resolve input samples and tests (~compartments)
  samples   <- rownames(samples_mat)
  n_samples <- length(samples)
  comps     <- colnames(samples_mat)
  n_comps   <- length(comps)
  
  ## Compute robustness rates per signal filter strength level
  out <- lapply( # iterate over signal filter strength
    robustness,
    function(r) {
      
      ## Filter robustness matricex
      robustness_mat <- r[
        samples, comps, drop = FALSE
      ]
      robustness_mat[!samples_mat] <- NA
      
      ## Count number of input samples per test
      ns <- colSums(samples_mat)
      
      ## Count number of input samples with robust signal per test
      rs <- colSums(robustness_mat, na.rm = TRUE)
      
      ## Compute robustness rate per test
      out <- rs/ns
      names(out) <- comps
      
      out
    }
  )
  
  if (tidy) {
    
    mcs     <- gsub('[ ].*$', '', names(out[[1]]))
    markers <- gsub('^MC[0-9]+[ ]', '', names(out[[1]]))
    
    out <- lapply(
      seq_along(out),
      function(i) {
        data.frame(
          'Metacluster' = mcs,
          'Marker'      = markers,
          'Compartment' = paste0(mcs, ' ', markers),
          'Rate'        = as.vector(out[[i]])
        )
      }
    )
  }
  
  out
}

## Function: collect DE results in tabular form ----

de_results <- function(
    res,                   # DA, DS-MFI or DS-Pheno results object
    predictor,             # biological predictor of interest
    confounder    = NULL,  # biological confounder of interest
    metacluster   = NULL,  # metacluster of interest (if not all)
    marker        = NULL,  # state marker of interest (if not all)
    robustness    = NULL,  # robustness info
    state_markers = NULL,  # all state markers
    wide_format   = FALSE
) {
  
  ## Resolve analysis type
  type       <- attributes(res)$AnalysisType
  
  ## Resolve whether confounder specified
  wconf <- !is.null(confounder) && confounder!='none'
  
  ## Resolve paths to results in `res` based on analysis type
  if (type=='Differential Abundance') {
    
    slot_stat   <- ifelse(wconf, 'confounders_joint', 'joint') # stats slot
    slot_change <- 'logFC'            # effect coefficient column
    type_change <- 'log2 fold change' # full name of effect coefficient
    analysis_type_short <- 'DA'       # short name of analysis type
  } else if (type=='Differential State (MFI)') {
    
    slot_stat <- ifelse(wconf, 'confounders_main', 'main') # stats slot
    slot_change <- 'logFC'            # effect coefficient column
    type_change <- 'log2 fold change' # full name of effect coefficient
    analysis_type_short <- 'DS-MFI'   # short name of analysis type
  } else if (type=='Differential State (Phenopositivity)') {
    
    slot_stat <- ifelse(wconf, 'confounders_main', 'main') # stats slot
    slot_change <- 'logodds'          # effect coefficient column
    type_change <- 'log2 odds'        # full name of effect coefficient
    analysis_type_short <- 'DS-Pheno' # short name of analysis type
  }
  
  ## Extract data frame with fitted model parameters
  if (wconf) { # confounder specified
    
    st <- res[[slot_stat]][[predictor]][[confounder]]
  } else { # no confounder specified
    
    st <- res[[slot_stat]][[predictor]]
  }
  
  ## Add metacluster and marker columns to the stats data frame
  if (type=='Differential Abundance') { # only metacluster relevant
    
    st$Metacluster <- st$Compartment
    st$Marker      <- NA
    if (!is.null(metacluster)) {
      st <- st[st$Metacluster==metacluster, , drop = FALSE]
    }
  } else { # both metacluster and marker relevant
    
    st$Metacluster <- gsub('[ ].*$', '', st$Compartment)
    st$Marker      <- gsub('^MC[0-9]+[ ]', '', st$Compartment)
    if (!is.null(metacluster)) {
      st <- st[st$Metacluster==metacluster, , drop = FALSE]
    }
    if (!is.null(marker)) {
      st <- st[st$Marker==marker, , drop = FALSE]
    }
    if (!is.null(state_markers)) {
      st <- st[st$Marker%in%state_markers, , drop = FALSE]
    }
  }
  
  ## Generate matrix indicating samples included in each test
  samples_mat <- samples_per_test(
    res        = res,
    predictor  = predictor,
    confounder = confounder,
    comps      = st$Compartment
  )
  
  ## Note numbers of included/excluded input samples to each test
  st$nIncludedSamples <- colSums(samples_mat)
  st$nExcludedSamples <- colSums(!samples_mat)
  
  ## Compute signal robustness rates per test (~compartment)
  if (
    type!='Differential Abundance' &&
    !is.null(robustness)
  ) {
    
    rob <- robustness_rate_per_test(
      res         = res,
      robustness  = robustness,
      samples_mat = samples_mat
    )
    st$RobustnessAt1MAD <- rob[[1]]
    st$RobustnessAt2MAD <- rob[[2]]
    st$RobustnessAt3MAD <- rob[[3]]
    st$RobustnessAt4MAD <- rob[[4]]
    st$RobustnessAt5MAD <- rob[[5]]
    
  } else {
    st$RobustnessAt1MAD <- NA
    st$RobustnessAt2MAD <- NA
    st$RobustnessAt3MAD <- NA
    st$RobustnessAt4MAD <- NA
    st$RobustnessAt5MAD <- NA
  }
  
  ## Extract p-values and effect coefficients
  p_pred_raw <- st$PValue
  p_pred_adj <- st$AdjPVal
  ch_pred    <- st[, slot_change]
  p_conf_raw <- NA
  p_conf_adj <- NA
  ch_conf    <- NA
  if (wconf) {
    p_conf_raw <- st$PValueConfounder
    p_conf_adj <- st$AdjPValConfounder
    ch_conf    <- st[, paste0(slot_change, 'Confounder')]
  }
  
  str_conf <- confounder
  if (!wconf) {
    confounder <- NA
    str_conf <- NA
  }
  
  ## Resolve identities of the model(s)
  name <- paste0(
    'STATS_', analysis_type_short,
    '_Pred', predictor,
    '_Conf', confounder
  )
  if (!is.null(metacluster)) {
    name <- paste0(name, '_', metacluster)
  }
  if (type!='Differential Abundance' && !is.null(marker)) {
    name <- paste0(name, '_Marker', marker)
  }
  
  ## Create output data frame
  d <- as.data.frame(
    cbind(
      c(
        'Analysis type',
        'Metacluster',
        'Marker',
        'Biological predictor',
        'Biological confounder',
        'Coefficient type',
        'Predictor coefficient',
        'Confounder coefficient',
        'Predictor p-value (raw)',
        'Predictor p-value (adjusted)',
        'Confounder p-value (raw)',
        'Confounder p-value (adjusted)',
        'Number of included samples',
        'Number of excluded samples',
        'Robustness rate at 1MAD',
        'Robustness rate at 2MAD',
        'Robustness rate at 3MAD',
        'Robustness rate at 4MAD',
        'Robustness rate at 5MAD'
      ),
      t(data.frame(
        type,
        st$Metacluster,
        st$Marker,
        predictor,
        confounder,
        type_change,
        ch_pred,
        ch_conf,
        p_pred_raw,
        p_pred_adj,
        p_conf_raw,
        p_conf_adj,
        st$nIncludedSamples,
        st$nExcludedSamples,
        st$RobustnessAt1MAD,
        st$RobustnessAt2MAD,
        st$RobustnessAt3MAD,
        st$RobustnessAt4MAD,
        st$RobustnessAt5MAD
      ))
    )
  )
  cn <- 'Name'
  if (ncol(d)==2) {
    cn <- c(cn, 'Value')
  } else {
    cn <- c(
      cn,
      paste0(
        'Value',
        formatC(
          seq_len(ncol(d)-1),
          width  = 5,
          format = 'd',
          flag   = '0'
        )
      )
    )
  }
  colnames(d) <- cn
  rownames(d) <- NULL
  
  ## Convert to long format if requested
  if (!wide_format) {
    
    res <- as.data.frame(
      t(d)[2:ncol(d), , drop = FALSE]
    )
    colnames(res) <- d$Name
    rownames(res) <- NULL
    res[, 'Predictor coefficient'] <-
      as.numeric(res[, 'Predictor coefficient'])
    res[, 'Confounder coefficient'] <-
      as.numeric(res[, 'Confounder coefficient'])
    res[, 'Predictor p-value (raw)'] <-
      as.numeric(res[, 'Predictor p-value (raw)'])
    res[, 'Predictor p-value (adjusted)'] <-
      as.numeric(res[, 'Predictor p-value (adjusted)'])
    res[, 'Confounder p-value (raw)'] <-
      as.numeric(res[, 'Confounder p-value (raw)'])
    res[, 'Confounder p-value (adjusted)'] <-
      as.numeric(res[, 'Confounder p-value (adjusted)'])
    res[, 'Number of excluded samples'] <-
      as.numeric(res[, 'Number of excluded samples'])
    res[, 'Robustness rate at 1MAD'] <-
      as.numeric(res[, 'Robustness rate at 1MAD'])
    res[, 'Robustness rate at 2MAD'] <-
      as.numeric(res[, 'Robustness rate at 2MAD'])
    res[, 'Robustness rate at 3MAD'] <-
      as.numeric(res[, 'Robustness rate at 3MAD'])
    res[, 'Robustness rate at 4MAD'] <-
      as.numeric(res[, 'Robustness rate at 4MAD'])
    res[, 'Robustness rate at 5MAD'] <-
      as.numeric(res[, 'Robustness rate at 5MAD'])
    d <- res
  }
  
  ## Specify name for CSV export
  attributes(d)$Name <- name
  
  d
}

## Define palette of distinctive colours ----

pal <- c(
  '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', '#006FA6', '#A30059',
  '#7A4900', '#0000A6', '#63FFAC', '#B79762', '#004D43', '#8FB0FF',
  '#997D87', '#5A0007', '#809693', '#1B4400', '#4FC601', '#3B5DFF',
  '#4A3B53', '#FF2F80', '#61615A', '#BA0900', '#6B7900', '#00C2A0',
  '#FFAA92', '#B903AA', '#D16100', '#DDEFFF', '#000035', '#7B4F4B',
  '#A1C299', '#300018', '#0AA6D8', '#013349', '#00846F', '#372101',
  '#FFB500', '#C2FFED', '#A079BF', '#CC0744', '#C0B9B2', '#C2FF99',
  '#001E09', '#00489C', '#6F0062', '#0CBD66', '#EEC3FF', '#456D75',
  '#B77B68', '#7A87A1', '#788D66', '#885578', '#FAD09F', '#FF8A9A',
  '#D157A0', '#BEC459', '#456648', '#0086ED', '#886F4C', '#34362D', 
  '#B4A8BD', '#00A6AA', '#452C2C', '#636375', '#A3C8C9', '#FF913F',
  '#938A81', '#575329', '#00FECF', '#B05B6F', '#8CD0FF', '#3B9700',
  '#04F757', '#C8A1A1', '#1E6E00', '#7900D7', '#A77500', '#6367A9',
  '#A05837', '#6B002C', '#772600', '#D790FF', '#9B9700', '#549E79',
  '#FFF69F', '#201625', '#72418F', '#BC23FF', '#99ADC0', '#3A2465',
  '#922329', '#5B4534', '#FDE8DC', '#404E55', '#0089A3', '#CB7E98',
  '#A4E804', '#324E72', '#6A3A4C'
)



