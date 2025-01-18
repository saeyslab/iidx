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

## iidx user-level module: 02_FlowSOMFeatureExtraction.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Uses FlowSOM to cluster data into subsets, extract salient
# features and diagnose remaining batch effects.


message(
  '## iidx: 02_FlowSOMFeatureExtraction.R STARTED ',
  as.character(round(Sys.time()))
)

## Get parameters and auxiliary functions ----

input_files_needed <- TRUE
source('99_SetUp.R')
source('02a_BatchEffectDiagnostics.R')

## Create results directory ----

if (!file.exists(fpath_res02)) {
  dir.create(fpath_res02)
}

## Create directory for feature matrices
if (!file.exists(fpath_features)) {
  dir.create(fpath_features)
}

## Construct large aggregate for feature extraction ----

message('Aggregating expression data')

## Sample from pre-processed (and possibly normalised) FCS files
set.seed(1); agg <-
  FlowSOM::AggregateFlowFrames(
    fileNames = fnames_fcs_preprocessed(),
    cTotal    = 6e6,
    channels  = names(channels) # keep all channels used in analysis
  )

## Save as FCS file
flowCore::write.FCS(agg, fname_agg)

## Train FlowSOM model for feature extraction -----

message('Training FlowSOM model using lineage markers')
t_fsom <- system.time(
  fsom <- FlowSOM::FlowSOM(
    input     = agg,
    colsToUse = names(channels)[idcs_channels_lineage],
    xdim      = xdim,
    ydim      = ydim,
    scale     = FALSE, # (we assume signal ranges are comparable)
    nClus     = n_metaclusters,
    seed      = 1
  )
)

## Make sure channel names available
i     <- !is.na(channels[names(fsom$prettyColnames)])
pc    <- fsom$prettyColnames
pc[i] <- channels[names(fsom$prettyColnames)][i]
fsom$prettyColnames <- pc

## Save time needed to train FlowSOM model
fsom$time <- t_fsom

## Save FlowSOM model and metaclustering ----

saveRDS(fsom, fname_fsom)
mcs <- FlowSOM::GetMetaclusters(fsom) # metacluster label per row of `agg`
saveRDS(mcs, fname_metaclustering)

## Create FlowSOMmary PDF to check quality of clustering
FlowSOM::FlowSOMmary(fsom, plotFile = fname_fsommary)

## Resolve features to extract ----

feature_types <- c(
  'counts',      # absolute abundance
  'percentages', # relative abundance
  'MFIs'         # state (by signal)
)

## If thresholds specified, compute phenopositivity rates
if (!is.null(thresholds) &&
    length(thresholds)>0
) {
  feature_types <- c(
    feature_types,
    'percentages_positive' # state (by phenopositivity)
  )
} 

## Extract features at metacluster level ----

## Prepare extraction for up to `feature_mapping_batch_size` files at a time
fnames <- fnames_fcs_preprocessed()
n_fcs  <- length(fnames)
s      <- seq(
  from = 1,
  to   = n_fcs,
  by   = feature_mapping_batch_size
)
idcs_start <- seq(
  from = 1,
  to   = n_fcs,
  by   = feature_mapping_batch_size
)
idcs_end <- c(idcs_start[-1]-1, n_fcs)
n_mb <- length(idcs_start) # number of iterations

## Determine whether to get phenopositivity rates ('percentages positive')
phenopos <- !is.null(thresholds) && length(thresholds)>0

if (run_parallel) { # multi-threading enabled
  
  message(
    'Running FlowSOM feature extraction'
  )
  
  ## Set up parallel processing
  cores <- detectCores()
  cl    <- makeCluster(cores[1]-1)
  registerDoSNOW(cl)
  
  pb <- utils::txtProgressBar(min = 0, max = n_mb, style = 3) # progress bar
  
  opts <- list(
    'progress' = function(i) utils::setTxtProgressBar(pb, i)
  )
  
  ## Map all analysis samples
  tmp <- foreach(
    i             = seq_len(n_mb),
    .combine      = c,
    .options.snow = opts
  ) %dopar% {
    
    i_start <- idcs_start[i]
    i_end   <- idcs_end[i]
    f       <- fnames[i_start:i_end]
    features <- FlowSOM::GetFeatures(
      fsom,
      files            = f,
      level            = 'metaclusters',
      type             = feature_types,
      MFI              = names(channels)[idcs_channels_state],
      positive_cutoffs = if (phenopos) { thresholds } else { NULL },
      filenames        = basename(f),
      silent           = TRUE
    )
    saveRDS(
      features,
      paste0(
        fname_features_prefix_postnorm, i_start, '-', i_end, '.RDS'
      )
    )
    NULL
  }
  close(pb)
  stopCluster(cl)
  
  if (normalise) { # if normalisation used, extract features for pre-norm also
    
    message(
      'Running FlowSOM feature extraction for data before normalisation'
    )
    
    ## Set up parallel processing
    cores <- detectCores()
    cl    <- makeCluster(cores[1]-1)
    registerDoSNOW(cl)
    
    pb <- utils::txtProgressBar(min = 0, max = n_mb, style = 3) # progress bar
    
    opts <- list(
      'progress' = function(i) utils::setTxtProgressBar(pb, i)
    )
    
    ## Map all pre-normalisation samples
    tmp <- foreach(
      i             = seq_len(n_mb),
      .combine      = c,
      .options.snow = opts
    ) %dopar% {
      
      i_start <- idcs_start[i]
      i_end   <- idcs_end[i]
      f       <- fnames_fcs_comptrans[i_start:i_end]
      features <- FlowSOM::GetFeatures(
        fsom,
        files            = f,
        level            = 'metaclusters',
        type             = feature_types,
        MFI              = names(channels)[idcs_channels_state],
        positive_cutoffs = if (phenopos) { thresholds } else { NULL },
        filenames        = basename(f),
        silent           = TRUE
      )
      saveRDS(
        features,
        paste0(
          fname_features_prefix_prenorm, i_start, '-', i_end, '.RDS'
        )
      )
      NULL
    }
    close(pb)
    stopCluster(cl)
  }
} else {
  
  message(
    'Running FlowSOM feature extraction'
  )
  
  pb <- utils::txtProgressBar(min = 0, max = n_mb, style = 3) # progress bar
  
  ## Map all analysis samples
  for (i in seq_len(n_mb)) {
    
    i_start <- idcs_start[i]
    i_end   <- idcs_end[i]
    f       <- fnames[i_start:i_end]
    features <- FlowSOM::GetFeatures(
      fsom,
      files            = f,
      level            = 'metaclusters',
      type             = feature_types,
      MFI              = names(channels)[idcs_channels_state],
      positive_cutoffs = if (phenopos) { thresholds } else { NULL },
      filenames        = basename(f),
      silent           = TRUE
    )
    saveRDS(
      features,
      paste0(
        fname_features_prefix_postnorm, i_start, '-', i_end, '.RDS'
      )
    )
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  if (normalise) {
    
    message(
      'Running FlowSOM feature extraction for data before normalisation'
    )
    
    pb <- utils::txtProgressBar(min = 0, max = n_mb, style = 3) # progress bar
    
    ## Map all pre-normalisation samples
    for (i in seq_len(n_mb)) {
      
      i_start <- idcs_start[i]
      i_end   <- idcs_end[i]
      f       <- fnames_fcs_comptrans[i_start:i_end]
      features <- FlowSOM::GetFeatures(
        fsom,
        files            = f,
        level            = 'metaclusters',
        type             = feature_types,
        MFI              = names(channels)[idcs_channels_state],
        positive_cutoffs = if (phenopos) { thresholds } else { NULL },
        filenames        = basename(f),
        silent           = TRUE
      )
      saveRDS(
        features,
        paste0(
          fname_features_prefix_prenorm, i_start, '-', i_end, '.RDS'
        )
      )
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
  }
}

## Gather full feature matrices ----

message('Saving extracted features as full feature matrices')

## Auxiliary functions to import and bind corresponding matrices of values...
rbind_by_slot <- function(fnames, slot) {
  do.call(
    rbind,
    lapply(fnames,  function(x) readRDS(x)[[slot]])
  )
}
rbind_by_all_slots <- function(fnames) {
  list(
    'metacluster_counts' =
      rbind_by_slot(fnames, 'metacluster_counts'),
    'metacluster_percentages' =
      rbind_by_slot(fnames, 'metacluster_percentages'),
    'metacluster_MFIs' =
      rbind_by_slot(fnames, 'metacluster_MFIs'),
    'metacluster_percentages_positive' =
      rbind_by_slot(fnames, 'metacluster_percentages_positive'),
    'cluster_counts' = # (currently NULL)
      rbind_by_slot(fnames, 'cluster_counts'),
    'cluster_percentages' = # (currently NULL)
      rbind_by_slot(fnames, 'cluster_percentages'),
    'cluster_MFIs' = # (currently NULL)
      rbind_by_slot(fnames, 'cluster_MFIs'),
    'cluster_percentages_positive' = # (currently NULL)
      rbind_by_slot(fnames, 'cluster_percentages_positive')
  )
}

fnames_pre   <- NULL
features_pre <- NULL

## If normalisation was not applied, features are saved as 'PostNorm' and no...
## ...'PreNorm' exists
fnames_post   <- paste0(
  fname_features_prefix_postnorm, idcs_start, '-', idcs_end, '.RDS'
)
features_post <- rbind_by_all_slots(fnames_post)
saveRDS(features_post, fname_features_postnorm)

if (normalise) { # if normalisation applied, gather pre-norm values separately
  fnames_pre   <- paste0(
    fname_features_prefix_prenorm, idcs_start, '-', idcs_end, '.RDS'
  )
  features_pre <- rbind_by_all_slots(fnames_pre)
  saveRDS(features_pre, fname_features_prenorm)
}

## Compute noise levels from background signal ----

if (unstained_samples_available) { # pre-processed unstained data exists
  
  message('Computing background signal in each batch using unstained samples')
  
  ## Resolve pre-processed unstained FCS file names and corresponding batches
  fnames  <- fnames_fcs_unstained_preprocessed()
  batches <- gsub('[.]fcs$', '', gsub('^.*Batch', '', fnames))
  
  ## Compute median and MAD of background signal per marker
  noise <- do.call(
    rbind,
    lapply(
      seq_along(batches),
      function(i) {
        batch <- batches[i]
        fname <- fnames[i]
        ff    <- flowCore::read.FCS(fname)
        ff    <- ff[, names(channels)[
          c(idcs_channels_lineage, idcs_channels_state)
        ]]
        data.frame(
          'Batch'   = batch,
          'Marker'  = channels[c(idcs_channels_lineage, idcs_channels_state)],
          'Median'  = apply(ff@exprs, 2, median),
          'MAD'     = apply(ff@exprs, 2, mad)
        )
      }
    )
  )
  noise$Batch  <- as.factor(noise$Batch)
  noise$Marker <- as.factor(noise$Marker)
  
  ## Save noise stats
  saveRDS(noise, fname_noise)
}

## Generate UMAPs for batch effect diagnostics ----

## (If normalisation was run, embeddings of feature vectors from both...
## ...pre-norm and post-norm data are generated and plotted side-by-side to...
## ...verify that batch effect was successfully reduced; if normalisation was...
## ...not run, only one embedding is generated)

## Resolve QC files to highlight them in embeddings
fnames_qc <- file.path(
  fpath_fcs_norm, annotation$FileName[annotation$QC]
)
if (length(fnames_qc)==0) {
  fnames_qc <- NULL
}

message('Visualising batch effects in lineage markers')

viz_norm(
  features_pre =
    if (normalise) {
      readRDS(fname_features_prenorm)$metacluster_percentages
    } else {
      readRDS(fname_features_postnorm)$metacluster_percentages
    },
  fname_pdf  = fname_viz_norm_lineage,
  annotation = annotation,
  pal        = pal,
  fnames_qc  = fnames_qc,
  features_post =
    if (normalise) {
      readRDS(fname_features_postnorm)$metacluster_percentages
    } else {
      NULL
    },
  seed = 1
)
  
message('Visualising batch effects in state markers using MFI values')

f_mfi <- handle_na_values( # handle potential NAs for UMAP first
  f_pre = if (normalise) {
      readRDS(fname_features_prenorm)$metacluster_MFIs
    } else {
      readRDS(fname_features_postnorm)$metacluster_MFIs
    },
  f_post = if (normalise) {
      readRDS(fname_features_postnorm)$metacluster_MFIs
    } else {
      NULL
    }
)
f_pre  <- f_mfi$pre
f_post <- f_mfi$post

viz_norm(
  features_pre  = f_pre,
  fname_pdf     = fname_viz_norm_state_mfi,
  annotation    = annotation,
  pal           = pal,
  fnames_qc     = fnames_qc,
  features_post = f_post
)

message('Visualising batch effects in state markers using phenopositivities')

f_pheno <- handle_na_values( # handle potential NAs for UMAP first
  f_pre = if (normalise) {
    readRDS(fname_features_prenorm)$metacluster_percentages_positive
  } else {
    readRDS(fname_features_postnorm)$metacluster_percentages_positive
  },
  f_post = if (normalise) {
    readRDS(fname_features_postnorm)$metacluster_percentages_positive
  } else {
    NULL
  }
)
f_pre  <- f_pheno$pre
f_post <- f_pheno$post

viz_norm(
  features_pre  = f_pre,
  fname_pdf     = fname_viz_norm_state_pheno,
  annotation    = annotation,
  pal           = pal,
  fnames_qc     = fnames_qc,
  features_post = f_post
)

message(
  '## iidx: 02_FlowSOMFeatureExtraction.R FINISHED ',
  as.character(round(Sys.time()))
)
