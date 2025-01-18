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

## iidx internal module: 01a_Normalisation.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines functions for training, deployment and evaluation of 
# CytoNorm models.


## Function: match FCS files by sample ID ----

sample_match <- function(
    ref, # sample identifiers
    real # real file names
) {
  
  ## Ignore parent paths
  ref  <- basename(ref)
  real <- basename(real)
  
  ## Initialise results
  unmatched_idcs <- c()
  matched_idcs   <- rep(NA, length(ref))
  
  ## Iterate over sample IDs
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
    
    ## If ambiguous, try exact match
    if (length(idcs) != 1) {
      
      idcs <- grep(one_ref, real, fixed = TRUE)
    }
    
    ## If still ambiguous, mark as unmatched
    if (length(idcs) != 1) {
      
      idcs <- NA
      unmatched_idcs <- c(unmatched_idcs, one_ref)
    }
    
    matched_idcs[idx] <- idcs
  }
  
  matched <- real[matched_idcs]
  
  data.frame(
    'Sample' = ref,
    'Real'   = matched,
    'Index'  = matched_idcs
  )
}

## Function: load marker signal per batch ----

load_batch_aggregate_data_for_marker <- function(
    channel,                # channel name
    marker,                 # corresponding marker name
    fpath_agg_batchwise,    # path to batch aggregates
    batch_ref,              # reference batch number(s)
    batch_norm,             # target batch number(s)
    annotation,             # sample-level annotation
    include_postnorm = TRUE # whether to include post-normalisation data
) {
  
  ## Initialise reference, pre-normalisation and post-normalisation data
  data_ref      <- vector(mode = 'list', length = length(batch_ref))
  data_prenorm  <- vector(mode = 'list', length = length(batch_norm))
  data_postnorm <- vector(mode = 'list', length = length(batch_norm))
  
  ## Collect signal for reference
  counter <- 1
  for (idx_batch in batch_ref) {
    
    str_batch <- formatC(idx_batch, width = 2, format = 'd', flag = '0')
    fname <- file.path(
      fpath_agg_batchwise,
      paste0('Aggregate_Ref_Batch', str_batch, '.fcs')
    )
    
    data_ref[[counter]] <- data.frame(
      'Signal'     = flowCore::read.FCS(fname)@exprs[, channel],
      'BatchFlag'  = paste0('Ref', str_batch),
      'BatchIndex' = idx_batch,
      'Marker'     = marker,
      'Type'       = 'Reference'
    )
    counter <- counter+1
  }
  
  ## Collect signal for target
  counter <- 1
  for (idx_batch in batch_norm) {
    
    str_batch <- formatC(idx_batch, width = 2, format = 'd', flag = '0')
    
    fname_pre <- file.path(
      fpath_agg_batchwise,
      paste0('Aggregate_PreNorm_Batch', str_batch, '.fcs')
    )
    fname_post <- file.path(
      fpath_agg_batchwise,
      paste0('Aggregate_PostNorm_Batch', str_batch, '.fcs')
    )
    
    data_prenorm[[counter]] <- data.frame(
      'Signal'     = flowCore::read.FCS(fname_pre)@exprs[, channel],
      'BatchFlag'  = paste0('Pre', str_batch),
      'BatchIndex' = idx_batch,
      'Marker'     = marker,
      'Type'       = 'Pre-norm'
    )
    if (include_postnorm) {
      data_postnorm[[counter]] <- data.frame(
        'Signal'     = flowCore::read.FCS(fname_post)@exprs[, channel],
        'BatchFlag'  = paste0('Post', str_batch),
        'BatchIndex' = idx_batch,
        'Marker'     = marker,
        'Type'       = 'Post-norm'
      )
    }
    counter <- counter + 1
  }
  
  ## Row-bind into a single data frame
  if (include_postnorm) {
    data <- rbind(
      do.call(rbind, data_ref),
      do.call(rbind, data_prenorm),
      do.call(rbind, data_postnorm)
    )
  } else {
    data <- rbind(
      do.call(rbind, data_ref),
      do.call(rbind, data_prenorm)
    )
  }
  rm(data_ref)
  rm(data_prenorm)
  rm(data_postnorm)
  
  ## Adjust factor levels
  b <- unique(data$BatchFlag)
  b <- b[
    c(
      grep('^Ref', b),
      as.vector(t(
        cbind(grep('^Pre', b), grep('^Post', b))
      ))
    )
  ]
  data$BatchFlag <- ordered(data$BatchFlag, levels = b)
  data$Type      <- as.factor(data$Type)
  data$Type      <- relevel(data$Type, ref = 'Pre-norm')
  data$Type      <- relevel(data$Type, ref = 'Reference')
  data$RefOrTarget <- as.character(data$Type)
  data$RefOrTarget[
    data$RefOrTarget %in% c('Pre-norm', 'Post-norm')
  ] <- 'Target'
  data$RefOrTarget <- ordered(
    data$RefOrTarget,
    levels = c('Reference', 'Target')
  )
  
  ## Attach colour palettes for plotting
  pal_ref  <- grDevices::colorRampPalette(
    colors = c('#C6DBEF', '#08306B')
  )(length(batch_ref))
  pal_pre  <- grDevices::colorRampPalette(
    colors = c('#FCBBA1', '#67000D')
  )(length(batch_norm))
  pal_post <- grDevices::colorRampPalette(
    colors = c('#C7E9C0', '#00441B')
  )(length(batch_norm))
  pal_norm <- as.vector(
    t(cbind(pal_pre, pal_post))
  )
  cc <- c(pal_ref, pal_norm)
  attributes(data)[['palette']]           <- cc
  attributes(data)[['palette_reference']] <- pal_ref
  attributes(data)[['palette_prenorm']]   <- pal_pre
  attributes(data)[['palette_postnorm']]  <- pal_post
  
  data
}

## Function: compute EMDs between batches ----

cross_batch_emd <- function(
    data,                   # from `load_batch_aggregate_data_for_marker`
    batch_ref,              # reference batch number(s)
    emd_breaks,             # breaks for binning signal
    include_postnorm = TRUE # whether to include post-normalisation data
) {
  
  ## Define EMD computation via histograms
  emd_between_batches <- function(
    batch1, type1, # type ~ ref/pre/post
    batch2, type2
  ) {
    ## Bin batch-1 signal data
    y1 <- data$Signal[
      data$BatchIndex == batch1 & data$Type == type1
    ]
    y1 <- y1[
      y1 >= min(emd_breaks) & y1 <= max(emd_breaks)
    ]
    y1 <- graphics::hist(
      y1,
      breaks = emd_breaks,
      plot   = FALSE
    )$counts
    
    ## Bin batch-2 signal data
    y2 <- data$Signal[
      data$BatchIndex == batch2 & data$Type == type2
    ]
    y2 <- y2[
      y2 >= min(emd_breaks) & y2 <= max(emd_breaks)
    ]
    y2 <- graphics::hist(
      y2,
      breaks = emd_breaks,
      plot   = FALSE
    )$counts
    
    ## Use iterative algorithm to compute EMD
    emdist::emd2d(
      A        = matrix(y1),
      B        = matrix(y2),
      max.iter = 5000
    )
  }
  
  ## Initialise EMDs for all pairs of batches
  bb <- sort(unique(annotation$Batch))
  nb <- length(bb)
  emd_prenorm  <- matrix(data = 0, nrow = nb, ncol = nb)
  emd_postnorm <- matrix(data = 0, nrow = nb, ncol = nb)
  rownames(emd_prenorm)    <-
    colnames(emd_prenorm)  <-
    rownames(emd_postnorm) <-
    colnames(emd_postnorm) <- bb
  
  niter <- floor( # number of computations required
    ((nrow(expand_grid('a' = 1:nb, 'b' = 1:nb)))-nb) / 2
  )
  
  pb <- utils::txtProgressBar(min = 0, max = niter, style = 3) # progress bar
  counter <- 0
  
  ## Iterate over pairs of batches
  for (idx in 1:(nb-1)) {
    
    for (jdx in (idx+1):nb) {
      
      i <- bb[idx]
      j <- bb[jdx]
      iref <- i %in% batch_ref
      jref <- j %in% batch_ref
      if (iref && jref) { # both i and j are reference
        
        ## Identical EMD for pre- and post-normalisation
        emd_prenorm[idx, jdx] <-emd_prenorm[jdx, idx] <-
          emd_postnorm[idx, jdx] <- emd_postnorm[jdx, idx] <-
          emd_between_batches(
            batch1 = i, type1 = 'Reference',
            batch2 = j, type2 = 'Reference'
          )
      } else if (iref) { # i is reference, j is target
        
        ## Pre-normalisation EMD
        emd_prenorm[idx, jdx] <- emd_prenorm[jdx, idx] <-
          emd_between_batches(
            batch1 = i, type1 = 'Reference',
            batch2 = j, type2 = 'Pre-norm'
          )
        if (include_postnorm) {
          
          ## Post-normalisation EMD
          emd_postnorm[idx, jdx] <- emd_postnorm[jdx, idx] <-
            emd_between_batches(
              batch1 = i, type1 = 'Reference',
              batch2 = j, type2 = 'Post-norm'
            )
        }
      } else if (jref) { # i is target, j is reference
        
        ## Pre-normalisation EMD
        emd_prenorm[idx, jdx] <- emd_prenorm[jdx, idx] <-
          emd_between_batches(
            batch1 = i, type1 = 'Pre-norm',
            batch2 = j, type2 = 'Reference'
          )
        if (include_postnorm) {
          
          ## Post-normalisation EMD
          emd_postnorm[idx, jdx] <- emd_postnorm[jdx, idx] <-
            emd_between_batches(
              batch1 = i, type1 = 'Post-norm',
              batch2 = j, type2 = 'Reference'
            )
        }
      } else { # both i and j are target
        
        ## Pre-normalisation EMD
        emd_prenorm[idx, jdx] <- emd_prenorm[jdx, idx] <-
          emd_between_batches(
            batch1 = i, type1 = 'Pre-norm',
            batch2 = j, type2 = 'Pre-norm'
          )
        if (include_postnorm) {
          
          ## Post-normalisation EMD
          emd_postnorm[idx, jdx] <- emd_postnorm[jdx, idx] <-
            emd_between_batches(
              batch1 = i, type1 = 'Post-norm',
              batch2 = j, type2 = 'Post-norm'
            )
        }
      }
      counter <- counter + 1
      utils::setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  
  list(
    'prenorm'  = emd_prenorm,
    'postnorm' =
      if (include_postnorm) { emd_postnorm } else { NULL }
  )
}

## Function: compress EMD heatmap ----

compress_heatmap <- function(
    x,              # EMD values from `cross_batch_emd`
    idcs,           # names of rows/columns to compress
    name_to = 'ref' # name for compressed rows/columns
) {
  
  ## Resolve heatmap rows/columns to compress or keep
  idcs_compress <- rownames(x)[rownames(x) %in% idcs]
  idcs_keep     <- which(!rownames(x) %in% idcs)
  
  ## Get average values of compressed rows/colums
  avg <- rowMeans(x[idcs_keep, idcs_compress])
  
  ## Create a compressed representation of the heatmap
  res <- rbind(
    cbind(
      x[idcs_keep, idcs_keep, drop = FALSE], avg
    ),
    c(avg, 0)
  )
  rownames(res)[nrow(res)]   <-
    colnames(res)[ncol(res)] <- name_to
  
  res
}

## Function: create training data per batch ----

create_downsampled_batch_aggregates <- function(
    batches,                         # batches to include
    annotation,                      # sample-level annotation
    cohorts,                         # cohorts to include
    fpath_fcs_comptrans,             # path to pre-processed FCS files
    fpath_training_files = NULL,     # path to dedicated CytoNorm training files
    channels             = channels, # channel-marker pairs
    cells_per_agg        = 1e6,      # training aggregate size
    seed                 = 1,        # seed for reproducibility
    only_return_fnames   = FALSE     # whether to only return file names, not...
                                     # ...generate them
) {
  
  ## Initialise file names per batch
  fnames <- vector(mode = 'list', length = length(batches))
  if (is.null(fpath_training_files)) { # train on outputs of pre-processing
    counter <- 1
    
    ## Locate sample files per batch
    for (idx_batch in batches) {
      
      fnames_batch <- file.path(
        fpath_fcs_comptrans,
        annotation$FileName[
          annotation$Cohort %in% cohorts & annotation$Batch == idx_batch
        ]
      )
      fnames_batch <- fnames_batch[file.exists(fnames_batch)]
      fnames[[counter]] <- fnames_batch
      counter <- counter + 1
    }
  } else { # train on dedicated, independently pre-processed samples
    
    counter <- 1
    
    ## Locate sample files per batch (string-matching via numeric part of name)
    for (idx_batch in batches) {
      
      fnames_train <- dir(fpath_training_files)
      fnames_ref   <- file.path(
        fpath_fcs_comptrans,
        annotation$FileName[
          annotation$Cohort %in% cohorts & annotation$Batch == idx_batch
        ]
      )
      fnames_ref        <- fnames_ref[file.exists(fnames_ref)]
      fnames_ref        <- basename(fnames_ref)
      matching          <- sample_match(ref = fnames_ref, real = fnames_train)
      fnames_train      <- file.path(fpath_training_files, matching$Real)
      fnames_train      <- fnames_train[basename(fnames_train) != 'NA']
      fnames[[counter]] <- fnames_train
      
      counter <- counter + 1
    }
  }
  
  ## Determine training data filenames
  str_batch <- formatC(batches, width = 2, format = 'd', flag = '0')
  fnames_norm_agg <- file.path(
    fpath_agg_normtrain_batchwise,
    paste0('Aggregate_Batch', str_batch, '.fcs')
  )
  if (only_return_fnames) { # only return the names
    return(fnames_norm_agg)
  }
  
  ## Create directory for training data files
  if (!file.exists(fpath_agg_normtrain_batchwise)) {
    dir.create(fpath_agg_normtrain_batchwise)
  }
  
  message('Creating aggregates training files per batch')
  
  counter <- 1
  
  ## Iterate over batch
  for (idx_batch in batches) {
    
    message('\tbatch #', idx_batch, ' (', counter, '/', length(batches), ')')
    f <- fnames[[counter]]
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    ## Generate and save data aggregate
    agg <- FlowSOM::AggregateFlowFrames(
      fileNames = f,
      cTotal    = cells_per_agg,
      channels  = names(channels)[
        c(idcs_channels_lineage, idcs_channels_state)
      ],
      silent    = TRUE
    )
    fname_norm_agg <- fnames_norm_agg[counter]
    flowCore::write.FCS(agg, fname_norm_agg)
    
    counter <- counter + 1
  }
  
  fnames_norm_agg
}

## Function: train a normalisation model ----

train_norm_model <- function(
    batches,               # batches to include
    batch_ref,             # reference batch(es)
    batch_norm,            # target batch(es)
    batch_ref_flag,        # flag to use for reference (not in `batches` nor 1)
    annotation,            # sample-level annotation
    cohorts,               # cohorts to include
    channels,              # channel-marker pairs
    idcs_channels_lineage, # indices of lineage markers in `channels`
    idcs_channels_state,   # indices of state markers in `channels`
    fname_norm_model,      # full path and name to save the trained model
    fpath_norm_model,      # directory to save trained model and intermediates
    fpath_training_files = NULL,  # path to dedicated CytoNorm training files
    train_on_qc          = FALSE, # whether to use QC files for training
    fnames_qc            = NULL,  # QC files per batch to use for training
    cells_per_agg        = 1e6,   # cells per batch training aggregate
    per_batch            = TRUE,  # whether to train a separate normalisation...
                                  # ...per target batch
    quantile_values      = ((1:n_quantiles)/n_quantiles)[-n_quantiles],
    # ^ quantile values for fitting CytoNorm splines
    quantile_limits      = c(-50, 260), # quantile limits for anchoring splines
    norm_xdim            = 12,          # SOM width for CytoNom FlowSOM model
    norm_ydim            = 12,          # SOM height for CytoNom FlowSOM model
    norm_n_metaclusters  = 5,     # metacluster count for CytoNom FlowSOM model
    train_on_lin         = TRUE,  # whether to use lineage marker to train...
                                  # ...CytoNorm FlowSOM model
    train_on_state       = FALSE, # whether to use state markers to train...
                                  # ...CytoNorm FlowSOM model
    transform_lin        = TRUE,# whether to train lineage marker normalisations
    transform_state      = TRUE,  # whether to train state marker normalisations
    training_data_exists = FALSE, # whether batch training aggregates exist
    seed                 = 1      # random seed for reproducibility
) {
  
  ## Resolve batch flags for CytoNorm
  batch_flags                             <- batches
  batch_flags[batch_flags %in% batch_ref] <- batch_ref_flag
  
  ## Resolve separate/joint normalisation
  if (!per_batch) {
    batch_flags[batch_flags %in% batch_norm] <- '1'
  }
  
  ## Resolve channels to use for clustering and to train normalisation for
  training_channels <- c()
  if (train_on_lin) {
    training_channels <- names(channels)[idcs_channels_lineage]
  }
  if (train_on_state) {
    training_channels <- c(
      training_channels,
      names(channels)[idcs_channels_state]
    )
  }
  target_channels <- c()
  if (transform_lin) {
    target_channels <- names(channels)[idcs_channels_lineage]
  }
  if (transform_state) {
    target_channels <- c(
      target_channels,
      names(channels)[idcs_channels_state]
    )
  }
  
  if (train_on_qc) { # QC files to be used for training
    
    if (is.null(fpath_training_files)) { # use pre-processed files
      
      fnames_training <- fnames_qc
    } else { # use dedicated training files
      
      fnames_training <- dir(fpath_training_files)
      matching <- sample_match(ref = fnames_qc, real = fnames_training)
      fnames_training <- file.path(fpath_training_files, matching$Real)
    }
  } else { # batch aggregates to be used for training
    
    ## Generate training data (or get file names if already generated)
    fnames_training <- create_downsampled_batch_aggregates(
      batches,
      annotation,
      cohorts,
      fpath_fcs_comptrans,
      fpath_training_files, # (can be NULL)
      channels,
      cells_per_agg,
      seed,
      only_return_fnames = training_data_exists
    )
  }
  
  ## Train normalisation model
  norm_model <- CytoNorm::CytoNorm.train(
    files          = fnames_training,
    labels         = batch_flags,
    channels       = target_channels,
    transformList  = NULL,
    outputDir      = fpath_norm_model,
    FlowSOM.params = list(
      nCells   = 1e6,
      xdim     = norm_xdim,
      ydim     = norm_ydim,
      nClus    = norm_n_metaclusters,
      channels = training_channels
    ),
    normMethod.train = CytoNorm::QuantileNorm.train,
    normParams       = list(
        quantileValues = quantile_values,
        goal           = as.character(batch_ref_flag),
        limit          = quantile_limits
    ),
    seed      = seed,
    clean     = TRUE,
    plot      = FALSE,
    verbose   = FALSE,
    recompute = TRUE
  )
  
  ## Save the trained model
  saveRDS(norm_model, fname_norm_model)
  
  norm_model
}

## Function: apply a trained normalisation model ----

apply_norm_model <- function(
    norm_model,           # trained normalisation model
    batches,              # batches to include
    annotation,           # sample-level annotation
    cohorts,              # cohorts to include
    fpath_fcs_comptrans,  # path to pre-processed FCS files
    fpath_output,         # path to save normalised FCS files
    batch_ref,            # reference batch(es)
    batch_norm,           # target batch(es)
    fpath_norm_tmp = NULL # path to save temporary files
) {
  
  ## Define function to locate target files and get batch of origin for each
  # (This is a stub to make normalisation of a separate, corresponding set of
  # files possible)
  get_target_files <- function(fpath) {
    
    if (is.null(fpath)) {
      fpath <- fpath_fcs_comptrans
    }
    
    fnames_target  <- dir(fpath) # candidates
    fnames_ref     <- file.path( # corresponding analysis files
      fpath_fcs_comptrans,
      annotation$FileName[
        annotation$Cohort %in% cohorts &
        annotation$Batch %in% batches
      ]
    )
    fnames_ref     <- fnames_ref[file.exists(fnames_ref)]
    fnames_ref     <- basename(fnames_ref)
    
    ## Get batches per target file
    batches_target <- annotation$Batch[
      match(fnames_ref, annotation$FileName)
    ]
    
    ## Match target files to the corresponding pre-processed files
    matching       <- sample_match(ref = fnames_ref, real = fnames_target)
    fnames_target  <- file.path(fpath, matching$Real)
    idcs_keep      <- which(basename(fnames_target) != 'NA')
    fnames_target  <- fnames_target[idcs_keep]
    batches_target <- batches_target[idcs_keep]
    n_files        <- length(fnames_target)
    attributes(fnames_target)$batches_target <- batches_target
    fnames_target
  }
    
  ## Resolve target files to normalise and corresponding batches
  fnames_target  <- get_target_files(NULL)
  batches_target <- attributes(fnames_target)$batches_target
  
  ## Iterate over target batches (to normalise against reference)
  for (idx in seq_along(batch_norm)) {
    
    batch <- batch_norm[idx]
    message('Normalising samples in batch ', batch)
    fnames <- fnames_target[as.integer(batches_target) == batch]
    
    ## Apply previously trained CytoNorm model
    CytoNorm::CytoNorm.normalize(
      model                 = norm_model,
      files                 = fnames,
      labels                =
        as.character(rep(batch, times = length(fnames))),
      transformList         = NULL,
      transformList.reverse = NULL,
      outputDir             = fpath_output,
      prefix                = '',
      normMethod.normalize  = CytoNorm::QuantileNorm.normalize,
      clean                 = TRUE,
      verbose               = FALSE
    )
  }
  
  ## Copy reference-batch files (no normalisation to be applied)
  fnames_refbatch <- fnames_target[batches_target %in% batch_ref]
  
  ## (Copy function called with up to 100 files at a time)
  n <- length(fnames_refbatch)
  s <- seq( 
    from = 1, to = n, by = 100
  ) 
  for (i_start in s) {
    
    i_end <- min(i_start + 99, n)
    message(
      'Copying reference-batch files ', i_start, '-', i_end, ' out of ', n
    )
    fnames <- fnames_refbatch[i_start:i_end]
    file.copy(
      from = fnames,
      to = file.path(fpath_output, basename(fnames))
    )
  }
}

