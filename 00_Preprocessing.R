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

## iidx user-level module: 00_Preprocessing.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Applies compensation (based on spillover) and transformation to 
# raw sample data and saves the resulting pre-processed FCS files.


message(
  '## iidx: 00_Preprocessing.R STARTED ',
  as.character(round(Sys.time()))
)

## Get parameters and auxiliary functions ----

input_files_needed <- TRUE
source('99_SetUp.R')
source('00a_CompensationAndTransformation.R')

## Create results directory ----

if (!file.exists(fpath_res00)) {
  dir.create(fpath_res00)
}

## Resolve channels to transform ----

ff <- flowCore::read.FCS(fnames_fcs_raw[1]) # representative sample
cn <- flowCore::colnames(ff)
cn <- cn[!cn %in% channels_remove]
if (!is.null(notrans) && length(notrans)>0) {
  cn_notrans <- rowSums(sapply(notrans, function(x) grepl(x, cn)))>0
} else {
  cn_notrans <- rep(FALSE, times = length(cn))
}
cn <- cn[!cn_notrans]

if (any(!cn%in%names(channels))) {
  stop(
    'Some channels lack marker names in the `channels` input parameter:\n\t',
    paste(cn[!cn%in%names(channels)], sep = '\n\t')
  )
}

## Resolve transformation formula ----

tf_listed          <- FALSE # flowCore::transformList provided?
tf_match_channels  <- FALSE
# ^ transform channel names matchable but not identical to flowFrame channels?

if (is.null(tf)) { # transform to be extracted from WSP
  
  if (is.null(fname_wsp)) {
    stop(
      'No transformation (`tf`) and no FlowJo workspace (`fname_wsp`) specified
    ')
  }
  
  message('Extracting transform function from FlowJo workspace: ', fname_wsp)
  
  ## Parse FlowJo workspace
  ws <- CytoML::open_flowjo_xml(fname_wsp)
  gs <- CytoML::flowjo_to_gatingset(ws, path = fpath_fcs_ref, name = 2)
  
  ## Extract transformation as flowCore::transformList
  tf <- flowWorkspace::gh_get_transformations(gs[[1]])
  tf <- flowCore::transformList(names(tf), tf)
  
  ## Match channel names in list to those from an FCS file
  m <- match_channel_names(ref = names(tf@transforms), real = cn)
  if (any(is.na(m$Real))) {
    unmatched <- m$Reference[is.na(m$Real)]
    stop(
      'Some channels listed in FlowJo transform cannot be matched:\n\t',
      paste(unmatched, collapse = '\n\t')
    )
  }
  tf_listed <- TRUE
  tf_match_channels <- TRUE
} else if (is.function(tf)) { # transform as simple R function
  
  message('Checking provided channel-independent transformation function')
  
  ## Check if function is valid
  e <- tryCatch(expr = {
      res <- tf(as.matrix(iris[, 1:4])-3.)
    },
    error = function(e) e
  )
  if ('error'%in%class(e)) {
    stop(
      'Provided `tf` function failed test with numeric matrix as input'
    )
  }
} else if ('transformList'%in%class(tf)) { # transform per channel via flowCore
  
  message('Checking provided flowCore::transformList')
  
  ## Match channel names in list to those from an FCS file
  m <- match_channel_names(ref = names(tf@transforms), real = cn)
  if (any(is.na(m$Real))) {
    unmatched <- m$Reference[is.na(m$Real)]
    stop(
      'Some channels listed in flowCore::transformList cannot be matched:\n\t',
      paste(unmatched, collapse = '\n\t')
    )
  }
  tf_listed <- TRUE
} else if (is.character(tf) && length(tf)==1) { # transform named
  
  message('Resolving requested transformation function: "', tf, '"')
  
  if (tf=='logicle') {
    
    if (
      is.null(tf_cofactor) ||
      !is.numeric(tf_cofactor) ||
      length(tf_cofactor)!=1
    ) {
      stop('Invalid `tf_cofactor` value for logicle transform')
    }
    tf <- flowCore::transformList(
      from = cn,
      tfun = flowCore::logicleTransform(w = tf_cofactor)
    )
    tf_listed <- TRUE
  } else if (tf=='arcsinh') {
    
    if (
      is.null(tf_cofactor) ||
      !is.numeric(tf_cofactor) ||
      length(tf_cofactor)!=1
    ) {
      stop('Invalid `tf_cofactor` value for arcsinh transform')
    }
    tf <- flowCore::transformList(
      from = cn,
      tfun = flowCore::arcsinhTransform(b = tf_cofactor)
    )
    tf_listed <- TRUE
  } else {
    
    stop('Requested transformation "', tf, '" is not implemented')
  }
} else {
  
  stop('Transformation function (`tf`) misconfigured')
}

## Compensate & transform stained samples ----

if (compensate) {
  
  message('Generating compensated and transformed FCS files')
} else {
  
  message('Generating transformed FCS files (no compensation applied)')
}

## Create directory for compensated and transformed FCS files
if (!file.exists(fpath_fcs_comptrans)) {
  dir.create(fpath_fcs_comptrans)
}

n_fcs <- length(fnames_fcs_raw) # `fnames_fcs_raw` generated by `99_Aux.R`

if (run_parallel) { # multi-threading enabled
  
  ## Set up parallel processing
  cores <- detectCores()
  cl    <- makeCluster(cores[1]-1)
  registerDoSNOW(cl)
  
  pb <- utils::txtProgressBar(min = 0, max = n_fcs, style = 3) # progress bar
  
  opts <- list(
    'progress' = function(i) utils::setTxtProgressBar(pb, i)
  )
  
  ## Process all FCS files
  tmp <- foreach(
    i             = seq_along(fnames_fcs_raw),
    .combine      = c,
    .options.snow = opts
  ) %dopar% {
    
    comptrans(
      idx_file            = i,
      fnames_fcs_raw      = fnames_fcs_raw,
      fpath_fcs_comptrans = fpath_fcs_comptrans,
      cn                  = cn,
      compensate          = compensate,
      tf                  = tf,
      tf_match_channels   = tf_match_channels,
      tf_listed           = tf_listed,
      channels            = channels,
      signal_limits       = signal_limits,
      fpath_subset_idcs   = fpath_subset_idcs,
      spillover           = NULL
    )
    
    NULL
  }
  close(pb)
  stopCluster(cl)
} else {
  
  pb <- utils::txtProgressBar(min = 0, max = n_fcs, style = 3) # progress bar
  
  ## Process all FCS files
  for (i in seq_along(fnames_fcs_raw)) {
    
    comptrans(
      idx_file            = i,
      fnames_fcs_raw      = fnames_fcs_raw,
      fpath_fcs_comptrans = fpath_fcs_comptrans,
      cn                  = cn,
      compensate          = compensate,
      tf                  = tf,
      tf_match_channels   = tf_match_channels,
      tf_listed           = tf_listed,
      channels            = channels,
      signal_limits       = signal_limits,
      fpath_subset_idcs   = fpath_subset_idcs,
      spillover           = NULL
    )
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
}

## Compensate & transform unstained samples ----

if (unstained_samples_available) {
  
  ## Check that all files per batch exist
  b       <- unique(annotation$Batch)
  batches <- b[!b%in%batch_remove]
  f       <- file.path(fpath_fcs_unstained, paste0('Batch', batches, '.fcs'))
  fe      <- file.exists(f)
  if (!all(fe)) {
    
    stop(
      'Unstained samples for the following batches do not exist: ',
      paste(batches[!fe], collapse = ', ')
    )
  }
  fnames_fcs_unstained <- f
  
  ## Create directory for compensated and transformed unstained FCS files
  if (!file.exists(fpath_fcs_unstained_comptrans)) {
    dir.create(fpath_fcs_unstained_comptrans)
  }
  
  ## Process all unstained FCS files
  for (i in seq_along(batches)) {
    batch <- batches[i]
    
    message(
      'Pre-processing batch ', batch,
      ' unstained sample (', i, '/', length(batches), ')'
    )
    
    fname <- fnames_fcs_unstained[i]
    
    ## Extract spillover matrix if compensation needed
    if (compensate) {
      
      ff <- flowCore::read.FCS(fname)
      spillover <- flowCore::keyword(ff)[['$SPILLOVER']]
      
      ## If not from the file itself, from stained sample from same batch
      if (is.null(spillover)) {
        
        fnames_ff_batch <- file.path(
          fpath_fcs, annotation$FileName[annotation$Batch==batch]
        )
        fname_ff_spill  <- fnames_ff_batch[file.exists(fnames_ff_batch)][1]
        ff_spill        <- flowCore::read.FCS(fname_ff_spill)
        spillover       <- flowCore::keyword(ff_spill)[['$SPILLOVER']]
        
        stopifnot(!is.null(spillover))
      }
    }
    
    comptrans(
      idx_file            = 1,
      fnames_fcs_raw      = fname,
      fpath_fcs_comptrans = fpath_fcs_unstained_comptrans,
      cn                  = cn,
      compensate          = compensate,
      tf                  = tf,
      tf_match_channels   = tf_match_channels,
      tf_listed           = tf_listed,
      channels            = channels,
      signal_limits       = signal_limits,
      fpath_subset_idcs   = NULL,
      spillover           = spillover
    )
  }
  
  ## Aggregate unstained samples ----
  
  message('Constructing aggregate of data from unstained samples')
  
  set.seed(1)
  agg <- FlowSOM::AggregateFlowFrames(
    fileNames = dir(
      fpath_fcs_unstained_comptrans,
      pattern = '[.]fcs',
      full.names = TRUE
    ),
    cTotal   = 6e6,
    channels = names(channels)
  )
  flowCore::write.FCS(
    x        = agg,
    filename = fname_agg_unstained
  )
}

message(
  '## iidx: 00_Preprocessing.R FINISHED ',
  as.character(round(Sys.time()))
)
