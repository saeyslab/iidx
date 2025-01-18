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

## iidx internal module: 00a_CompensationAndTransformation.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines functions for compensation and transformation of large
# amounts of FCS data.


## Function: make flowFrame metadata valid ----
## make_fcs_valid adopted from github.com/stuchly/tviblindi

make_fcs_valid <- function(
  exprs,        # expression matrix
  fcs   = NULL, # flowCore::flowFrame to use as template
  desc1 = NULL  # `desc` slot
) {
  
  if(!is.null(fcs)) {
    
    mn  <- flowCore::markernames(fcs)
    fcs <- flowCore::flowFrame(exprs, parameters = flowCore::parameters(fcs))
    flowCore::markernames(fcs) <- mn
  } else {
    
    mn  <- NULL
    fcs <- flowCore::flowFrame(exprs)
  }
  
  params <- flowCore::parameters(fcs)
  pd     <- NULL
  cols   <- as.vector(pd$name)
  idcs   <- match(cols, pd$name)
  
  if (any(is.na(idcs))) {
    stop('Invalid column specifier')
  }
  
  keyval <- list()
  for (channel_number in seq_len(ncol(exprs))) {
    
    channel_name <- colnames(exprs)[channel_number]
    if (is.null(desc1)) {
      desc1 <- colnames(exprs)[channel_number]
    }
    channel_id    <- paste('$P', channel_number, sep = '')
    channel_range <- max(exprs[, channel_number]) + 1
    channel_min   <- min(0, min(exprs[, channel_number])-1)
    
    plist <- matrix(
      c(channel_name, desc1[channel_number], channel_range, channel_min, channel_range - 1)
    )
    rownames(plist) <- c('name', 'desc', 'range', 'minRange', 'maxRange')
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    
    keyval[[paste('$P', channel_number, 'B', sep = '')]] <- '32'
    keyval[[paste('$P', channel_number, 'R', sep = '')]] <- toString(channel_range)
    keyval[[paste('$P', channel_number, 'E', sep = '')]] <- '0,0'
    keyval[[paste('$P', channel_number, 'N', sep = '')]] <- channel_name
    keyval[[paste('$P', channel_number, 'S', sep = '')]] <- channel_name
  }
  
  params@data <- data.frame(pd)
  if(!is.null(fcs)) {
    
    fcs <- flowCore::flowFrame(exprs, parameters = params)
  } else {
    
    fcs <- flowCore::flowFrame(exprs, parameters = params)
  }
  flowCore::keyword(fcs) <- keyval
  
  if (!is.null(mn)) {
    flowCore::markernames(fcs) <- mn
  }
  
  fcs
}

## Function: compensate and transform FCS file ----

comptrans <- function(
    idx_file,                 # index of file in `fnames_fcs_raw`
    fnames_fcs_raw,           # all full paths to raw FCS files
    fpath_fcs_comptrans,      # directory to save outputs
    cn,                       # channel names
    compensate,               # whether to apply signal compensation
    tf,                       # transformation function
    tf_match_channels,        # whether channel names need to be matched
    tf_listed,                # whether `tf` is a `flowCore::transformList`
    channels,                 # channel-marker pairs
    signal_limits,            # signal range for trimming
    fpath_subset_idcs = NULL, # path to subsets indices (for pre-gating)
    spillover         = NULL  # spillover matrix (if not in the FCS file)
) {
  
  fname <- fnames_fcs_raw[idx_file]
  
  ## Import FCS file
  ff <- flowCore::read.FCS(fname)
  
  ## Check for missing channels in FCS file
  if (any(!cn%in%flowCore::colnames(ff))) {
    
    missing <- cn[!cn%in%flowCore::colnames(ff)]
    stop(
      'Sample "', basename(fname),
      '" lacks channels specified by transformation function:\n\t',
      paste(missing, collapse = '\n\t')
    )
  }
  
  ## Check for missing channels in `channels` input parameter
  if (any(!names(channels)%in%flowCore::colnames(ff))) {
    
    missing <- names(channels)[!names(channels)%in%flowCore::colnames(ff)]
    stop(
      'Sample "', basename(fname),
      '" lacks channels specified by `channels` input parameter:\n\t',
      paste(missing, collapse = '\n\t')
    )
  }
  
  ## Load subsetting indices if pre-gating applied to this sample
  if (!is.null(fpath_subset_idcs)) {
    
    sn <- file.path(fpath_subset_idcs, sub('.fcs', '.RDS', basename(fname)))
    if (file.exists(sn)) {
      idcs <- readRDS(sn)
      ff <- ff[idcs, ]
    }
  }
  
  ## Compensate using spillover matrix from FCS file or manually provided one
  if (compensate) {
    
    if (is.null(spillover)) {
      spillover <- flowCore::keyword(ff)[['$SPILLOVER']]
    }
    ff <- flowCore::compensate(ff, spillover)
  }
  
  ## Transform expression data
  if (is.function(tf)) { # `tf` is simple R function
    
    ff@exprs[, cn] <- tf(ff@exprs[, cn])
  } else if (tf_listed) { # `tf` as `flowCore::transformList`
    
    if (tf_match_channels) {
      # ^ channel names need to be matched (e.g., handle 'Comp-' prefix)
      
      cn_ff <- flowCore::colnames(ff)
      m <- match_channel_names(
        ref  = names(tf@transforms),
        real = cn_ff
      )
      
      ## Change channel names, transform and revert back
      flowCore::colnames(ff)[m$Index] <- m$Reference
      ff                              <- flowCore::transform(ff, tf)
      flowCore::colnames(ff)          <- cn_ff
    } else { # no matching needed
      
      ff <- flowCore::transform(ff, tf)
    }
  } else {
    
    stop('Transformation function in `tf` is not valid')
  }
  
  ## Exclude channels not used in analysis
  ff <- ff[, names(channels)]
  
  ## Make sure marker names are correct
  flowCore::markernames(ff) <- channels[names(flowCore::markernames(ff))]
  
  ### Trim by signal limits
  if (!is.null(signal_limits)) {
    
    signal_min <- signal_limits[1]
    signal_max <- signal_limits[2]
    
    ## Only keep events within signal range
    mask <- rowSums(
      ff@exprs[, cn]<signal_min)==0 &
      rowSums(ff@exprs[, cn]>signal_max
    )==0
    ff <- ff[mask, , drop = FALSE]
  }
  
  ## Add 'sampleid' channel with `idx_file`
  id <- matrix(
    rep(idx_file),
    nrow = flowCore::nrow(ff),
    ncol = 1
  )
  colnames(id) <- 'sampleid'
  ff <- flowCore::fr_append_cols(ff, id)
  
  ## Ensure flowCore::flowFrame metadata is valid
  ff <- make_fcs_valid(ff@exprs, ff)
  
  ## Save result
  flowCore::write.FCS(
    x         = ff,
    filename = file.path(fpath_fcs_comptrans, basename(fname))
  )
}
