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

## iidx user-level module: 03_OutlierAndNoiseDetection.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Uses FlowSOM-derived features to flag outlier samples and 
# computes background signal levels per marker.


message(
  '## iidx: 03_OutlierAndNoiseDetection.R STARTED ',
  as.character(round(Sys.time()))
)

## Get parameters and auxiliary functions ----

input_files_needed <- FALSE
source('99_SetUp.R')

## Create results directory ----

if (!file.exists(fpath_res03)) {
  dir.create(fpath_res03)
}

## Flag outlier samples ----

if (is.null(dev_type)||n_dev==Inf) { # no outlier detection
  
  message('Skipping outlier detection step since it is disabled')
  
  ## Leave vector of sample-level outlier names empty
  outliers <- c()
  saveRDS(outliers, fname_sample_outliers)
} else {
  
  message('Flagging outlier samples for incompatibility with clustering model')
  
  if (dev_type=='mad') { # median absolute deviation from median
    avg <- median
    dev <- mad
  } else if (dev_type=='sd') { # standard deviation from mean
    avg <- mean
    dev <- sd
  } else {
    stop('Invalid `dev_type`')
  }
  
  ## Load metacluster sizes as percentages of each sample
  features <- readRDS(fname_features_postnorm)
  mc_perc  <- features$metacluster_percentages
  
  ## Plot cut-offs for sample outliers under chosen input parameters
  plot_dim <- ceiling(sqrt(n_metaclusters)) # dimension for square grid of plots
  tryCatch(expr = {dev.off()}, error = function(e) {}) # ensure no open device
  pdf(
    file   = fname_viz_sample_outliers,
    width  = plot_dim*2,
    height = plot_dim*2
  )
  
  ## Initialise plots of sorted percentage distributions
  par(
    mfrow = c(plot_dim, plot_dim),
    mar   = c(2, 2, 2, 2)
  )
  
  ## Iterate over metaclusters
  for (mc in colnames(mc_perc)) {
    
    y <- mc_perc[, mc]
    
    ## Plot distribution
    plot(
      sort(y),
      pch  = 20,
      cex  = 0.3,
      main = mc
    )
    
    ## Draw outlier cut-off
    mu <- avg(y)
    s  <- n_dev*dev(y)
    omin <- mu-s
    omax <- mu+s
    abline(h = omin, lwd = 1.2, col = 'darkred')
    abline(h = omax, lwd = 1.2, col = 'darkred')
  }
  
  ## Initialise plots of unsorted percentage distributions
  par(
    mfrow = c(plot_dim, plot_dim),
    mar   = c(2, 2, 2, 2)
  )
  
  ## Iterate over metaclusters
  for (mc in colnames(mc_perc)) {
    
    y <- mc_perc[, mc]
    
    ## Plot distribution
    plot(
      y,
      pch  = 20,
      cex  = 0.3,
      main = mc
    )
    
    ## Draw outlier cut-off
    mu <- avg(y)
    s  <- n_dev*dev(y)
    omin <- mu-s
    omax <- mu+s
    abline(h = omin, lwd = 1.2, col = 'darkred')
    abline(h = omax, lwd = 1.2, col = 'darkred')
  }
  dev.off()
  
  ## Compute outliers
  sample_outliers_per_mc <-
    lapply(
      seq_len(ncol(mc_perc)), # iterate over metaclusters
      function(idx) {
        
        y    <- mc_perc[, idx]
        mu   <- avg(y)
        s    <- n_dev*dev(y)
        omin <- mu-s
        omax <- mu+s
        
        ## Get indices of extreme values
        idcs <- as.vector(which(y < omin | y > omax))
        
        ## Return corresponding sample names
        basename(rownames(mc_perc)[idcs])
      }
    )
  
  ## Flag samples that are outliers for at least 1 metacluster
  outliers <- sort(unlist(sample_outliers_per_mc))
  saveRDS(outliers, fname_sample_outliers)
}

## Determine signal robustness ----

if (unstained_samples_available) { # noise levels were computed
  
  message('Determining signal robustness based on background signal data')
  
  ## Load noise stats per state marker and batch
  noise <- readRDS(fname_noise)
  noise <- noise[
    noise$Marker%in%as.vector(channels)[idcs_channels_state], , drop = FALSE
  ]
  
  ## Load state marker MFIs of non-outlier samples
  mfis <- readRDS(fname_features_postnorm)$metacluster_MFIs
  mfis <- mfis[!rownames(mfis)%in%outliers, , drop = FALSE]
  
  ## Separate marker names from MFI feature matrix columns
  markers <- as.factor(gsub('^MC[0-9]+[ ]', '', colnames(mfis)))
  
  ## Align sample names with batch numbers
  samples <- rownames(mfis)
  batches <- as.factor(
    annotation$Batch[match(samples, annotation$FileName)]
  )
  
  ## Compute robustness (FALSE/TRUE) per marker of each sample at 5 levels...
  ## ...of noise cut-offs
  n <- nlevels(markers)*nlevels(batches) # number of combinations
  counter <- 0
  
  res <- lapply(
    c(1, 2, 3, 4, 5), # cut-off as number of MADs above background signal median
    function(n_mad) {
      
      ## Set cut-offs for each combination
      noise$Cutoff <- noise$Median+noise$MAD*n_mad
      
      ## Iterate over batches
      res <- lapply(
        levels(batches),
        function(batch) {
          
          i_batch   <- which(batches==batch)
          
          ## Initialise matrix of FALSE/TRUE for robustness status
          batch_res <- matrix(
            NA,
            nrow = sum(batches==batch),
            ncol = length(markers)
          )
          colnames(batch_res) <- markers
          rownames(batch_res) <- rownames(mfis)[i_batch]
          
          ## Get noise cut-offs per combination
          batch_noise   <- noise[noise$Batch==batch, ]
          batch_cutoffs <- batch_noise$Cutoff[
            match(markers, batch_noise$Marker)
          ]
          
          ## Return robustness status per combination
          t(apply(
            X = mfis[i_batch, , drop = FALSE],
            MARGIN = 1,
            FUN = function(x) x>batch_cutoffs
          ))
        }
      )
      res <- do.call(rbind, res)
      res <- res[rownames(mfis), , drop = FALSE]
      res[is.na(res)] <- FALSE
      
      res
    }
  )
  
  ## Save signal robustness results
  robustness <- res
  saveRDS(robustness, fname_robustness)
}

message(
  '## iidx: 03_OutlierAndNoiseDetection.R FINISHED ',
  as.character(round(Sys.time()))
)

