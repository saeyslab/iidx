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

## iidx internal module: 99_CheckInputs.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Resolves inputs FCS data.


## Function: check if a variable is specified ----

specified <- function(x) {
  !is.null(x) && length(x) > 0
}

## Function: gather paths to all valid FCS inputs ----

resolve_input_fcs <- function(
    annotation,          # sample-level annotation
    fpath_fcs,           # path to raw FCS files directory
    fpath_fcs_comptrans, # path to compensated & transformed FCS files directory
    fpath_fcs_norm,      # path to normalised FCS files directory
    normalise,           # whether normalisation is enabled
    fpath_fcs_unstained = NULL, # path to unstained FCS files directory
    cohorts             = NULL, # cohorts to include
    batch_remove        = NULL, # batches to exclude
    verbose             = TRUE  # print messages about exclusion of files
) {
  
  ## Match annotation to existing files
  in_annotation <- annotation$FileName
  in_path       <- dir(fpath_fcs)
  fnames        <- file.path(
    fpath_fcs,
    in_annotation[in_annotation %in% in_path]
  )
  
  ## Report missing files
  n_unmatched <- length(in_annotation)-length(fnames)
  if (n_unmatched > 0) {
    
    if (verbose) {
      message(
        n_unmatched, '/', nrow(annotation),
        ' samples excluded because of missing FCS files'
      )
    }
  }
  
  ## Exclude missing files from annotation and paths
  a <- annotation[
    match(basename(fnames), annotation$FileName), , drop = FALSE 
  ]
  fnames_fcs_raw <- file.path(fpath_fcs, a$FileName)
  
  ## Filter samples by cohorts
  if (specified(cohorts)) {
    
    ## Flag and count samples from excluded cohorts
    idcs_wrong_cohort <- which(!a$Cohort %in% cohorts)
    n_wrong_cohort    <- length(idcs_wrong_cohort)
    
    ## Report files excluded by cohort filtering
    if (n_wrong_cohort > 0) {
      
      if (verbose) {
        message(
          n_wrong_cohort, '/', length(fnames_fcs_raw),
          ' samples excluded by cohort filtering'
        )
      }
      
      ## Exclude flagged files from annotation and paths
      fnames_fcs_raw <- fnames_fcs_raw[-idcs_wrong_cohort]
      a <- a[
        -idcs_wrong_cohort, , drop = FALSE
      ]
    }
  }
  
  ## Exclude samples by batch
  if (specified(batch_remove)) {
    
    ## Flag and count samples from excluded batches
    idcs_wrong_batch <- which(a$Batch %in% batch_remove)
    n_wrong_batch    <- length(idcs_wrong_batch)
    
    ## Report files excluded by batch filtering
    if (n_wrong_batch > 0) {
      
      if (verbose) {
        message(
          n_wrong_batch, '/', length(fnames_fcs_raw),
          ' samples excluded by batch filtering'
        )
      }
      
      ## Exclude flagged files from annotation and paths
      fnames_fcs_raw <- fnames_fcs_raw[-idcs_wrong_batch]
      a <- a[
        -idcs_wrong_batch, , drop = FALSE
      ]
    }
  }
  
  ## Check for missing compensated and transformed files
  fnames_fcs_comptrans <- NULL
  if (file.exists(fpath_fcs_comptrans)) { # compensation and transformation done
    
    # Report missing compensated and transformed files
    fnames_fcs_comptrans <- file.path( 
      fpath_fcs_comptrans,
      basename(fnames_fcs_raw)
    )
    idcs_missing <- which(!file.exists(fnames_fcs_comptrans))
    n_missing    <- length(idcs_missing)
    if (n_missing > 0) {
      
      if (verbose) {
        message(
          n_missing, ifelse(n_missing == 1, ' file', ' files'),
          ' missing from compensated & transformed'
        )
      }
    }
  }
  
  ## Check for missing normalised files
  fnames_fcs_norm <- NULL
  if (file.exists(fpath_fcs_norm)) { # normalisation done
    
    # Report missing normalised files
    fnames_fcs_norm <- file.path(
      fpath_fcs_norm,
      basename(fnames_fcs_comptrans)
    )
    idcs_missing <- which(!file.exists(fnames_fcs_comptrans))
    n_missing    <- length(idcs_missing)
    if (n_missing > 0) {
      
      if (verbose) {
        message(
          n_missing, ifelse(n_missing == 1, ' file', ' files'),
          ' missing from normalised'
        )
      }
    }
  }
  
  ## Locate unstained sample files
  if (!is.null(fpath_fcs_unstained)) { # unstained samples available
    
    ## Resolve raw unstained file names
    b       <- unique(a$Batch)
    batches <- b[!b%in%batch_remove]
    f       <- file.path(fpath_fcs_unstained, paste0('Batch', batches, '.fcs'))
    fe      <- file.exists(f)
    batches <- batches[fe]
    
    ## Resolve compensated and transformed unstained file names
    fnames_fcs_unstained_raw <- f[fe]
    f <- file.path(fpath_fcs_unstained_comptrans, basename(fnames_fcs_unstained_raw))
    f <- f[file.exists(f)]
    fnames_fcs_unstained_comptrans <- if (length(f)==0) { NULL } else { f }
  }
  
  list(
    'fnames_fcs_raw'                 = fnames_fcs_raw,
    'fnames_fcs_comptrans'           = fnames_fcs_comptrans,
    'fnames_fcs_norm'                = fnames_fcs_norm,
    'fnames_fcs_unstained_raw'       = fnames_fcs_unstained_raw,
    'fnames_fcs_unstained_comptrans' = fnames_fcs_unstained_comptrans,
    'annotation'                     = a
  )
}
