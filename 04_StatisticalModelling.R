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

## iidx user-level module: 04_StatisticalModelling.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Performs a differential expression analysis using previously
# identified cell subsets and sample-level conditions.


message(
  '## iidx: 04_StatisticalModelling.R STARTED ',
  as.character(round(Sys.time()))
)

## Get parameters and auxiliary functions ----

input_files_needed <- FALSE
source('99_SetUp.R')
source('04b_FullExperiments.R')
library(tidyverse)

## Create results directory ----

if (!file.exists(fpath_res04)) {
  dir.create(fpath_res04)
}

## Collect FlowSOM-derived features ----

inputs <- prep_inputs(
  fname_features_postnorm, fname_sample_outliers, annotation, batch_remove
)

## Count number of analysis types to run
n_out <- 2+(length(thresholds)>0) # DA, DS-MFI, (DS-Pheno)

## Run differential analyses and save results ----

message('Performing differential expression analyses')

message('(1/', n_out, ') Differential Abundance tests')
res_fsom_da <- test_da(
  counts      = inputs$counts,
  annotation  = annotation,
  predictors  = predictors,  # potential drivers of DE
  confounders = confounders, # effects to disentangle from predictors
  verbose     = TRUE
)
saveRDS(res_fsom_da, fname_fsom_da)

message('(2/', n_out, ') Differential State (MFI) tests')
res_fsom_ds_mfi <- test_ds(
  annotation  = annotation,
  mfi         = inputs$MFIs,
  counts      = inputs$counts,
  predictors  = predictors,    # potential drivers of DE
  confounders = confounders,   # effects to disentangle from predictors
  state_markers =
    as.vector(channels)[idcs_channels_state],
  parallel    = run_parallel,  # (multi-threading gives a big speed-up here)
  verbose     = TRUE
)
saveRDS(res_fsom_ds_mfi, fname_fsom_ds_mfi)

res_fsom_ds_pheno <- NULL
if (n_out==3) {
  message('(3/3) Differential State (Phenopositivity) tests')
  res_fsom_ds_pheno <- test_ds( 
    annotation  = annotation,
    phenopos    = inputs$PercPos,
    counts      = inputs$counts,
    predictors  = predictors,   # potential drivers of DE
    confounders = confounders,  # effects to disentangle from predictors
    state_markers =
      as.vector(channels)[idcs_channels_state],
    parallel    = run_parallel, # (multi-threading gives a big speed-up here)
    verbose     = TRUE
  )
}
saveRDS(res_fsom_ds_pheno, fname_fsom_ds_pheno)

message(
  '## iidx: 04_StatisticalModelling.R FINISHED ',
  as.character(round(Sys.time()))
)

