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

## iidx user-level module: 99_SetUp.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)

# This script specifies inputs to the analysis pipeline and values of tunable
# hyperparameters.

library(tidyverse)

## Multi-threading ----
run_parallel <- TRUE
# ^ whether to allow parallelising some for loops

## Input FCS files ----
# fpath_fcs <- '/auto/net/fs4-Flowcytometry/Storage/u_ysa/Projects/twin_study/52e4c244-0e9b-4573-8e11-d10c249a871c-twinstudy2/TNK/fcs files_comped and pregated'
fpath_fcs <- '../ts/TwinStudy/TwinStudyInputFCSFiles/'

## Data annotation ----
annotation <- read.csv('Annotation.csv')
# ^ required columns:
## -> FileName (string)
## -> Batch (integer)
## -> QC (Boolean)
## -> Cohort (string)
## -> phenotype conditions (predictors to use in differential analyses)

## Condition specs ----
## -> phenotype conditions to use as predictors in statistical models
predictors <- c('Sex', 'Age')
## -> phenotype conditions to use as potential confounders in statistical models
confounders <- c('Sex', 'Age')

## FlowJo workspace specification ----

fname_wsp <- NULL
# ^ name of FlowJo workspace to extract transform from, or NULL
fpath_fcs_ref <- NULL
# ^ directory w/ single FCS file in the workspace (for flowWorkspace to work)

## Panel specs ----

### Channel-marker pairs ----
channels <- c(
  'FSC-A'  = 'FSC-A',    #  1
  'FSC-H'  = 'FSC-H',    #  2
  'SSC-A'  = 'SSC-A',    #  3
  'B515-A' = 'RSV Cav1', #  4
  'B610-A' = 'CD141',    #  5
  'B660-A' = 'CD123',    #  6
  'B710-A' = 'CD16',     #  7
  'B780-A' = 'IgD',      #  8
  'G575-A' = 'CD32',     #  9
  'G610-A' = 'CD40',     # 10
  'G660-A' = 'CD85j',    # 11
  'G710-A' = 'CD11c',    # 12
  'G780-A' = 'CXCR3',    # 13
  'R670-A' = 'IgA',      # 14
  'R730-A' = 'CD27',     # 15
  'R780-A' = 'CD19',     # 16
  'U390-A' = 'CD1c',     # 17
  'U450-A' = 'Viability',# 18
  'U500-A' = 'CD21',     # 19
  'U570-A' = 'TACI',     # 20
  'U660-A' = 'HLA-DR',   # 21
  'U740-A' = 'IgG',      # 22
  'U785-A' = 'CD20',     # 23
  'V450-A' = 'IL-21R',   # 24
  'V510-A' = 'CD14',     # 25
  'V570-A' = 'IgM',      # 26
  'V605-A' = 'BAFF-R',   # 27
  'V655-A' = 'CD10',     # 28
  'V710-A' = 'CD23',     # 29
  'V750-A' = 'CXCR5',    # 30
  'V785-A' = 'CD64',     # 31
  'Time'   = 'Time'      # 32
)


### Ad hoc panel corrections ----
channels_remove <- c('remove_from_FS_FM', 'FJ Sample QC')
channels_rename <- c( # from - to
  'U395-A' = 'U390-A', 
  'U670-A' = 'U660-A'
)

### Lineage vs. state marker division ----
# idcs_channels_lineage <- c(4, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 24, 25, 27, 28, 29, 30)
# idcs_channels_state   <- c(5, 6, 7, 8, 9, 12, 15, 21, 26, 31)
idcs_channels_lineage <- c(8, 14, 15, 16, 19, 21, 22, 23, 26, 28)
idcs_channels_state   <- c(9, 10, 11, 12, 13, 17, 20, 24, 27, 29, 30)

### Phenopositivity thresholds ----
thresholds <- c(
  'IgD' = 1.70,
  'CD1c' = 1.85,
  'CD21' = 1.40,
  'HLA-DR' = 2.00,
  'IgG' = 1.70,
  'IL-21R' = 1.50,
  'TACI' = 1.50,
  #'CD14' = 1.60,
  'IgM' = 1.50,
  'BAFF-R' = 2.30,
  'CD10' = 1.80,
  'CD23' = 1.60,
  #'RSV Cav1' = 1.75,
  #'CD141' = 1.85,
  #'CD123' = 1.70,
  #'CD16' = 1.70,
  'CD32' = 1.90,
  'CD40' = 2.30,
  'CD85j' = 1.55,
  'CD11c' = 1.75,
  'CXCR3' = 2.00,
  'CXCR5' = 1.60,
  'IgA' = 2.20,
  'CD27' = 1.75,
  'CD19' = 1.15
)

## Pre-processing specs ----

### Indices of cell subsets ----
fpath_subset_idcs <- NULL
# ^ optional path to cell (row) indices per file to use exclusively

### Compensation from spillover ----
compensate <- TRUE

### Transformation ----
tf <- readRDS('TransformList.RDS')
# ^ options:
# NULL:                    transformation is taken from WSP file
# flowCore::transformList: a flowCore transformList that can be used directly
# function:                function to use directly on numeric matrix (for requested channels)
# string:                  'logicle' or 'arcsinh'
tf_cofactor <- 1.
# ^ specify in case `tf` is 'logicle' (argument `w`) or 'arcsinh' (argument `b`)
notrans <- c('^FSC', '^SSC', '^Time')
# ^ patterns for channels that should not be transformed (can be empty or NULL)

signal_limits <- c(-1.0, 4.5)

## Unstained sample handling ----

fpath_fcs_unstained <- '00_UnstainedSamples'
# ^ path to samples without antibodies per batch, named "Batch{NUMBER}.fcs", or
#   NULL if they are not available

## Normalisation specification ----

### Basic batch effect correction specs ----
normalise <- TRUE # whether to normalise
batch_remove <- c(1, 2, 3, 4) # which batches (if any) to remove from analysis
fname_norm_model_precomputed <- NULL
# ^ if specified, the section 'Normalisation training set' below is irrelevant

### Normalisation training set ----
fpath_fcs_normtrain <- NULL
# ^ path to files to use for training a new normalisation model (matched
#   via numeric part of filename) or NULL if the pre-processed input files
#   should be used
train_norm_on_qc <- FALSE
# ^ use QC files for training, or downsampled aggregates per batch

### Reference vs. target batch(es) ----
# batch_ref  <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
# batch_norm <- c(1, 2, 3, 4)
batch_ref  <- c(5)
batch_norm <- c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
per_batch  <- TRUE
# ^ train separate norm model for each target batch or correct jointly

### Spline-fitting hyperparameters ----
n_quantiles     <- 15
quantile_values <- c(0.0001, ((1:(n_quantiles-1))/(n_quantiles-1))[-(n_quantiles-1)], 0.9999)
quantile_limits <- NULL


### Breaks for histogram to compute EMD ----
emd_breaks <- seq(from = min(signal_limits), to = max(signal_limits), length.out = 10)
# ^ for computing differences between batches based on binned signal (more
#   conservative)

### Normalisation FlowSOM (meta)clustering ----
norm_xdim <- 12
norm_ydim <- 12
norm_n_metaclusters <- 5
# ^ fewer than for the feature-extraction FlowSOM model

### CytoNorm train/transform markers ----
train_on_lin    <- TRUE # use lineage markers in the norm FlowSOM
train_on_state  <- FALSE # use state markers in the norm FlowSOM
transform_lin   <- TRUE # normalise lineage markers by CytoNorm
transform_state <- TRUE # normalise state markers by CytoNorm

## Feature-extraction FlowSOM model specs ----

### Main FlowSOM clustering params ----
xdim           <- 12
ydim           <- 12
n_metaclusters <- 40

### FlowSOM feature extraction specs ----
feature_mapping_batch_size <- 100
# ^ how many FCS files to process at a time

## FlowSOM outlier sample exclusion ----

n_dev <- 20
# ^ max deviations from average metacluster percentage before sample is removed
#   as outlier
dev_type <- 'mad'
# ^ for deviations, use MADs from median ('mad') or SDs from average ('sd');
#   or NULL for no outlier sample exclusion

## Profiling ----

gating_available <- TRUE
# ^ whether a FlowJo workspace with gating on the FlowSOM aggregate is available
fname_gating_wsp <- './05_ProfilingWorkspace.wsp'

gate_names <- list(
  'DN memory B' = c(
    '/Mature B cell/Class-switched/DN memory/IgM+ Memory',
    '/Mature B cell/Class-switched/DN memory/IgM+ Memory/IgM_AM',
    '/Mature B cell/Class-switched/DN memory/IgM+ Memory/IgM_IM',
    '/Mature B cell/Class-switched/DN memory/IgM+ Memory/IgM_RM',
    '/Mature B cell/Class-switched/DN memory/IgM+ Memory/IgM_TLM',
    '/Mature B cell/Class-switched/DN memory/Other memory'
  ),
  'IgA+ memory B' = c(
    '/Mature B cell/Class-switched/IgA+ memory',
    '/Mature B cell/Class-switched/IgA+ memory/IgA_AM',
    '/Mature B cell/Class-switched/IgA+ memory/IgA_IM',
    '/Mature B cell/Class-switched/IgA+ memory/IgA_RM',
    '/Mature B cell/Class-switched/IgA+ memory/IgA_TLM'
  ),
  'IgG+ memory B' = c(
    '/Mature B cell/Class-switched/IgG+ memory',
    '/Mature B cell/Class-switched/IgG+ memory/IgG_AM',
    '/Mature B cell/Class-switched/IgG+ memory/IgG_IM',
    '/Mature B cell/Class-switched/IgG+ memory/IgG_RM',
    '/Mature B cell/Class-switched/IgG+ memory/IgG_TLM'
  ),
  'IgD+CD27+ memory B' = c(
    '/Mature B cell/IgD+CD27+/IgD-only memory'
  ),
  'IgD+CD27+ marginal-zone B' = c(
    '/Mature B cell/IgD+CD27+/MZ B cells',
    '/Mature B cell/IgD+CD27+/MZ B cells/CD21+ MZ',
    '/Mature B cell/IgD+CD27+/MZ B cells/CD21- MZ'
  ),
  'CD21+ naive B' = c(
    '/Mature B cell/Naive/CD21+ naive'
  ),
  'CD21- naive B' = c(
    '/Mature B cell/Naive/CD21- naive'
  ),
  'Transitional B' = c(
    '/Transitional B cells'
  )
)
# ^ scheme for renaming and/or merging gates using gating matrix column names
#   or indices (or NULL); if specified, unnamed columns will be dropped (!)

require(parallel)
require(foreach)
require(doParallel)
require(doSNOW)

source('99_Aux.R')

