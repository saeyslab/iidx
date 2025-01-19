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
# hyperparameters. The function of each input parameter is explained in the
# README.md file of the iidx repository.

# Input parameter values for our T- and NK-cell analysis are specified here.


library(tidyverse)

## Global ----

run_parallel <- TRUE
fpath_fcs <- '00_StainedSamples'
fpath_fcs_unstained <- '00_UnstainedSamples'
annotation <- read.csv('Annotation.csv')
batch_remove <- c(1, 2, 3, 4)

channels <- c(
  'FSC-A'  = 'FSC-A',             
  'FSC-H'  = 'FSC-H',             
  'SSC-A'  = 'SSC-A',             
  'B515-A' = 'TCR Vd1',           
  'B610-A' = 'CD57',              
  'B660-A' = 'CD244',             
  'B710-A' = 'CD127',             
  'B780-A' = 'CXCR5',             
  'V450-A' = 'CD122',             
  'V510-A' = 'CD3',               
  'V570-A' = 'CD8a',              
  'V605-A' = 'CD158',             
  'V655-A' = 'CD28',              
  'V710-A' = 'TCR Va7_2',         
  'V750-A' = 'CD38',              
  'V785-A' = 'CD27',              
  'U390-A' = 'CCR7',              
  'U450-A' = 'Viability',         
  'U500-A' = 'CD16',              
  'U570-A' = 'CD56',              
  'U660-A' = 'PD-1',              
  'U740-A' = 'CD95',              
  'U785-A' = 'CD4',               
  'R670-A' = 'CD1d-PBS57 tetramer',
  'R730-A' = 'CD45RA',            
  'R780-A' = 'CCR5',              
  'G575-A' = 'TCR Vg9',           
  'G610-A' = 'TCR Vd2',           
  'G660-A' = 'CD161',             
  'G710-A' = 'HLA-DR',            
  'G780-A' = 'CD314',             
  'Time'   = 'Time'               
)
idcs_channels_lineage <-
  c(4, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 24, 25, 27, 28, 29, 30)
idcs_channels_state   <- c(5, 6, 7, 8, 9, 12, 15, 21, 26, 31)

## 00_Preprocessing.R ----

fname_wsp <- NULL
fpath_fcs_ref <- NULL
channels_remove <- c(
  'remove_from_FS_FM',
  'FJ Sample QC'
)
channels_rename <- c( # from - to
  'U395-A' = 'U390-A', 
  'U670-A' = 'U660-A'
)
thresholds <- c(
  'CD3'       = 1.60,
  'CCR5'      = 1.15,
  'CCR7'      = 1.35,
  'CD1d-PBS57 tetramer' = 1.60,
  'CD4'       = 1.50,
  'CD8a'      = 1.70,
  'CD16'      = 1.50,
  'CD27'      = 1.65,
  'CD28'      = 1.45,
  'CD38'      = 1.80,
  'CD45RA'    = 1.70,
  'CD56'      = 1.70,
  'CD57'      = 1.50,
  'CD95'      = 1.50,
  'CD122'     = 1.85,
  'CD244'     = 1.50,
  'CD127'     = 1.40,
  'CD158'     = 1.85,
  'CD161'     = 1.35,
  'CD314'     = 1.70,
  'CXCR5'     = 1.60,
  'HLA-DR'    = 1.75,
  'PD-1'      = 1.50,
  'TCR Va7_2' = 1.80,
  'TCR Vd1'   = 1.65,
  'TCR Vd2'   = 2.30,
  'TCR Vg9'   = 2.40
)
fpath_subset_idcs <- '99_SubsetIdcs'
compensate <- TRUE
tf <- readRDS('TransformList.RDS')
tf_cofactor <- 1.
notrans <- c('^FSC', '^SSC', '^Time')
signal_limits <- c(-1.0, 4.2)

## 01_BatchEffect.R ----

normalise <- FALSE
fname_norm_model_precomputed <- NULL
fpath_fcs_normtrain <- NULL
train_norm_on_qc <- FALSE
batch_ref       <- c(5)
batch_norm      <-c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
per_batch       <- TRUE
n_quantiles     <- 15
quantile_values <-
  c(0.0001, ((1:(n_quantiles-1))/(n_quantiles-1))[-(n_quantiles-1)], 0.9999)
quantile_limits <- NULL
emd_breaks          <-
  seq(from = min(signal_limits), to = max(signal_limits), length.out = 10)
norm_xdim           <- 12
norm_ydim           <- 12
norm_n_metaclusters <- 5
train_on_lin        <- TRUE
train_on_state      <- FALSE
transform_lin       <- TRUE
transform_state     <- TRUE

## 02_FlowSOMFeatureExtraction.R ----

xdim           <- 12
ydim           <- 12
n_metaclusters <- 50
feature_mapping_batch_size <- 100

## 03_OutlierAndNoiseDetection.R ----

n_dev <- 15
dev_type <- 'sd'

## 04_StatisticalModelling.R ----

predictors  <- c('Sex', 'Age')
confounders <- c('Sex', 'Age')

## 05_Profiling.R ----

gating_available <- TRUE
fname_gating_wsp <- './05_ProfilingWorkspace.wsp'

gate_names <- list(
  'iNKT' = c(
    '/CD3+/iNKT'
  ),
  'Vd1+Vg9+ T' = c(
    '/CD3+/Non-NKT/Vd1+Vg9+'
  ),
  'Vd1+Vg9- T' = c(
    '/CD3+/Non-NKT/Vd1+Vg9-'
  ),
  'Vd1-Vg9+ T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9+',
    '/CD3+/Non-NKT/Vd1-Vg9+/Vg9+Vd2-',
    '/CD3+/Non-NKT/Vd1-Vg9+/Vg9+Vd2dim',
    '/CD3+/Non-NKT/Vd1-Vg9+/Vg9+Vd2high'
  ),
  'Vd1-Vg9- T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-'
  ),
  'CD4+ central memory T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4+/CD4:CM'
  ),
  'CD4+ effector memory T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4+/CD4:EM'
  ),
  'CD4+ naive T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4+/CD4:R7+RA+/CD4:Naive'
  ),
  'CD4+ stem cell-like memory T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4+/CD4:R7+RA+/CD4:TSCM'
  ),
  'CD4+ terminal effector T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4+/CD4:TE'
  ),
  'CD4+CD8+ T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4+CD8+'
  ),
  'CD4-CD8- T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD4-CD8-'
  ),
  'CD8+ central memory T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD8+/CD8:CM'
  ),
  'CD8+ effector memory T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD8+/CD8:EM'
  ),
  'CD8+ naive T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD8+/CD8:R7+RA+/CD8:Naive'
  ),
  'CD8+ stem cell-like memory T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD8+/CD8:R7+RA+/CD8:TSCM'
  ),
  'CD8+ terminal effector T' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/Conventional/CD8+/CD8:TE'
  ),
  'MAIT' = c(
    '/CD3+/Non-NKT/Vd1-Vg9-/MAIT'
  ),
  'CD3-' = c(
    '/CD3-'
  ),
  'Early NK' = c(
    '/CD3-/Early NK'
  ),
  'Mature NK' = c(
    '/CD3-/Mature NK'
  ),
  'NK1' = c(
    '/CD3-/NK1'
  ),
  'NK2' = c(
    '/CD3-/NK2'
  ),
  'Terminal NK' = c(
    '/CD3-/Terminal NK'
  )
)

## Additional input parameters and functions ----

source('99_Aux.R')
