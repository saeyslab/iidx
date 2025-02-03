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

# Input parameter values for our B-cell analysis are specified here.


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
  'B515-A' = 'RSV Cav1', 
  'B610-A' = 'CD141',    
  'B660-A' = 'CD123',    
  'B710-A' = 'CD16',     
  'B780-A' = 'IgD',      
  'G575-A' = 'CD32',     
  'G610-A' = 'CD40',     
  'G660-A' = 'CD85j',    
  'G710-A' = 'CD11c',    
  'G780-A' = 'CXCR3',    
  'R670-A' = 'IgA',      
  'R730-A' = 'CD27',     
  'R780-A' = 'CD19',     
  'U390-A' = 'CD1c',     
  'U450-A' = 'Viability',
  'U500-A' = 'CD21',     
  'U570-A' = 'TACI',     
  'U660-A' = 'HLA-DR',   
  'U740-A' = 'IgG',      
  'U785-A' = 'CD20',     
  'V450-A' = 'IL-21R',   
  'V510-A' = 'CD14',     
  'V570-A' = 'IgM',      
  'V605-A' = 'BAFF-R',   
  'V655-A' = 'CD10',     
  'V710-A' = 'CD23',     
  'V750-A' = 'CXCR5',    
  'V785-A' = 'CD64',     
  'Time'   = 'Time'      
)
idcs_channels_lineage <- c(8, 14, 15, 16, 19, 21, 22, 23, 26, 28)
idcs_channels_state <- c(9, 10, 11, 12, 13, 17, 20, 24, 27, 29, 30)

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
  'IgD'    = 1.70,
  'CD1c'   = 1.85,
  'CD21'   = 1.40,
  'HLA-DR' = 2.00,
  'IgG'    = 1.70,
  'IL-21R' = 1.50,
  'TACI'   = 1.50,
  'IgM'    = 1.50,
  'BAFF-R' = 2.30,
  'CD10'   = 1.80,
  'CD23'   = 1.60,
  'CD32'   = 1.90,
  'CD40'   = 2.30,
  'CD85j'  = 1.55,
  'CD11c'  = 1.75,
  'CXCR3'  = 2.00,
  'CXCR5'  = 1.60,
  'IgA'    = 2.20,
  'CD27'   = 1.75,
  'CD19'   = 1.15
)
fpath_subset_idcs <- NULL
compensate <- TRUE
tf <- readRDS('TransformList.RDS')
tf_cofactor <- 1.
notrans <- c('^FSC', '^SSC', '^Time')
signal_limits <- c(-1.0, 4.5)

## 01_BatchEffect.R ----

normalise <- TRUE
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

n_dev <- NULL
dev_type <- 'mad'

## 04_StatisticalModelling.R ----

predictors  <- c('Sex', 'Age')
confounders <- c('Sex', 'Age')

## 05_Profiling.R ----

gating_available <- TRUE
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

## Additional input parameters and functions ----

source('99_Aux.R')
