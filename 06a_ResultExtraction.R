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

## iidx internal module: 06a_ResultExtraction.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Formats results of differential expression analysis to make them
# readable.


## Function: collect DE results in a data frame ----

get_stats <- function(
  res,                   # DA, DS-MFI or DS-Pheno results
  predictor,             # biological predictor of interest
  confounder    = NULL,  # biological confounder of interest
  metacluster   = NULL,  # metacluster of interest (if not all)
  marker        = NULL,  # state marker of interest (if not all)
  robustness    = NULL,  # robustness info
  state_markers = NULL,  # all state markers to make sure DS is limited to them
  wide_format   = TRUE
) {
  
  de_results(
    res           = res,
    predictor     = predictor,
    confounder    = confounder,
    metacluster   = metacluster,
    marker        = marker,
    robustness    = robustness,
    state_markers = state_markers,
    wide_format   = wide_format
  )
}


