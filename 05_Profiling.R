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

## iidx user-level module: 05_Profiling.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Creates profiles (descriptions) of FlowSOM metasclusters, using
# manually defined population gates and distributions of signal per marker.


message(
  '## iidx: 05_Profiling.R STARTED ',
  as.character(round(Sys.time()))
)

## Get parameters and auxiliary functions ----

input_files_needed <- FALSE
source('99_SetUp.R')

## Create results directory ----

if (!file.exists(fpath_res05)) {
  dir.create(fpath_res05)
}

## Gather aggregate data, labels and features ----

inputs <- prep_inputs(
  fname_features_postnorm, fname_sample_outliers, annotation, batch_remove
)
agg  <- flowCore::read.FCS(fname_agg)  # aggregate expression data
fsom <- readRDS(fname_fsom)            # trained FlowSOM model
mcs  <- FlowSOM::GetMetaclusters(fsom) # metacluster labels per row of `agg`

if (gating_available) { # FlowJo workspace with gating hierarchy available
  
  ## Map all gates onto aggregate expression data ----
  
  ## Parse FlowJo workspace
  ws <- CytoML::open_flowjo_xml(fname_gating_wsp)
  gs <- CytoML::flowjo_to_gatingset(
    ws,
    path = fpath_res02, # use aggregate FCS as reference sample for parsing
    name = 1
  )
  
  ## Prepare expression data for applying the gating
  gs_data <- flowWorkspace::GatingSet(flowCore::flowSet(agg))
  
  ## Extract names (paths) of all gates
  gates <- flowWorkspace::gs_get_pop_paths(gs[[1]])[-1]
  
  ## Iterate over gates
  for (gate in gates) {
    
    ## Resolve the parent gate
    pre <- gsub('/[^/]*$', '', gate)
    if (pre=='') {
      pre <- 'root'
    }
    
    ## Extract the geometric definition of gate
    def <- flowWorkspace::gh_pop_get_gate(
      obj = gs[[1]],
      y   = gate
    )
    
    ## Add gate to the hierarchy
    flowWorkspace::gs_pop_add(
      gs     = gs_data,
      gate   = def,
      parent = pre
    )
  }
  
  ## Resolve gate memberships of each cell ----
  
  suppressMessages({
    flowWorkspace::recompute(x = gs_data, y = 'root')
  })
  
  ## Compute matrix of Boolean assignments ----
  
  ## Get matrix for all nodes of gating tree
  gm <- flowWorkspace::gh_pop_get_indices_mat(
    gh = gs_data,
    y  = gates
  )
  
  ## Resolve renamings and mergers of gates
  if (!is.null(gate_names)) {
    gm <- do.call(
      cbind,
      lapply(gate_names, function(g) {
        rowSums(gm[, g, drop = FALSE])>0
      }))
  }
  
  ## Save the gating matrix
  saveRDS(gm, fname_agg_gm)
  
  ## Create vector of labels per cell ----
  
  ## Extract labels from matrix, and 'Unlabeled' for cells outside all gates
  agg_labs <- FlowSOM::ManualVector(
    manualMatrix = gm,
    cellTypes    = colnames(gm)
  )
  
  ## Save as factor vector
  saveRDS(agg_labs, fname_agg_labels)
  
  ## Compute metacluster compositions from labels ----
  
  ## Compute cell counts per metacluster per label
  mc_prof <- do.call(
    rbind,
    lapply(
      levels(mcs),
      function(mc) table(agg_labs[mcs==mc])
    )
  )
  rownames(mc_prof) <- levels(mcs)
  
  ## Indicate clustering model for forward compatibility
  attributes(mc_prof)$BaseModel <- 'FlowSOM'
  
  ## Save all metacluster compositions (profiles)
  saveRDS(mc_prof, fname_mc_profiles)
}

## Compute signal quartiles across metaclusters ----

## Get expression data with markers (not channels) as column names
cn   <- flowCore::colnames(agg) # channels
mask <- cn%in%names(channels)   # channels used in analysis
cn[mask] <- channels[cn[mask]]  # rename channels to markers
flowCore::colnames(agg) <- cn   # use as column names
d <- agg@exprs # extract expression data as matrix
markers <- as.vector(colnames(d))

## Load metacluster labels per cell
mc      <- readRDS(fname_metaclustering)
mcs     <- levels(mc)

## Initialise signal quartiles per metacluster
mc_quart <- vector(mode = 'list', length = length(mcs))

## Iterate over metaclusters
for (idx_mc in seq_along(mcs)) {
  
  this_mc <- mcs[idx_mc]
  
  ## Compute quartiles of metacluster-specific signal
  x <- apply(
    X = d[
      mc==this_mc,
      channels[c(idcs_channels_lineage, idcs_channels_state)],
      drop = FALSE
    ],
    MARGIN = 2,
    FUN = function(dd) {
      quantile(
        x     = dd,
        probs = c(.25, .50, .75)
      )
    }
  )
  mc_quart[[idx_mc]] <- data.frame(
    'Compartment' =  paste0('MC', this_mc),
    'Marker'      = rep(colnames(x), each = 3),
    'Quantile'    = c(.25, .50, .75), 
    'Value'       = as.vector(x)
  )
}

## Bind together in a single data frame
mc_quart <- do.call(rbind, mc_quart)

## Compute quartiles of overall signal
x <- apply(
  X = d[
    ,
    channels[c(idcs_channels_lineage, idcs_channels_state)],
    drop = FALSE
  ],
  MARGIN = 2,
  FUN = function(dd) {
    quantile(
      x     = dd,
      probs = c(.25, .50, .75)
    )
  }
)
agg_quart <- data.frame(
  'Compartment' = 'All',
  'Marker'      = rep(colnames(x), each = 3),
  'Quantile'    = c(.25, .50, .75), 
  'Value'       = as.vector(x)
)

## Bind metacluster-specific and overall signal together
mc_quart <- rbind(agg_quart, mc_quart)
rownames(mc_quart) <- NULL

if (unstained_samples_available) { # unstained FCS files were processed
  
  message('Extracting background signal from unstained samples')
  
  ## Extract aggregate unstained expression data
  agg_u    <- flowCore::read.FCS(fname_agg_unstained)
  cn       <- flowCore::colnames(agg_u) # channels
  mask     <- cn%in%names(channels)     # channels used in analysis
  cn[mask] <- channels[cn[mask]]        # rename channels to markers
  flowCore::colnames(agg_u) <- cn       # use as column names
  d_u      <- agg_u@exprs               # extract expression data as matrix
  
  ## Compute quartiles of unstained signal
  x <- apply(
    X = d_u[
      ,
      channels[c(idcs_channels_lineage, idcs_channels_state)],
      drop = FALSE
    ],
    MARGIN = 2,
    FUN = function(dd) {
      quantile(
        x     = dd,
        probs = c(.25, .50, .75)
      )
    }
  )
  agg_u_quart <- data.frame(
    'Compartment' = 'Unstained',
    'Marker'      = rep(colnames(x), each = 3),
    'Quantile'    = c(.25, .50, .75), 
    'Value'       = as.vector(x)
  )
  
  ## Bind stained and unstained signal quartile values together
  mc_quart <- rbind(agg_u_quart, mc_quart)
}

## Convert metacluster and marker flags to factors for easier plotting
mc_quart$Compartment <- as.factor(mc_quart$Compartment)
mc_quart$Marker      <- as.factor(mc_quart$Marker)

## Save signal quartile data
saveRDS(mc_quart, fname_mc_quartiles)

message(
  '## iidx: 05_Profiling.R FINISHED ',
  as.character(round(Sys.time()))
)
