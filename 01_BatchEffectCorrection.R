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

## iidx user-level module: 01_BatchEffectCorrection.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Uses CytoNorm normalisation to reduce batch effects in FCS data.

message(
  '## iidx: 01_BatchEffectCorrection.R STARTED ',
  as.character(round(Sys.time()))
)

## Get parameters and auxiliary functions ----

input_files_needed <- TRUE
source('99_SetUp.R')
source('01a_Normalisation.R')

## Create results directory ----

if (!file.exists(fpath_res01)) {
  dir.create(fpath_res01)
}

## Enumerate data cohorts and batches

cohorts <- unique(annotation$Cohort)
batches <- sort(unique(annotation$Batch))

## Enumerate channels and markers used in analysis

ch <- channels[c(idcs_channels_lineage, idcs_channels_state)]
m  <- as.vector(ch)
ch <- names(ch)

## Create directory for normalised FCS files
if (!file.exists(fpath_fcs_norm)) {
  dir.create(fpath_fcs_norm)
}

## Resolve batch per sample
idcs_ref              <- annotation$Batch %in% batch_ref
idcs_norm             <- annotation$Batch %in% batch_norm

## Resolve reference-batch and target-batch FCS files
fnames_ref_comptrans  <- fnames_fcs_comptrans[idcs_ref]
fnames_norm_comptrans <- fnames_fcs_comptrans[idcs_norm]
batch_ref_perfile     <- annotation$Batch[idcs_ref]
batch_norm_perfile    <- annotation$Batch[idcs_norm]

## Resolve QC files if present
fnames_qc             <- file.path(
  fpath_fcs_comptrans,
  annotation$FileName[annotation$QC]
)
batches_qc            <- annotation$Batch[annotation$QC]

## Create small pre-norm batch aggregates ----

message('Sampling representative data per batch')

## Create directory for batch aggregates
if (!file.exists(fpath_agg_batchwise)) {
  dir.create(fpath_agg_batchwise)
}

pb <- utils::txtProgressBar(
  min = 0, max = length(batch_ref)+length(batch_norm), style = 3
) # progress bar
counter <- 1
set.seed(1)

## Iterate over reference batch(es)
for (idx_batch in batch_ref) {
  
  str_batch <- formatC(idx_batch, width = 2, format = 'd', flag = '0')
  idcs      <- which(annotation$Batch == idx_batch)
  fname_ref <- file.path( # output filepath
    fpath_agg_batchwise,
    paste0('Aggregate_Ref_Batch', str_batch, '.fcs')
  )
  
  ## Aggregate files in batch
  ff <- FlowSOM::AggregateFlowFrames(
    fileNames = fnames_fcs_comptrans[idcs],
    cTotal    = 1e4,
    channels  = ch,
    silent    = TRUE
  )
  flowCore::write.FCS(x = ff, filename = fname_ref)
  utils::setTxtProgressBar(pb, counter)
  counter <- counter + 1
}

## Iterate over target batch(es)
for (idx_batch in batch_norm) {
  
  str_batch     <- formatC(idx_batch, width = 2, format = 'd', flag = '0')
  fname_prenorm <- file.path( # pre-norm output filepath
    fpath_agg_batchwise,
    paste0('Aggregate_PreNorm_Batch', str_batch, '.fcs')
  )
  idcs <- which(
    annotation$Batch[
      match(basename(fnames_fcs_comptrans), annotation$FileName)
    ] == idx_batch
  )
  
  ## Aggregate files in batch
  ff <- FlowSOM::AggregateFlowFrames(
    fileNames = fnames_fcs_comptrans[idcs],
    cTotal    = 1e4,
    channels  = ch,
    silent    = TRUE
  )
  flowCore::write.FCS(x = ff, filename = fname_prenorm)
  utils::setTxtProgressBar(pb, counter)
  counter <- counter + 1
}
close(pb)

## Compute and plot pre-norm EMDs ----

## Initialise EMD values for batch pairs and heatmaps per marker
markers_to_plot  <- m
channels_to_plot <- ch
nm     <- length(markers_to_plot)
d_pre  <- vector(mode = 'list', length = nm)
p_pre  <- vector(mode = 'list', length = nm)

## Iterate over markers
counter <- 1
for (idx in seq_along(markers_to_plot)) {
  
  this_channel <- channels_to_plot[idx]
  this_marker  <- markers_to_plot[idx]
  
  ## Load data from small batch aggregates
  data <- load_batch_aggregate_data_for_marker(
    channel             = this_channel,
    marker              = this_marker,
    fpath_agg_batchwise = fpath_agg_batchwise,
    batch_norm          = batch_norm,
    batch_ref           = batch_ref,
    annotation          = annotation,
    include_postnorm    = FALSE # normalised data does not exist at this point
  )
  
  message(
    'Computing EMD between batches for ', this_marker,
    ' (', idx, '/', length(markers_to_plot), ')'
  )
  
  ## Compute EMD matrix for each pair of batches
  emd <- cross_batch_emd(
    data             = data,
    batch_ref        = batch_ref,
    emd_breaks       = emd_breaks,
    include_postnorm = FALSE
  )
  emd_prenorm <- emd$prenorm
  emd_prenorm[upper.tri(emd_prenorm)] <- NA
  # ^ make distance matrix lower triangular
  d_pre[[counter]]  <- emd_prenorm
  
  ## Resolve colour scale for plotting
  ymax   <- max(emd_prenorm, na.rm = TRUE)
  breaks <- seq(from = 0, to = ymax, by = ymax/100)
  
  ## Generate heatmap from the distance matrix
  plot_emd_pre <- pheatmap::pheatmap(
    emd_prenorm,
    main            = paste0(this_marker, ' (pre-norm)'),
    cluster_rows    = FALSE,
    cluster_cols    = FALSE,
    breaks          = breaks,
    display_numbers = TRUE,
    fontsize        = 6,
    angle_col       = 0,
    legend          = TRUE,
    border_color    = NA,
    na_col          = NA,
    silent          = TRUE
  )
  attributes(plot_emd_pre)$Channel <- this_channel
  attributes(plot_emd_pre)$Marker  <- this_marker
  
  p_pre[[counter]]  <- plot_emd_pre
  counter <- counter + 1
}

## Save PDF of all EMD heatmaps
tryCatch(expr = {dev.off()}, error = function(e) {}) # ensure no open device
pdf(fname_viz_emd_prenorm, width = 10, height = 10)
for (idx in seq_along(p_pre)) { # iterate over heatmaps
  print(p_pre[[idx]])
  if (idx < length(p_pre)) {
    plot.new()
  }
}
dev.off()

if (normalise) {
  
  ## Normalise compensated & transformed FCS files ----
  
  ## Set up input parameters for training normalisation model
  train_params <- list(
    'batches'               = batches,
    'batch_ref'             = batch_ref,
    'batch_norm'            = batch_norm,
    'batch_ref_flag'        = batch_ref_flag,
    'annotation'            = annotation,
    'cohorts'               = cohorts,
    'channels'              = channels,
    'idcs_channels_lineage' = idcs_channels_lineage,
    'idcs_channels_state'   = idcs_channels_state,
    'train_on_qc'           = train_norm_on_qc,
    'fnames_qc'             = fnames_qc,
    'cells_per_agg'         = 1e6,
    'per_batch'             = per_batch,
    'quantile_values'       = quantile_values,
    'quantile_limits'       = quantile_limits,
    'norm_xdim'             = norm_xdim,
    'norm_ydim'             = norm_ydim,
    'norm_n_metaclusters'   = norm_n_metaclusters,
    'seed'                  = 1
  )
  
  ## Set up parameters for applying trained normalisation model
  transform_params <- list(
    'batches'             = batches,
    'annotation'          = annotation,
    'cohorts'             = cohorts,
    'fpath_fcs_comptrans' = fpath_fcs_comptrans,
    'batch_ref'           = batch_ref,
    'batch_norm'          = batch_norm
  )
  
  ## Load pre-computed normalisation model or train a new one
  if (!is.null(fname_norm_model_precomputed)) {
    
    norm_model     <- readRDS(fname_norm_model_precomputed)
  } else {
    
    ## Train normalisation model
    norm_model <- do.call(
      train_norm_model,
      c(
        train_params,
        list(
          'fname_norm_model'     = fname_norm_model,
          'fpath_norm_model'     = fpath_norm_model,
          'fpath_training_files' = fpath_fcs_normtrain, # (can be NULL)
          'train_on_lin'         = train_on_lin,
          'train_on_state'       = train_on_state,
          'transform_lin'        = transform_lin,
          'transform_state'      = transform_state
        )
      )
    )
  }
  
  ## Normalise previously compensated and transformed FCS files
  do.call(
    apply_norm_model,
    c(
      transform_params,
      list(
        'norm_model'         = norm_model,
        'fpath_output'       = fpath_fcs_norm
      )
    )
  )
  
  ## Get names and corresponding batches of normalised FCS files
  
  fnames_fcs_norm <- file.path(
    fpath_fcs_norm, basename(fnames_fcs_comptrans)
  )
  fnames_fcs_norm <- fnames_fcs_norm[
    file.exists(fnames_fcs_norm)
  ]
  batches <- annotation$Batch[
    match(basename(fnames_fcs_norm), annotation$FileName)
  ]
  
  ## Create small post-norm batch aggregates ----
  
  message('Sampling normalised data per target batch')
  
  pb <- utils::txtProgressBar(
    min = 0, max = length(batch_norm), style = 3
  ) # progress bar
  counter <- 1
  set.seed(1)
  
  ## Iterate over target batch(es)
  for (idx_batch in batch_norm) {
    
    str_batch <- formatC(idx_batch, width = 2, format = 'd', flag = '0')
    idcs <- which(
      annotation$Batch[
        match(basename(fnames_fcs_norm), annotation$FileName)
      ] == idx_batch
    )
    fname_postnorm <- file.path(
      fpath_agg_batchwise,
      paste0('Aggregate_PostNorm_Batch', str_batch, '.fcs')
    ) # output filepath
    
    ## Aggregate files in batch
    ff <- FlowSOM::AggregateFlowFrames(
      fileNames = fnames_fcs_norm[idcs],
      cTotal    = 1e4,
      channels  = ch,
      silent    = TRUE
    )
    flowCore::write.FCS(x = ff, filename = fname_postnorm)
    utils::setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  close(pb)
  
  ## Compute and plot pre- and post-norm EMDs ----
  
  markers_to_plot  <- m
  channels_to_plot <- ch
  
  ## (Pre-norm EMDs recomputed to adjust colour scale if needed; the pre-norm
  ## EMD plots PDF is overwritten, in case the scale changes)
  
  ## Initialise EMD values for batch pairs and heatmaps per marker
  nm     <- length(markers_to_plot)
  p_pre  <- vector(mode = 'list', length = nm)
  p_post <- vector(mode = 'list', length = nm)
  d_pre  <- vector(mode = 'list', length = nm)
  d_post <- vector(mode = 'list', length = nm)
  
  ## Iterate over markers
  counter <- 1
  for (idx in seq_along(markers_to_plot)) {
    
    this_channel <- channels_to_plot[idx]
    this_marker  <- markers_to_plot[idx]
    
    ## Load data from small batch aggregates
    data <- load_batch_aggregate_data_for_marker(
      channel             = this_channel,
      marker              = this_marker,
      fpath_agg_batchwise = fpath_agg_batchwise,
      batch_norm          = batch_norm,
      batch_ref           = batch_ref,
      annotation          = annotation,
      include_postnorm    = TRUE # get both pre- and post-norm data
    )
    
    message(
      'Computing signal changes for ', this_marker,
      ' (', idx, '/', length(markers_to_plot), ')'
    )
    
    ## Compute pre- and post-norm EMD matrices for each pair of batches
    emd <- cross_batch_emd(
      data             = data,
      batch_ref        = batch_ref,
      emd_breaks       = emd_breaks,
      include_postnorm = TRUE
    )
    emd_prenorm <- emd$prenorm
    emd_prenorm[upper.tri(emd_prenorm)] <- NA
    # ^ make pre-norm distance matrix lower triangular
    d_pre[[counter]]  <- emd_prenorm
    emd_postnorm <- emd$postnorm
    emd_postnorm[upper.tri(emd_postnorm)] <- NA
    # ^ make post-norm distance matrix lower triangular
    d_post[[counter]]  <- emd_postnorm
    
    ## Resolve colour scale for plotting
    ymax   <- max(emd_prenorm, emd_postnorm, na.rm = TRUE)
    breaks <- seq(from = 0, to = ymax, by = ymax/100)
    
    ## Generate heatmaps from the distance matrices
    plot_emd_pre <- pheatmap::pheatmap(
      emd_prenorm,
      main            = paste0(this_marker, ' (pre-norm)'),
      cluster_rows    = FALSE,
      cluster_cols    = FALSE,
      breaks          = breaks,
      display_numbers = TRUE,
      fontsize        = 6,
      angle_col       = 0,
      legend          = TRUE,
      border_color    = NA,
      na_col          = NA,
      silent          = TRUE
    )
    attributes(plot_emd_pre)$Channel <- this_channel
    attributes(plot_emd_pre)$Marker  <- this_marker
    
    plot_emd_post <- pheatmap::pheatmap(
      emd_postnorm,
      main            = paste0(this_marker, ' (post-norm)'),
      cluster_rows    = FALSE,
      cluster_cols    = FALSE,
      breaks          = breaks,
      display_numbers = TRUE,
      fontsize        = 6,
      angle_col       = 0,
      legend          = TRUE,
      border_color    = NA,
      na_col          = NA,
      silent          = TRUE
    )
    attributes(plot_emd_post)$Channel <- this_channel
    attributes(plot_emd_post)$Marker  <- this_marker
    
    p_pre[[counter]]  <- plot_emd_pre
    p_post[[counter]] <- plot_emd_post
    counter <- counter + 1
  }
  
  ## Save PDF of all pre-norm EMD heatmaps
  tryCatch(expr = {dev.off()}, error = function(e) {}) # ensure no open device
  pdf(fname_viz_emd_prenorm, width = 10, height = 10)
  for (idx in seq_along(p_pre)) {
    
    print(p_pre[[idx]])
    if (idx < length(p_pre)) {
      plot.new()
    }
  }
  dev.off()
  
  ## Save PDF of all post-norm EMD heatmaps
  tryCatch(expr = {dev.off()}, error = function(e) {}) # ensure no open device
  pdf(fname_viz_emd_postnorm, width = 10, height = 10)
  for (idx in seq_along(p_post)) {
    
    print(p_post[[idx]])
    if (idx < length(p_post)) {
      plot.new()
    }
  }
  dev.off()
  
  ## Compare signal distribution pre- and post-norm -----
  
  message('Plotting signal distributions before and after normalisation')
  
  markers_to_plot  <- m
  channels_to_plot <- ch
  nm               <- length(markers_to_plot)
  
  tryCatch(expr = {dev.off()}, error = function(e) {}) # ensure no open device
  pdf(file = fname_viz_signal_distributions)
  
  ## Iterate over markers
  p <- vector(mode = 'list', length = nm)
  for (idx in seq_along(markers_to_plot)) {
    
    this_channel <- channels_to_plot[idx]
    this_marker  <- markers_to_plot[idx]
    
    ## Load data from small batch aggregates
    data <- load_batch_aggregate_data_for_marker(
      channel             = this_channel,
      marker              = this_marker,
      fpath_agg_batchwise = fpath_agg_batchwise,
      batch_norm          = batch_norm,
      batch_ref           = batch_ref,
      annotation          = annotation,
      include_postnorm    = TRUE
    )
    
    ## Determine quantiles to highlight for easy comparison
    quantiles_to_plot <- seq(from = 0.1, to = 0.9, by = 0.2)
    
    ## Get quantile values for reference, pre-norm target and post-norm target
    q_ref <- quantile(
      data$Signal[data$Type=='Reference'],
      probs = quantiles_to_plot
    )
    q_pre <- quantile(
      data$Signal[data$Type=='Pre-norm'],
      probs = quantiles_to_plot
    )
    q_post <- quantile(
      data$Signal[data$Type=='Post-norm'],
      probs = quantiles_to_plot
    )
    
    ## Add quantile values to the data
    data$Quantiles <- NA
    data$Quantiles[
      data$Type=='Reference'
    ][
      1:length(quantiles_to_plot)
    ] <- q_ref
    data$Quantiles[
      data$Type=='Pre-norm'
    ][
      1:length(quantiles_to_plot)
    ] <- q_pre
    data$Quantiles[
      data$Type=='Post-norm'
    ][
      1:length(quantiles_to_plot)
    ] <- q_post
    
    subt <- paste0(
      'Target ',
      if (length(batch_norm) > 1) { 'batches' } else { 'batch' },
      ': ',
      paste(sort(batch_norm), collapse = ', ')
    ) 
    
    ## Extract colour palette associated with the loaded data
    pal <- attributes(data)$palette
    pal[1] <- 'darkblue'
    
    ## Plot signal densities
    plot_dens <-
      ggplot2::ggplot(
        data    = data,
        mapping = ggplot2::aes(x = .data$Signal, colour = .data$BatchFlag)
      ) +
      ggplot2::geom_vline(
        mapping   = aes(xintercept = data$Quantiles),
        linewidth = 0.35,
        alpha     = 0.7
      ) +
      ggplot2::geom_density(
        size = .2
      ) +
      ggplot2::scale_colour_manual(
        values = pal
      ) +
      ggplot2::facet_wrap(
        ~.data$Type,
        ncol = 1
      ) +
      ggplot2::xlim(
        signal_limits
      ) +
      ggplot2::theme_minimal() +
      ggplot2::ylab(
        'Density'
      ) +
      ggplot2::ggtitle(
        m
      ) +
      ggplot2::theme(
        legend.position = 'none',
        axis.title.y    = ggplot2::element_blank(),
        axis.text.y     = ggplot2::element_blank(),
        axis.ticks.y    = ggplot2::element_blank()
      )
    plot(plot_dens)
    
    ## Store plot (for viewing in RStudio)
    attributes(plot_dens)$Marker  <- m
    attributes(plot_dens)$Channel <- ch
    p[[idx]] <- plot_dens
  }
  dev.off()
  
  ## Plot CytoNorm model splines ----
  
  ## Simplify channel names
  nmod <- norm_model
  nmod$fsom$prettyColnames <- gsub(' .* <.*>', '', nmod$fsom$prettyColnames)
  
  ## Set x- and y-axis range for the quantile plots
  if (!is.null(signal_limits)) { # use `signal_limits` if available
    spline_lim <- signal_limits
  } else if (!is.null(quantile_limits)) { # or `quantile_limits` if available
    spline_lim <- quantile_limits
  } else { # or `quantile_values` range
    spline_lim <- c(min(quantile_values), max(quantile_values))
  }
  
  ## Widen the plotting range by some margin
  margin <- diff(spline_lim)/3
  spline_lim[1] <- spline_lim[1]-margin
  spline_lim[2] <- spline_lim[2]+margin
  
  ## Generate plot of splines for each metacluster and channel
  p_splines <- CytoNorm::plotSplines(
    model    = nmod,
    channels = nmod$clusterRes[[1]]$channels
  )
  
  ## Adjust axis limits
  p_splines <- lapply(
    p_splines,
    function(p) {p + ggplot2::xlim(spline_lim) + ggplot2::ylim(spline_lim)}
  )
  
  nc <- length(unique(p_splines[[1]]$data$channel)) # number of columns
  nr <- length(unique(p_splines[[1]]$data$cluster)) # number of rows
  
  ## Save PDF of all CytoNorm model splines
  tryCatch(expr = {dev.off()}, error = function(e) {}) # ensure no open device
  pdf(
    file = fname_viz_splines,
    width = nc*1, # multiply by bigger number for wider plots
    height = nr*1 # multiply by bigger number for taller plots
  )
  for (p in p_splines) {
    plot(p)
  }
  dev.off()
} else {
  
  message('Skipping normalisation step since it is disabled')
}

message(
  '## iidx: 01_BatchEffectCorrection.R FINISHED ',
  as.character(round(Sys.time()))
)
