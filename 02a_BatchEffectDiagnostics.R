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

## iidx internal module: 02a_BatchEffectDiagnostics.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines functions for generation of diagnostic UMAP plots for
# detection of batch effects.


## Function: resolve NAs in feature matrices ----

handle_na_values <- function(
    f_pre,        # feature matrix for pre-normalisation data
    f_post = NULL # feature matrix for post-normalisation data
) {
  
  norm <- !is.null(f_post)
  
  ## Resolve potential NA values in MFI/phenopositivity matrices; impute...
  ## ...where possible (this is only done for the purposes of UMAP...
  ## ...visualisation)
  
  nr <- nrow(f_pre)
  remove_cols <- c() # initialise vector of columns with only NAs
  for (idx_col in seq_len(ncol(f_pre))) {
    
    n <- c(
      sum(is.na(f_pre[, idx_col]))/nr, # proportion of NAs in pre-norm
      if (norm) {
        sum(is.na(f_post[, idx_col]))/nr # proportion of NA in post-norm
      } else {
        c()
      }
    )
    if (1. %in% n) { # entire column only NA in at least one of the matrices
      
      remove_cols <- c(remove_cols, idx_col) # mark column for removal
    } else { # entire column not fully NA for either pre- or post-norm matrix
      
      f_pre[
        is.na(f_pre[, idx_col]), idx_col
      ] <- mean(f_pre[, idx_col], na.rm = TRUE) # impute pre-norm using mean
      if (norm) {
        f_post[
          is.na(f_post[, idx_col]), idx_col
        ] <- mean(f_post[, idx_col], na.rm = TRUE) # impute post-norm using mean
      }
    }
  }
  if (length(remove_cols)>0) { # columns marker for removal
    f_pre  <- f_pre[, -remove_cols] # filtered pre-norm matrix
    if (norm) {
      f_post <- f_post[, -remove_cols] # filtered post-norm matrix
    }
  }
  
  list(
    'pre'  = f_pre,
    'post' = f_post
  )
}

## Function: evaluate batch effect using UMAP ----

viz_norm <- function(
    features_pre,         # pre-normalisation features
    fname_pdf,            # output PDF file name
    annotation,           # sample-level annotation
    pal,                  # manual colour palette
    fnames_qc     = NULL, # QC sample file names
    features_post = NULL, # post-normalisation features
    seed          = 1     # random seed for reproducibility
) {
  
  ## Resolve whether comparison between pre-norm and post-norm can be plotted
  include_postnorm <- !is.null(features_post)
  
  ## Auxiliary function for dimensionality reduction to 2d
  f_dimred <- function(x) {
    x <- umap::umap(x)$layout
    colnames(x) <- c('X', 'Y')
    x
  }    
  
  ## Embed pre-norm feature vectors in 2d
  set.seed(seed)
  l_pre  <- f_dimred(features_pre) 
  
  ## Embed post-norm feature vectors in 2d
  l_post <- NULL
  if (include_postnorm) { # post-norm data available
    set.seed(seed)
    l_post <- f_dimred(features_post)
  }
  
  ## Label points by batch
  l_pre <- data.frame(
    l_pre,
    'Batch' = as.factor(
      annotation$Batch[
        match(basename(rownames(features_pre)), annotation$FileName)
      ]
    )
  )
  if (include_postnorm) {
    l_post <- data.frame(
      l_post,
      'Batch' = as.factor(
        annotation$Batch[
          match(basename(rownames(features_post)), annotation$FileName)
        ]
      )
    )
  }
  
  ## Label points as pre- or post-norm, if applicable
  if (include_postnorm) {
    
    layout <- rbind(l_pre, l_post)
    cond   <- as.factor(c(
      rep('Before Normalisation', nrow(l_pre)),
      rep('After Normalisation', nrow(l_post))
    ))
    cond   <- factor(cond, levels = c(
      'Before Normalisation', 'After Normalisation'
    ))
    layout <- cbind(layout, Condition = cond)
  } else {
    
    layout <- l_pre
    layout$Condition <- as.factor('Pre-processed')
  }
  
  ## Generate combined plot with colour-coding per batch
  p_all <-
    ggplot2::ggplot(
      data    = layout,
      mapping = ggplot2::aes(x = .data$X, y = .data$Y, colour = .data$Batch)
    ) +
    ggplot2::facet_grid(
      ~.data$Condition
    ) +
    ggplot2::geom_point(
      size = 1.2, alpha = 0.85
    ) +
    ggplot2::scale_colour_manual(
      values = pal
    ) +
    ggplot2::theme_minimal()
  
  ## Generate plots with each batch highlighted separately
  batches <- sort(unique(layout$Batch))
  one_hot_labels <- lapply(
    batches,
    function(idx) as.factor(as.integer(layout$Batch==idx))
  )
  one_hot_plots <- lapply(
    batches,
    function(idx) {
      l <- data.frame(
        'X'         = layout$X,
        'Y'         = layout$Y,
        'col'       = one_hot_labels[[idx]],
        'Condition' = layout$Condition
      )
      ggplot2::ggplot(
        data    = l,
        mapping = ggplot2::aes(x = .data$X, y = .data$Y, colour = .data$col),
        alpha   = 0.95
      ) +
        ggplot2::facet_grid(
          ~.data$Condition
        ) +
        ggplot2::geom_point(
          data = l[l$col==0, ], size = .1, col = 'grey'
        ) +
        ggplot2::geom_point(
          data = l[l$col==1, ], size = .3, col = 'darkred'
        ) +
        ggplot2::ggtitle(
          paste0('Batch ', idx)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position  = 'none',
          strip.background = ggplot2::element_rect(colour = 'black', fill = NA),
          plot.title       = ggplot2::element_text(size = 10),
          axis.text.x      = ggplot2::element_blank(),
          axis.text.y      = ggplot2::element_blank(),
          axis.title.x     = ggplot2::element_blank(),
          axis.title.y     = ggplot2::element_blank(),
          axis.ticks       = ggplot2::element_blank()
        )
    })
  
  ## Highlight batch-wise QC samples if all given
  if (
    !is.null(fnames_qc)&&
    length(fnames_qc)==length(batches)
  ) {
    
    ## Generate labels for QC points by batch
    labs_qc <- rep(0, times = nrow(l_pre))
    labs_qc[
      basename(rownames(features_pre)) %in% basename(fnames_qc)
    ] <- batches
    
    labs_qc <- as.factor(labs_qc)
    levels(labs_qc)[1] <- 'not QC' # (rest in grey)
    
    ## Generate combined plot with colour-coding of QC samples per batch
    p <-
      ggplot2::ggplot(
        data    = layout,
        mapping = ggplot2::aes(x = .data$X, y = .data$Y, colour = labs_qc)
      ) +
      ggplot2::facet_grid(
        ~.data$Condition
      ) +
      ggplot2::geom_point(
        size = 1.0
      ) +
      ggplot2::scale_colour_manual(
        values = c('grey', pal)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        colour = 'QC batch label'
      ) +
      ggplot2::geom_point(
        data = layout[
          rep(
            labs_qc!='not QC',
            times =
              if (include_postnorm) {
                2 # double amount of points if post-norm embedding included
              } else {
                1
              } 
          ), ,
          drop = FALSE
        ],
        colour = pal[
          rep(
            labs_qc[labs_qc!='not QC'],
            times = if (include_postnorm) { 2 } else { 1 }
          )
        ],
        size = 2.0
      )
  } else {
    fnames_qc <- NULL
  }
  
  ## Save all plots in PDF
  pdf(fname_pdf, width = 7, height = 5.5)
  plot(p_all)
  if (!is.null(fnames_qc)) {
    plot(p)
  }
  for (idx in seq_along(one_hot_plots)) {
    plot(one_hot_plots[[idx]])
  }
  dev.off()
}


