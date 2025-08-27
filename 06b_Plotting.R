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

## iidx internal module: 06b_Plotting.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Generates (interactive) plots of differential expression results.


## Function: volcano plot of leads from DE results ----

plot_volcano <- function(
    res,                       # DA, DS-MFI or DS-Pheno results
    marker            = NULL,  # state marker of interest for filtering leads
    predictor         = NULL,  # biological predictor of interest
    confounder        = NULL,  # biological confounder of interest
    robustness_rates  = NULL,  # signal robustness rates
    robustness_cutoff = 0.,    # robustness cut-off for plotted points
    alpha_p           = 0.05,  # significance level for adjusted p-values
    alpha_eff         = NULL,  # cut-off for absolute value of effect coefficient
    size_factor       = 1.0,   # scaling factor for plotting
    label_hits        = TRUE,  # whether to label significant leads
    interactive       = FALSE, # whether to use interactive plot elements
    state_markers     = NULL   # all state markers
) {
  
  ## Resolve analysis type
  type       <- attributes(res)$AnalysisType
  
  ## Resolve base feature-extraction model (left in for future compatibility)
  base_model <- attributes(res)$BaseModel
  
  ## Resolve whether to apply robustness filter for plotted points
  filter_by_robustness <-
    type!='Differential Abundance' &&
    !is.null(robustness_rates) &&
    robustness_cutoff>0.
  
  ## Extract stats data frame
  if (type=='Differential Abundance') {
    
    if (is.null(confounder) || confounder=='none') {
      input <- res$joint[[predictor]]
    } else {
      input <- res$confounders_joint[[predictor]][[confounder]]
    }
  } else {
    
    if (is.null(confounder) || confounder=='none') {
      input <- res$main[[predictor]]
    } else {
      input <- res$confounders_main[[predictor]][[confounder]]
    }
  }
  
  ## Resolve whether predictor is continuous (not binary)
  cont <- attributes(res)$Continuous[predictor]
  
  ## Resolve cut-offs on log scales
  alpha_neglog10p  <- -log10(alpha_p)
  alpha_logeff     <- NULL
  if (!is.null(alpha_eff)) {
    
    alpha_logeff  <- log2(1.0+alpha_eff)
    alphas_eff    <- c(1.-alpha_eff, 1.+alpha_eff)
    alphas_logeff <- c(-alpha_logeff, alpha_logeff)
  }
  
  ## Flag non-significant results
  idcs_insignif <- input$AdjPVal>=alpha_p
  
  ## Generate labels per point (result)
  labs            <- input$Compartment
  input$AllLabels <- labs
  input$Label     <- labs
  
  ## Calculate neg log p-values
  input$neglog10p <- -log10(input$AdjPVal)
  
  ## Resolve responder variable
  resp <- if (type=='Differential Abundance') {
    'Abundance'
  } else if (type=='Differential State (MFI)') {
    'Signal median'
  } else if (type=='Differential State (Phenopositivity)') {
    'Proportion of phenopositives'
  }
  
  ## Filter out minimal-effect & lineage-marker DS-Pheno hits
  if (type=='Differential State (Phenopositivity)') {
    
    mask    <- input$neglog10p<Inf
    input   <- input[mask, , drop = FALSE]
    
    markers <- gsub('^MC[0-9]+[ ]', '', input$Compartment)
    mcs     <- gsub(' .*$', '', input$Compartment)
    comps   <- input$Compartment
    
    input   <- input[markers%in%state_markers, , drop = FALSE]
  }
  
  ## Resolve description of predictor
  if (!cont) { # binary predictor
    
    ## Specify which group is the positive
    str_predictor <- paste0(
      predictor, '(', levels(as.factor(annotation[,predictor]))[2], ')'
    )
  } else { # continuous predictor
    
    str_predictor <- predictor
  }
  
  ## Resolve description of model
  desc <- paste0(resp, ' ~ ', str_predictor)
  if (
    !is.null(confounder) && confounder!='none'
  ) {
    desc <- paste0(
      desc, ', corrected for ', confounder, ' (fixed effect)'
    )
  }
  
  ## Resolve description of batch modelling
  if (
    grepl('Differential State', type)
  ) { # modelling based on limma/Beta-GLMMs
    desc <- paste0(
      desc, '\nBatch modelled via random intercepts'
    )
  } else { # modelling based on edgeR
    desc <- paste0(
      desc, '\nBatch modelled via additional fixed effects'
    )
  }
  
  ## Resolve size of points for significant results
  idcs_signif <- input$AdjPVal<=alpha_p
  point_sizes <- ifelse(idcs_signif, 2.8, 1.4)*size_factor
  
  ## Resolve interactive/standard point plotting function
  f_point <- ggplot2::geom_point
  if (interactive) {
    f_point <- ggiraph::geom_point_interactive
  }
  
  ## Resolve effect coefficient label
  if (type=='Differential Abundance') {
    xmap <- 'logFC'
  } else if (type=='Differential State (MFI)') {
    xmap <- 'logFC'
  } else if (type=='Differential State (Phenopositivity)') {
    xmap <- 'logodds'
  }
  
  ## Determine axis limits
  xlims <- c(
    min(input[[xmap]], na.rm = TRUE), max(input[[xmap]], na.rm = TRUE)
  )
  ylims <- c(
    0, max(input$neglog10p, na.rm = TRUE)
  )
  
  ## Resolve robustness labelling and filtering
  if (!is.null(robustness_rates) && type!='Differential Abundance') {
    
    ## Align robustness values with DE results
    r <- robustness_rates$Rate[
      match(robustness_rates$Compartment, input$Compartment)
    ]
    
    ## Include robustness rates in point labels
    input$AllLabels <- paste0(
      input$AllLabels, ': ',
      round(x = r*100, digits = 2),
      '% robustness rate'
    )
    
    ## Filter results by robustness if requested
    if (filter_by_robustness) {
      input <- input[r>=robustness_cutoff, , drop = FALSE]
    }
  }
  
  ## Filter results by marker if requested
  if (
    type!='Differential Abundance'&&!is.null(marker)
  ) {
    
    markers <- gsub('^MC[0-9]+[ ]', '', input$Compartment)
    input   <- input[markers==marker, , drop = FALSE]
  }
  
  ## Generate volcano plot
  p <-
    ggplot2::ggplot(
      data    = input,
      mapping = ggplot2::aes(
        x       = .data[[xmap]],
        y       = .data$neglog10p,
        col     = .data$neglog10p,
        tooltip = .data$AllLabels,
        data_id = .data$AllLabels
      )
    ) +
    ggplot2::xlim(
      xlims
    ) +
    ggplot2::ylim(
      ylims
    ) +
    ggplot2::geom_vline( # zero-effect line
      xintercept = 0.,
      lwd        = 0.4*size_factor,
      col        = 'grey',
      linetype   = 2
    ) +
    ggplot2::scale_color_gradient(
      low  = '#7d7d7d',
      high = '#0a0054'
    ) +
    f_point( # results per compartment
    ) + 
    ggplot2::geom_hline( # significance level line
      yintercept = alpha_neglog10p,
      lwd        = 0.4*size_factor
    )
  
  if (!is.null(alpha_eff)) { # indicate effect size cut-off
    p <- p +
      ggplot2::geom_vline(
        xintercept = alphas_logeff,
        lwd        = 0.4*size_factor
      )
  }
  
  if (label_hits) { # label points
    p <- p +
      ggrepel::geom_label_repel(
        mapping = ggplot2::aes(
          label = .data$Label
        ),
        alpha        = 0.8,
        force        = 10,
        max.overlaps = 15
      )
  }
  
  ## Resolve x-axis label
  xlab <-
    if (type=='Differential Abundance') {
      bquote(log[2]*' fold change')
    } else if (type=='Differential State (MFI)') {
      bquote(log[2]*' fold change')
    } else if (type=='Differential State (Phenopositivity)') {
      bquote(ln*' odds')
    }
  
  ## Add labels
  p <- p +
    ggplot2::labs(
      x        = xlab,
      y        = bquote('-'*log[10]*' p'),
      title    = paste0(type, ' results by compartment'),
      subtitle = desc
    ) +
    ggplot2::theme_bw(
    ) +
    ggplot2::theme(
      legend.position = 'none',
      plot.title      = ggplot2::element_text(size = 12),
      plot.subtitle   = ggplot2::element_text(size = 10),
      axis.title.x    = ggplot2::element_text(size = 10),
      axis.title.y    = ggplot2::element_text(size = 10)
    )
  
  p
}

## Function: plot feature values against predictor ----

plot_single_hit <- function(
    res,                        # DA, DS-MFI or DS-Pheno results
    mc,                         # metacluster of interest
    annotation,                 # sample-level annotation
    data,                       # feature values
    marker             = NULL,  # marker of interest (for DS)
    noise              = NULL,  # noise cut-offs (for DS-MFI)
    robustness_rates   = NULL,  # robustness rate (for DS)
    show_noise_cutoffs = FALSE, # whether to plot noise cut-offs (for DS-MFI)
    predictor          = NULL,  # predictor of interest
    confounder         = NULL,  # confounder of interest
    view               = c(     # plot view
      'By batch',
      'By cohort',
      'Smooth',
      'Binned'
    ),
    n_bins             = 50,    # number of bins if view is 'Binned'
    lim_quantiles      = NULL,  # outcome quantile limits to handle extremes
    interactive        = FALSE  # whether to use interactive plot elements
) {
  
  ## Resolve analysis type
  type       <- attributes(res)$AnalysisType
  
  ## Resolve base feature-extraction model (left in for future compatibility)
  base_model <- attributes(res)$BaseModel
  
  ## Resolve whether predictor is continuous (not binary)
  cont  <- attributes(res)$Continuous[predictor]
  
  ## Resolve data view
  view <- view[1]
  
  ## Resolve name of compartment (metacluster-marker combination)
  if (type=='Differential Abundance') {
    comp <- mc
    marker <- NULL
  } else {
    comp <- paste0(mc, ' ', marker)
  }
  
  ## Resolve which samples to include
  include_samples <- samples_per_test(
    res        = res,
    predictor  = predictor,
    confounder = confounder
  )[, comp, drop = TRUE]
  samples <- names(include_samples)[include_samples]
  
  ## Extract summary stats for the hit
  st <- get_stats(
    res         = res,
    predictor   = predictor,
    confounder  = confounder,
    metacluster = mc,
    marker      = marker
  )
  
  ## Resolve whether confounder specified
  wconf <- !is.null(confounder) && confounder!='none'
  
  ## Resolve interactive/standard point plotting function
  if (interactive) {
    f_point <- ggiraph::geom_point_interactive
  } else {
    f_point <- ggplot2::geom_point
  }
  
  ## Resolve inputs and annotation
  if (type=='Differential Abundance') {
    
    ## Collect abundance data
    y <- data[
      rownames(data)%in%samples,
      gsub('^%', '', colnames(data))==comp,
      drop = TRUE
    ]
    
    ## Resolve axis labels
    lab_pred <- predictor
    lab_resp <- paste0(comp, ' abundance')
  } else if (type=='Differential State (MFI)') {
    
    ## Collect MFI data
    y <- data[
      rownames(data)%in%samples,
      colnames(data)==comp,
      drop = TRUE
    ]
    
    ## Resolve axis labels
    lab_pred <- predictor
    lab_resp <- paste0(comp, ' signal median')
  } else if (type=='Differential State (Phenopositivity)') {
    
    ## Collect phenopositivity data
    y <- data[
      rownames(data)%in%samples,
      colnames(data)==comp,
      drop = TRUE
    ]
    
    ## Resolve axis labels
    lab_pred <- predictor
    lab_resp <- paste0(comp, ' % of phenopositives')
  }
  
  ## Collect predictor, batch and cohort values
  idcs   <- match(samples, annotation$FileName)
  x      <- annotation[idcs, predictor, drop = TRUE]
  batch  <- as.factor(annotation[idcs, 'Batch', drop = TRUE])
  cohort <- as.factor(annotation[idcs, 'Cohort', drop = TRUE])
  
  ## Count included and excluded samples
  n_samples  <- length(include_samples)
  n_included <- sum(include_samples)
  
  ## Collect all filtered inputs in a data frame
  d <- data.frame(
    'Sample'    = samples,
    'Response'  = y,
    'Predictor' = x,
    'Batch'     = batch,
    'Cohort'    = cohort
  )
  
  ## Generate caption with summary stats
  caption <-
    'Excluded samples have missing outcome, predictor or confounder values'
  str_coef <-
    round(
      x      = as.numeric(st$Value[st$Name=='Predictor coefficient']),
      digits = 4
    )
  if (sign(str_coef)==1 && str_coef<0.0001) {
    str_coef <- 'small positive (<0.0001)'
  } else if (sign(str_coef)==(-1) && str_coef>(-0.0001)) {
    str_coef <- 'small negative (>-0.0001)'
  } else {
    str_coef <- as.character(str_coef)
  }
  str_pval <-
    round(
      x = as.numeric(st$Value[st$Name=='Predictor p-value (adjusted)']),
      digits = 4
    )
  if (str_pval<0.0001) {
    str_pval <- as.character('<0.0001')
  } else {
    str_pval <- as.character(str_pval)
  }
  caption <- paste0(
    caption,
    '<br>Predictor coefficient: <b>', str_coef, '</b>, ',
    'predictor <i>p</i>-value (adjusted): <b>', str_pval, '</b>' 
  )
  
  ## Resolve predictor axis limits
  if (is.null(lim_quantiles)) { # quantile limits not specified
    
    if (type=='Differential Abundance') {
      lim_quantiles <- c(.00, .98)
    } else if (type=='Differential State (MFI)') {
      lim_quantiles <- c(.01, .99)
    } else if (type=='Differential State (Phenopositivity)') {
      lim_quantiles <- c(.01, .99)
    }
  }
  caption <- paste0(
    caption, '<br><i>y</i>-axis quantile limits: ',
    paste0(
      '<b>[', paste(lim_quantiles, collapse = ', '), ']</b>'
    )
  )
  
  if (cont) { # continuous predictor
    
    ## Generate subtitle with sample counts
    subtitle <- paste0(
      '<b>N<sub>total</sub>=', n_samples,
      '</b>, N<sub>used</sub>=', n_included
    )
    
    ## Add robustness rate to subtitle if applicable
    if (!is.null(robustness_rates) && !is.null(marker)) {
      
      r <- robustness_rates$Rate[
        robustness_rates$Metacluster==mc & robustness_rates$Marker==marker
      ]
      subtitle <- paste0(
        subtitle, ', robustness rate=',
        round(x = r*100, digits = 2),
        '%'
      )
    }
    
    ## Initialise plot of feature values against a continuous predictor
    p <-
      ggplot2::ggplot(
        data    = d,
        mapping = ggplot2::aes(
          x = .data$Predictor,
          y = .data$Response
        )
      ) +
      ggplot2::xlab(
        lab_pred
      ) +
      ggplot2::ylab(
        lab_resp
      ) +
      ggplot2::theme_minimal(
      ) +
      ggplot2::theme(
        plot.caption.position = 'panel',
        plot.caption          = ggtext::element_markdown(hjust = 0.)
      ) +
      ggplot2::labs(
        title    = paste0('Sample-level data for ', comp),
        subtitle = subtitle,
        caption  = caption
      )
    
    ## Add points/tiles/contours for feature values to the plot
    if (view=='Binned') { # tile per bin
      
      p <- p +
        ggnewscale::new_scale_color(
        ) +
        ggplot2::geom_bin_2d(
          bins = n_bins
        )
      
    } else if (view=='Smooth') { # filled contourss
      
      p <- p +
        ggnewscale::new_scale_color(
        ) +
        ggplot2::geom_density2d_filled(
          bins = n_bins
        ) +
        ggplot2::theme(
          legend.position = 'none'
        )
      
    } else { # point per sample
      
      if (view=='By batch') {
        strat <- 'Batch'
      } else {
        strat <- 'Cohort' 
      }
      
      p <- p +
        ggnewscale::new_scale_color(
        ) +
        ggplot2::geom_point(
          size    = .0,
          alpha   = .0,
          mapping = ggplot2::aes(col = .data[[strat]])
        ) +
        f_point(
          size = .5,
          mapping = ggplot2::aes(
            col     = .data[[strat]],
            tooltip = paste0(
              .data$Sample,
              ' (batch ', .data$Batch, ')'
            )
          )
        )
      
      ## If not interactive, add marginal density plots
      if (!interactive) {
        p <- ggExtra::ggMarginal(
          p,
          type        = 'density',
          size        = 5,
          groupColour = TRUE,
          groupFill   = FALSE
        )
      }
    }
    
    if (type=='Differential State (MFI)') {
      
      ## Specify y-axis scale and limits
      p <- p +
        ggplot2::scale_y_continuous(
          limits = stats::quantile(
            x     = d$Response,
            probs = lim_quantiles
          )
        )
      
      ## Plot noise cuf-offs in terms of MFI if requested
      if (
        show_noise_cutoffs&&!is.null(noise)
      ) {
        
        ## Extract marker-specific cut-offs
        dn <- noise[noise$Marker==marker, , drop = FALSE]
        
        ## Generate labels per cut-off
        dn$Label <- paste0(
          marker, ' Batch ', dn$Batch, ' background noise cut-off: ',
          round(x = dn$Cutoff, digits = 2)
        )
        
        ## Draw noise cut-offs as lines per batch
        if (nrow(dn)>0) {
          
          if (interactive) {
            fline <- ggiraph::geom_hline_interactive
          } else {
            fline <- ggplot2::geom_hline
          }
          p <- p +
            ggnewscale::new_scale_color() +
            fline(
              data = dn,
              mapping = ggplot2::aes(
                yintercept = .data$Cutoff,
                col        = .data$Batch,
                tooltip    = .data$Label
              ),
              show.legend = FALSE
            )
        }
      }
    } else { # DA or DS-Pheno results
      
      ## Specify y-axis scale and limits
      p <- p +
        ggplot2::scale_y_continuous(
          labels = scales::percent,
          limits = stats::quantile(d$Response, lim_quantiles)
        )
    }
  } else { # binary predictor
    
    ## Resolve groups and number of samples per group
    groups    <- levels(as.factor(annotation[, predictor]))
    neg_group <- groups[1]
    pos_group <- groups[2]
    neg_n     <- sum(x==neg_group)
    pos_n     <- sum(x==pos_group)
    
    ## Generate subtitle with excluded/included sample count and positive group
    subtitle <- paste0(
      '<b>N<sub>total</sub>=', n_samples, '</b>, ',
      'N<sub>used</sub>=', n_included, ', ',
      'positive group: <b>', pos_group, '</b>'
    )
    
    ## Add robustness rate to subtitle if applicable
    if (!is.null(robustness_rates) && !is.null(marker)) {
      
      r <- robustness_rates$Rate[
        robustness_rates$Metacluster==mc & robustness_rates$Marker==marker
      ]
      subtitle <- paste0(
        subtitle, ', robustness rate=',
        round(x = r*100, digits = 2),
        '%'
      )
    }
    
    ## Add sample counts by group to subtitle
    subtitle <- paste0(
      subtitle,
      '<br>',
      'N<sub>pos</sub>=', pos_n, ', ',
      'N<sub>neg</sub>=', neg_n
    )
    
    ## Initialise plot of feature values against a binary predictor
    p <-
      ggplot2::ggplot(
        data = d,
        mapping = ggplot2::aes(
          x       = .data$Response,
          y       = .data$Predictor
        )
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::percent,
        limits = stats::quantile(d$Response, lim_quantiles)
      ) +
      ggplot2::xlab(
        lab_resp
      ) +
      ggplot2::ylab(
        lab_pred
      ) +
      ggplot2::theme_minimal(
      ) +
      ggplot2::theme(
        plot.caption.position = 'panel',
        plot.caption          = ggtext::element_markdown(hjust = 0.)
      ) +
      ggplot2::labs(
        title    = paste0('Sample-level data for ', comp),
        subtitle = subtitle,
        caption  = caption
      )
    
    ## Add points/tiles/contours for feature values to the plot
    if (view=='Smooth') { # filled contours
      
      p <- p +
        ggnewscale::new_scale_color(
        ) +
        ggplot2::geom_violin(
          scale = 'count',
          mapping = ggplot2::aes(
            fill = .data$Predictor
          )
        ) +
        ggplot2::theme(
          legend.position = 'none'
        )
      
    } else if (view=='Binned') { # tile per bin
      
      p <- p +
        ggnewscale::new_scale_color(
        ) +
        ggplot2::geom_bin_2d(
        )
      
    } else { # point per sample
      
      if (view=='By batch') {
        strat <- 'Batch'
      } else {
        strat <- 'Cohort' 
      }
      
      set.seed(1) # seed for reproducibility of jittered point positions
      p <- p +
        ggnewscale::new_scale_color(
        ) +
        ggridges::geom_density_ridges(
          panel_scaling = FALSE,
          fill          = scales::alpha('white', alpha = .0),
          alpha         = .3,
          scale         = .5,
          mapping = ggplot2::aes(
            col = .data[[strat]]
          )
        ) +
        f_point(
          mapping = ggplot2::aes(
            col     = .data[[strat]],
            tooltip = paste0(.data$Sample, ' (', .data$Batch, ')')
          ),
          size  = .55,
          alpha = .6,
          position = ggplot2::position_nudge(
            y = sample(
              x       = seq(-0.3, -0.01, length.out = 1000),
              size    = nrow(d),
              replace = TRUE
            )
          )
        )
    }
    
    ## Plot noise cuf-offs in terms of MFI if requested
    if (
      type=='Differential State (MFI)'&&show_noise_cutoffs&&!is.null(noise)
    ) {
      
      ## Extract marker-specific cut-offs
      dn <- noise[noise$Marker==marker, , drop = FALSE]
      
      ## Generate labels per cut-off
      dn$Label <- paste0(
        marker, ' Batch ', dn$Batch, ' background noise cut-off: ',
        round(x = dn$Cutoff, digits = 2)
      )
      
      ## Draw noise cut-offs as lines per batch
      if (nrow(dn)>0) {
        
        if (interactive) {
          fline <- ggiraph::geom_vline_interactive
        } else {
          fline <- ggplot2::geom_vline
        }
        
        p <- p +
          ggnewscale::new_scale_color(
          ) +
          fline(
            data = dn,
            mapping = ggplot2::aes(
              xintercept = .data$Cutoff,
              col        = .data$Batch,
              tooltip    = .data$Label
            ),
            show.legend = FALSE
          )
      }
    }
    
    ## Specify x-axis scale and limits
    if (type=='Differential State (MFI)') {
      
      p <- p +
        ggplot2::scale_x_continuous(
          limits = stats::quantile(
            x     = d$Response, 
            probs = lim_quantiles
          )
        )
    } else {
      p <- p +
        ggplot2::scale_x_continuous(
          limits = stats::quantile(
            x     = d$Response,
            probs = lim_quantiles
          ),
          labels = scales::percent
        )
    }
  }
  
  ## Configure font sizes and Markdown use
  p <- p +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 12),
      plot.subtitle = ggtext::element_markdown(size = 10),
      axis.title.x  = ggplot2::element_text(size = 10),
      axis.title.y  = ggplot2::element_text(size = 10)
    )
  
  p
}

## Function: plot model batch and confounder stats ----

plot_batch_stats <- function(
    res,              # DA, DS-MFI or DS-Pheno results
    mc,               # metacluster of interest
    annotation,       # sample-level annotation
    predictor,        # biological predictor of interest
    marker = NULL,    # state marker of interest (for DS)
    confounder = NULL # biological confounder of interest
) {
  
  ## Resolve analysis type
  type       <- attributes(res)$AnalysisType
  
  ## Resolve base feature-extraction model (left in for future compatibility)
  base_model <- attributes(res)$BaseModel
  
  ## Resolve whether confounder specified
  wconf      <- !is.null(confounder)&&confounder!='none'
  
  ## Resolve whether predictor and confounder is continuous/binary
  cont      <- attributes(res)$Continuous[predictor]
  conf_cont <- is.numeric(annotation[, confounder]) && length(unique(annotation[, confounder]))>2
  
  ## Resolve name of compartment (metacluster-marker combination)
  if (type=='Differential Abundance') {
    comp <- mc
  } else {
    comp <- paste0(mc, ' ', marker)
  }
  
  ## Extract fitted model parameters
  if (type=='Differential Abundance') {
    
    if (wconf) { # confounder specified
      
      ## Extract joint-model and batch-wise model stats
      x_joint <- res$confounders_joint[[predictor]][[confounder]]
      x_batch <- lapply(
        res$confounders_batch,
        function(x) x[[predictor]][[confounder]]
      )
      
      ## Create mask for selected compartment
      idx_comp <- x_joint$Compartment==comp
      
      ## Get predictor coefficient & p-val in joint model
      fc_joint <- x_joint$FC[idx_comp]
      ch_joint <- sign(fc_joint)*(abs(fc_joint)-1.)
      p_joint  <- x_joint$AdjPVal[idx_comp] 
      s_joint  <- sig_asterisks(p_joint)
      
      ## Get predictor coefficients & p-vals in batch models
      fc_batch <- sapply(x_batch, function(x) x$FC[idx_comp])
      ch_batch <- sapply(fc_batch, function(x) sign(x)*(abs(x)-1.))
      p_batch  <- sapply(x_batch, function(x) x$AdjPVal[idx_comp])
      s_batch  <- sig_asterisks(p_batch)
      
      ## Get confounder coefficient & p-val in joint model
      c_fc_joint <- x_joint$FCConfounder[idx_comp]
      c_ch_joint <- sign(c_fc_joint)*(abs(c_fc_joint)-1.)
      c_p_joint  <- x_joint$AdjPValConfounder[idx_comp] 
      c_s_joint  <- sig_asterisks(c_p_joint)
      
      ## Get confounder coefficients & p-vals in batch models
      c_fc_batch <- sapply(x_batch, function(x) x$FCConfounder[idx_comp])
      c_ch_batch <- sapply(c_fc_batch, function(x) sign(x)*(abs(x)-1.))
      c_p_batch  <- sapply(x_batch, function(x) x$AdjPValConfounder[idx_comp])
      c_s_batch  <- sig_asterisks(c_p_batch)
      
      ## Determine positions of labels in the plot
      bp        <- c(ch_batch, ch_joint)
      lp        <- bp+sign(bp)*diff(range(bp))/30
      lp[lp>=0] <- lp[lp>=0]-diff(range(bp))/50
      lp[lp<0]  <- lp[lp<0]-diff(range(bp))/20
      c_bp      <- c(c_ch_batch, c_ch_joint)
      c_lp      <- c_bp-diff(range(c_bp))/20+sign(c_bp)*diff(range(c_bp))/30
      
      ## Gather predictor stats in a data frame
      pred_title <- paste0('Predictor: ', predictor)
      d_pred <- data.frame(
        'Batch'         = as.factor(
          sub('Batch', '', c(names(ch_batch), 'ALL'))
        ),
        'LogFC'         = c(fc_batch, fc_joint),
        'Change'        = c(ch_batch, ch_joint),
        'AdjPVal'       = c(p_batch, p_joint),
        'Significance'  = c(s_batch, s_joint),
        'LabelPosition' = lp,
        'IsJoint'       = c(rep(FALSE, times = length(ch_batch)), TRUE),
        'Effect'        = pred_title
      )
      
      ## Gather confounder stats in a data frame
      conf_title <- paste0('Confounder: ', confounder)
      d_conf <- data.frame(
        'Batch'         = as.factor(
          sub('Batch', '', c(names(c_ch_batch), 'ALL'))
        ),
        'LogFC'         = c(c_fc_batch, c_fc_joint),
        'Change'        = c(c_ch_batch, c_ch_joint),
        'AdjPVal'       = c(c_p_batch, c_p_joint),
        'Significance'  = c(c_s_batch, c_s_joint),
        'LabelPosition' = c_lp,
        'IsJoint'       = c(rep(FALSE, times = length(c_ch_batch)), TRUE),
        'Effect'        = conf_title
      )
      
      ## Make sure batch flags remain ordered in the plot
      b  <- unique(d_pred$Batch)
      bn <- as.integer(gsub('ALL', '999', gsub('Batch', '', b)))
      for (bl in b[order(bn, decreasing = TRUE)]) {
        
        d_pred$Batch <- relevel(d_pred$Batch, ref = bl)
        d_conf$Batch <- relevel(d_conf$Batch, ref = bl)
      }
      
      ## Combine predictor and confounder data frame
      d <- rbind(d_pred, d_conf)
      
      ## Make sure predictor comes before confounder in the plot
      d$Effect <- as.factor(d$Effect)
      d$Effect <- relevel(d$Effect, ref = pred_title)
      
      ## Generate plot of effect size and significance by batch
      p <-
        ggplot2::ggplot(
          data = d,
          mapping = ggplot2::aes(
            x    = .data$Batch,
            y    = .data$Change,
            col  = .data$Batch,
            fill = .data$IsJoint
          )
        ) +
        ggplot2::facet_wrap(
          ~.data$Effect,
          nrow   = 2,
          scales = 'free_y'
        ) +
        ggplot2::xlab(
          ''
        ) +
        ggplot2::ylab(
          'Slope coefficient'
        ) +
        ggplot2::geom_bar( # bars for effect sizes
          stat       = 'identity',
          width      = .6,
          linewidth  = .5
        ) +
        ggplot2::geom_hline( # zero-effect line
          yintercept = 0.,
          size       = .5
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::geom_text( # asterisk labels per significance
          size     = 7,
          fontface = 'bold',
          vjust    = .5,
          mapping = ggplot2::aes(
            x     = .data$Batch,
            y     = .data$LabelPosition,
            label = .data$Significance
          )
        ) +
        ggplot2::scale_fill_manual(
          values = c('white', 'lightgrey')
        ) +
        ggplot2::theme(
          legend.position = 'none',
          axis.text.x     = ggplot2::element_text(
            angle = 0,vjust = 0., hjust = .5, size = 10
          ),
          axis.text.y     = ggplot2::element_text(size = 10),
          plot.title      = ggplot2::element_text(size = 12),
          strip.text      = ggplot2::element_text(size = 10)
        ) +
        ggplot2::ggtitle(
          'edgeR model coefficients by batch and jointly'
        )
    } else { # no confounder specified
      
      ## Extract joint-model and batch-wise model stats
      x_joint <- res$joint[[predictor]]
      x_batch <- lapply(
        res$batch,
        function(x) x[[predictor]]
      )
      
      ## Create mask for selected compartment
      idx_comp <- x_joint$Compartment==comp
      
      ## Get predictor coefficient & p-val in joint model
      fc_joint <- x_joint$FC[idx_comp]
      ch_joint <- sign(fc_joint)*(abs(fc_joint)-1.)
      p_joint  <- x_joint$AdjPVal[idx_comp] 
      s_joint  <- sig_asterisks(p_joint)
      
      ## Get predictor coefficients & p-vals in batch models
      fc_batch <- sapply(x_batch, function(x) x$FC[idx_comp])
      ch_batch <- sapply(fc_batch, function(x) sign(x)*(abs(x)-1.))
      p_batch  <- sapply(x_batch, function(x) x$AdjPVal[idx_comp])
      s_batch  <- sig_asterisks(p_batch)
      
      ## Determine positions of labels in the plot
      bp        <- c(ch_batch, ch_joint)
      lp        <- bp+sign(bp)*diff(range(bp))/30
      lp[lp>=0] <- lp[lp>=0]-diff(range(bp))/50
      lp[lp<0]  <- lp[lp<0]-diff(range(bp))/20
      
      ## Gather predictor stats in a data frame
      pred_title <- paste0('Predictor: ', predictor)
      d <- data.frame(
        'Batch'         = as.factor(c(paste0(names(ch_batch)), 'ALL')),
        'LogFC'         = c(fc_batch, fc_joint),
        'Change'        = c(ch_batch, ch_joint),
        'AdjPVal'       = c(p_batch, p_joint),
        'Significance'  = c(s_batch, s_joint),
        'LabelPosition' = lp,
        'IsJoint'       = c(rep(FALSE, times = length(ch_batch)), TRUE),
        'Effect'        = as.factor(pred_title)
      )
      
      ## Make sure batch flags remain ordered in the plot
      b <- unique(d$Batch)
      bn <- as.integer(gsub('ALL', '999', gsub('Batch', '', b)))
      for (bl in b[order(bn, decreasing = TRUE)]) {
        d$Batch <- relevel(d$Batch, ref = bl)
      }
      
      ## Generate plot of effect size and significance by batch
      p <-
        ggplot2::ggplot(
          data = d,
          mapping = ggplot2::aes(
            x    = .data$Batch,
            y    = .data$Change,
            col  = .data$Batch,
            fill = .data$IsJoint
          )
        ) +
        ggplot2::facet_wrap(
          ~.data$Effect
        ) +
        ggplot2::xlab(
          ''
        ) +
        ggplot2::ylab(
          'Change per unit increment'
        ) +
        ggplot2::geom_bar( # bars for effect sizes
          stat       = 'identity',
          width      = .6,
          linewidth  = .5
        ) +
        ggplot2::geom_hline(
          yintercept = 0.,
          size       = .5
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::geom_text( # asterisk labels per significance
          size     = 7,
          fontface = 'bold',
          vjust    = .5,
          mapping = ggplot2::aes(
            x     = .data$Batch,
            y     = .data$LabelPosition,
            label = .data$Significance
          )
        ) +
        ggplot2::scale_fill_manual(
          values = c('white', 'lightgrey')
        ) +
        ggplot2::theme(
          legend.position = 'none',
          axis.text.x     = ggplot2::element_text(
            angle = 0, vjust = 0., hjust = .5, size = 10
          ),
          axis.text.y     = ggplot2::element_text(size = 10),
          plot.title      = ggplot2::element_text(size = 12),
          strip.text      = ggplot2::element_text(size = 10)) +
        ggplot2::ggtitle(
          'edgeR model coefficients by batch and jointly'
        )
    }
    
  } else { # DS-MFI or DS-Pheno results
    
    ## Resolve outcome and axis labels
    if (type=='Differential State (MFI)') {
      
      outcome       <- 'mfi'
      lab_yaxis_ri  <- paste0(comp, ' signal median')
      lab_yaxis_rsq <- 'batch-level R<sup>2</sup> estimate'
    } else if (type=='Differential State (Phenopositivity)') {
      
      outcome       <- 'phenopos'
      lab_yaxis_ri  <- paste0(comp, ' % of phenopositives')
      lab_yaxis_rsq <- 'batch-level R<sup>2</sup> estimate'
    }
    
    if (wconf) { # confounder specified
      
      ## Extract predictor and confounder effect coefficients
      if (outcome=='mfi') {
        ch_pred <- res$confounders_main[[predictor]][[confounder]]$change[
          res$confounders_main[[predictor]][[confounder]]$Compartment==comp
        ]
        ch_conf <- res$confounders_main[[predictor]][[confounder]]$changeConfounder[
          res$confounders_main[[predictor]][[confounder]]$Compartment==comp
        ]
      } else if (outcome=='phenopos') {
        ch_pred <- res$confounders_main[[predictor]][[confounder]]$odds[
          res$main[[predictor]]$Compartment==comp
        ]-1.
        ch_conf <- res$confounders_main[[predictor]][[confounder]]$oddsConfounder[
          res$confounders_main[[predictor]][[confounder]]$Compartment==comp
        ]
      }
      
      ## Extract predictor and confounder p-values
      p_pred <- res$confounders_main[[predictor]][[confounder]]$AdjPVal[
        res$confounders_main[[predictor]][[confounder]]$Compartment==comp
      ]
      s_pred <- sig_asterisks(p_pred)
      p_conf <- res$confounders_main[[predictor]][[confounder]]$AdjPValConfounder[
        res$confounders_main[[predictor]][[confounder]]$Compartment==comp
      ]
      s_conf <- sig_asterisks(p_conf)
      
      ## Extract batch-level random intercept values and confidence intervals
      x_ri <-
        res$confounders_random_intercepts[[predictor]][[confounder]]
      ri_avg_pred <-
        x_ri$Intercept[x_ri$Compartment==comp]
      ri_min_pred <-
        x_ri$ConfidenceIntervalMin[x_ri$Compartment==comp]
      ri_max_pred <-
        x_ri$ConfidenceIntervalMax[x_ri$Compartment==comp]
      batches <-
        x_ri$Batch[x_ri$Compartment==comp]
      
      d_ri <- data.frame(
        'Label'   = factor(batches, levels = batches),
        'Batch'   = factor(batches, levels = batches),
        'RI'      = ri_avg_pred,
        'CIMin'   = ri_min_pred,
        'CIMax'   = ri_max_pred
      )
      
      ## Extract batch-level R-squared estimates
      x_rsq <-
        res$confounders_main[[predictor]][[confounder]][, c('Compartment', 'Rsq')]
      d_rsq <- data.frame(
        'Batch' = as.factor('Joint'),
        'Rsq'   = x_rsq[, 'Rsq'][x_rsq$Compartment==comp]
      )
      x_brsq <-
        res$confounders_batch_r_squared[[predictor]][[confounder]]
      d_brsq <- data.frame(
        'Label'    = factor(batches, levels = batches),
        'Batch'    = factor(batches, levels = batches),
        'BatchRsq' = x_brsq$Rsq[x_brsq$Compartment==comp]
      )
      
      ## Generate plot of random intercept ranges per batch
      p_ri <-
        ggplot2::ggplot(
          data = d_ri,
          mapping = ggplot2::aes(
            x   = .data$Label,
            col = .data$Batch
          )
        ) +
        ggplot2::xlab(
          'Batch'
        ) +
        ggplot2::ylab(
          lab_yaxis_ri
        ) +
        ggplot2::geom_point(
          mapping = ggplot2::aes(y = .data$RI)
        ) +
        ggplot2::geom_linerange(
          mapping = ggplot2::aes(
            ymin = .data$CIMin,
            ymax = .data$CIMax
          )
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::theme(
          legend.position = 'none',
          axis.text.x = ggplot2::element_text(
            angle = 0, vjust = 0.5, hjust = .5, size = 10
          ),
          axis.title.y = ggplot2::element_text(size = 10)
        ) +
        ggplot2::ggtitle(
          'Batch-level intercepts'
        )
      
      ## Generage plot of R-squared estimates per batch
      p_brsq <- ggplot2::ggplot(
        data    = d_brsq,
        mapping = ggplot2::aes(
          x   = .data$Label,
          y   = .data$BatchRsq,
          col = .data$Batch
        )
      ) +
        ggplot2::geom_hline( # overall R^2 value line
          yintercept = d_rsq$Rsq,
          lwd        = 0.6,
          linetype   = 2,
          col        = 'darkgrey'
        ) +
        ggplot2::geom_hline( # R^2~0 line
          yintercept = 0,
          lwd        = 0.6,
          linetype   = 1,
          col        = '#2b2b2b'
        ) +
        ggplot2::geom_segment(
          mapping = ggplot2::aes(
            x    = .data$Label,
            xend = .data$Label,
            y    = .data$BatchRsq,
            yend = 0.
          )
        ) +
        ggplot2::geom_point(
          mapping = ggplot2::aes(y = .data$BatchRsq)
        ) +
        ggplot2::xlab(
          'Batch'
        ) + 
        ggplot2::ylab(
          lab_yaxis_rsq
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::theme(
          legend.position = 'none',
          axis.text.x = ggplot2::element_text(
            angle = 0, vjust = 0.5, hjust = .5, size = 10
          ),
          axis.title.y = ggtext::element_markdown(size = 10)
        ) +
        ggplot2::ggtitle(
          'Batch-level goodness-of-fit'
        )
      
      ## Resolve predictor and confounder titles
      if (cont) { # continuous predictor
        pred_title <- predictor
      } else { # binary predictor
        pred_title <- paste0( # specify positive group
          predictor, ' (', levels(as.factor(annotation[, predictor]))[2], ')'
        )
      }
      if (conf_cont) { # continuous confounder
        conf_title <- confounder
      } else { # binary confounder
        conf_title <- paste0( # specify positive group
          confounder, ' (', levels(as.factor(annotation[, confounder]))[2], ')'
        )
      }
      
      ## Collect predictor and confounder effect and significant stats
      cp <- c(ch_pred, ch_conf)
      d_eff <- data.frame(
        'Effect'    = as.factor(c(pred_title, conf_title)),
        'LabelType' = as.factor(c(
          'adjusted <i>p</i>-value',
          'adjusted <i>p</i>-value',
          'effect direction',
          'effect direction'
        )),
        'Label'     = c(
          sig_level(c(p_pred, p_conf), prefix = FALSE),
          ifelse(c(ch_pred, ch_conf)>0, '↑', '↓')
        ),
        'Asterisks' = sig_asterisks(c(p_pred, p_conf)),
        'Signif'    = c(p_pred, p_conf)<.05
      )
      
      ## Generate a grid of predictor and confounder effect and significance
      p_eff <-
        ggplot2::ggplot(
          data    = d_eff,
          mapping = ggplot2::aes(
            x   = .data$Effect,
            y   = .data$LabelType,
            col = .data$Signif)
        ) +
        ggplot2::scale_size_manual(
          values = c(5, 7.5)
        ) +
        ggplot2::scale_x_discrete(
          position = 'top'
        ) +
        ggplot2::geom_text(
          mapping = ggplot2::aes(
            y     = .data$LabelType,
            x     = .data$Effect,
            label = .data$Label,
            size  = .data$LabelType
          ),
          color = dplyr::if_else(
            d_eff$Signif, 'black', '#858585'
          )
        ) +
        ggplot2::geom_vline(
          xintercept = 1.5,
          color      = 'darkgrey'
        ) +
        ggplot2::geom_hline(
          yintercept = 1.5,
          color      = 'darkgrey'
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::ggtitle(
          'Effects'
        ) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect('#f5f5f5', 'white', 5),
          legend.position  = 'none',
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.spacing    = ggplot2::unit(0, 'cm'),
          axis.ticks       = ggplot2::element_blank(),
          axis.text.y      = ggtext::element_markdown(
            size = 10, colour = 'black', angle = 90, hjust = .5
          ),
          axis.text.x      = ggtext::element_markdown(
            size = 10, colour = 'black'
          ),
          axis.title.x     = ggplot2::element_blank(),
          axis.title.y     = ggplot2::element_blank()
        )
      
      ## Combine RI, goodness-of-fit and predictor-vs-confounder plots
      p <- cowplot::plot_grid(
        cowplot::plot_grid(
          p_ri,
          p_brsq,
          nrow = 2
        ),
        p_eff,
        ncol = 2,
        rel_widths = c(length(batches)/2, 3)
      )
    } else { # no confounder specified
      
      ## Extract predictor effect coefficients
      if (outcome=='mfi') {
        ch_pred <- res$main[[predictor]]$change[
          res$main[[predictor]]$Compartment==comp
        ]
      } else if (outcome=='phenopos') {
        ch_pred <- res$main[[predictor]]$odds[
          res$main[[predictor]]$Compartment==comp
        ]-1.
      }
      
      ## Extract predictor p-values
      p_pred <- res$main[[predictor]]$AdjPVal[
        res$main[[predictor]]$Compartment==comp
      ]
      s_pred <- sig_asterisks(p_pred)
      
      ## Extract batch-level random intercept values and confidence intervals
      x_ri        <- res$random_intercepts[[predictor]]
      ri_avg_pred <- x_ri$Intercept[x_ri$Compartment==comp]
      ri_min_pred <- x_ri$ConfidenceIntervalMin[x_ri$Compartment==comp]
      ri_max_pred <- x_ri$ConfidenceIntervalMax[x_ri$Compartment==comp]
      batches     <- x_ri$Batch[x_ri$Compartment==comp]
      d_ri <- data.frame(
        'Label'   = factor(batches, levels = batches),
        'Batch'   = factor(batches, levels = batches),
        'RI'      = ri_avg_pred,
        'CIMin'   = ri_min_pred,
        'CIMax'   = ri_max_pred
      )
      
      ## Extract batch-level R-squared estimates
      x_rsq <- res$main[[predictor]][, c('Compartment', 'Rsq')]
      d_rsq <- data.frame(
        'Batch' = as.factor('Joint'),
        'Rsq'   = x_rsq[, 'Rsq'][x_rsq$Compartment==comp]
      )
      x_brsq <- res$batch_r_squared[[predictor]]
      d_brsq <- data.frame(
        'Label'    = factor(batches, levels = batches),
        'Batch'    = factor(batches, levels = batches),
        'BatchRsq' = x_brsq$Rsq[x_brsq$Compartment==comp]
      )
      
      ## Generate plot of random intercept ranges per batch
      p_ri <- ggplot2::ggplot(
        data    = d_ri,
        mapping = ggplot2::aes(
          x   = .data$Label,
          col = .data$Batch
        )
      ) +
        ggplot2::xlab(
          ''
        ) +
        ggplot2::ylab(
          lab_yaxis_ri
        ) +
        ggplot2::geom_point(
          mapping = ggplot2::aes(
            y = .data$RI
          )
        ) +
        ggplot2::geom_linerange(
          mapping = ggplot2::aes(
            ymin = .data$CIMin,
            ymax = .data$CIMax
          )
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::theme(
          legend.position = 'none',
          axis.text.x = ggplot2::element_text(
            angle = 0, vjust = 0.5, hjust = .5, size = 10
          ),
          axis.title.y = ggplot2::element_text(size = 10)
        ) +
        ggplot2::ggtitle(
          'Batch-level intercepts'
        )
      
      ## Generage plot of R-squared estimates per batch
      p_brsq <-
        ggplot2::ggplot(
          data = d_brsq,
          mapping = ggplot2::aes(
            x   = .data$Label,
            y   = .data$BatchRsq,
            col = .data$Batch
          )
        ) +
        ggplot2::geom_hline( # overall R^2 value line
          yintercept = d_rsq$Rsq,
          lwd        = 0.6,
          linetype   = 2,
          col        = 'darkgrey'
        ) +
        ggplot2::geom_hline( # R^2~0 line
          yintercept = 0,
          lwd        = 0.6,
          linetype   = 1,
          col        = '#2b2b2b'
        ) +
        ggplot2::geom_segment(
          mapping = ggplot2::aes(
            x    = .data$Label,
            xend = .data$Label,
            y    = 0.,
            yend = .data$BatchRsq
          )
        ) +
        ggplot2::geom_point(
          mapping = ggplot2::aes(
            y = .data$BatchRsq
          )
        ) +
        ggplot2::xlab(
          'Batch'
        ) + 
        ggplot2::ylab(
          lab_yaxis_rsq
        ) +
        ggplot2::theme_minimal(
        ) +
        ggplot2::theme(
          legend.position = 'none',
          axis.text.x = ggplot2::element_text(
            angle = 0, vjust = 0.5, hjust = .5, size = 10
          ),
          axis.title.y = ggtext::element_markdown(size = 10)
        ) +
        ggplot2::ggtitle(
          'Batch-level goodness-of-fit'
        )
      
      ## Combine RI and goodness-of-fit plots
      p <- cowplot::plot_grid(p_ri, p_brsq, nrow = 2)
    }
  }
  p
}

## Function: binned biaxial expression plot ----

plot_comp_scatter <- function(
    comp,                        # metacluster of interest
    m1,                          # x-axis marker
    m2,                          # y-axis marker
    channels,                    # channel-marker mapping
    fsom,                        # FlowSOM model
    thresholds        = NULL,    # phenopositivity thresholds per marker
    noise             = NULL,    # noise levels per marker and batch
    view =                 # whether to plot phenopositivity thresholds or...
      c('pheno', 'noise'), # ...noise cut-offs
    agg_unstained     = NULL,    # aggregate unstained expression data
    bg_unstained      = FALSE,   # whether background should be unstained
    bg_only           = FALSE,   # whether to plot background only
    n_bins            = 150,     # number of bins for aggregation
    n_downsample_comp = 1e3,     # max subsampling count for metacluster
    n_downsample_rest = 1e5,     # max subsampling count for other cells
    seed              = 1        # random seed for subsampling
) {
  
  ## Resolve whether to draw phenopositivity thresholds or noise cut-offs
  view <- view[1]
  
  ## Resolve whether to plot background signal
  bg_unstained <- bg_unstained&&!is.null(agg_unstained)
  
  ## Resolve channel names
  ch1 <- names(channels)[channels==m1]
  ch2 <- names(channels)[channels==m2]
  
  ## Resolve metacluster name and index
  str_comp <- comp
  n        <- nrow(fsom$data)
  meta     <- grepl('^MC', comp)
  comp     <- as.integer(gsub('^MC', '', gsub('^C', '', comp)))
  
  ## Get mask for cells in metacluster
  idcs <- which(FlowSOM::GetMetaclusters(fsom)==comp)
  mask <- rep(FALSE, times = n)
  mask[idcs] <- TRUE
  
  ## Determine axis limits from all available data (incl. unstained if given)
  margin <- 1e-5
  xl <- stats::quantile(
    x = if (!is.null(agg_unstained)) {
      c(fsom$data[, ch1], agg_u@exprs[, ch1])
    } else {
      fsom$data[, ch1]
    },
    probs = c(margin, 1.-margin),
    na.rm = TRUE
  )
  yl <- stats::quantile(
    x = if (!is.null(agg_unstained)) {
      c(fsom$data[, ch2], agg_u@exprs[, ch2])
    } else {
      fsom$data[, ch2]
    },
    probs = c(margin, 1.-margin),
    na.rm = TRUE
  )
  
  ## Subsample data from metacluster if needed
  if (!bg_only) {
    
    idcs_comp <- which(mask)
    if (
      length(idcs_comp)>n_downsample_comp
    ) {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      idcs_comp <- idcs_comp[
        sample(seq_along(idcs_comp), size = n_downsample_comp)
      ]
    }
    d_comp <- fsom$data[
      idcs_comp, c(ch1, ch2), drop = FALSE
    ]
    colnames(d_comp) <- c('Channel1', 'Channel2')
  }
  
  ## Subsample data from the rest if needed
  if (bg_unstained) { # unstained signal in the background
    
    bg        <- agg_unstained@exprs[, c(ch1, ch2)]
    idcs_rest <- seq_len(nrow(bg))
    if (
      length(idcs_rest)>n_downsample_rest
    ) {
      
      if (!is.null(seed)) {
        set.seed(seed)
      }
      idcs_rest <- idcs_rest[
        sample(seq_along(idcs_rest), size = n_downsample_rest)
      ]
    }
    d_rest <- bg[idcs_rest, ]
  } else { # stained signal in the background
    
    idcs_rest <- which(!mask)
    if (
      length(idcs_rest)>n_downsample_rest
    ) {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      idcs_rest <- idcs_rest[
        sample(seq_along(idcs_rest), size = n_downsample_rest)
      ]
    }
    d_rest <- fsom$data[ # FlowSOM aggregate data
      idcs_rest, c(ch1, ch2), drop = FALSE
    ]
  }
  colnames(d_rest) <- c('Channel1', 'Channel2')
  
  ## Initialise biaxial expression plot
  p <-
    ggplot2::ggplot(
    ) +
    ggplot2::theme(
      legend.position = 'none'
    )
  
  ## Generate plot background
  p1 <- p +
    ggplot2::geom_bin2d(
      data    = d_rest,
      mapping = ggplot2::aes(
        x = .data$Channel1,
        y = .data$Channel2
      ),
      bins =
        if (bg_unstained) {
          ceiling(n_bins/2)
        } else {
          n_bins
        }
    )
  
  ## Change colour scale if background is unstained
  if (bg_unstained) {
    
    p1 <- p1 +
      ggplot2::scale_fill_gradient(
        low  = 'aquamarine4',
        high = 'black'
      )
  }
  
  ## Generate plot foreground
  if (!bg_only) { 
    
    p2 <- p + 
      ggplot2::geom_bin2d(
        data = d_comp,
        mapping = ggplot2::aes(
          x = .data$Channel1,
          y = .data$Channel2,
        ),
        bins = ceiling(n_bins/3),
        alpha = .8
      ) +
      ggplot2::scale_fill_gradient(
        low  = 'gold',
        high = 'brown'
      )
  }
  foreground <- if (!bg_only) { # metacluster to be highlighted
    ggplot2::geom_tile(
      data = ggplot2::layer_data(p2),
      linejoin = 'round'
    )
  } else { # no metacluster highlighting
    NULL
  }
  
  ## Combine background and foreground into one plot
  pp <- ggplot2::ggplot(
    data    = ggplot2::layer_data(p1),
    mapping = ggplot2::aes(
      x    = .data$x,
      y    = .data$y,
      fill = .data$fill
    )
  ) +
    ggplot2::geom_tile(
      linejoin = 'round'
    ) +
    foreground +
    ggplot2::scale_fill_identity(
    ) +
    ggplot2::xlab(
      m1
    ) +
    ggplot2::ylab(
      m2
    ) +
    ggplot2::theme_grey(
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(c(0.5, 1.5, 0.5, 1.5), 'cm')
    ) +
    ggplot2::ggtitle(
      NULL
    )
  
  if (
    view=='pheno'&&!is.null(thresholds)
  ) { # phenopositivity thresholds to be drawn
    
    if (m1%in%names(thresholds)) {
      
      ## Add x-axis marker phenopositivity threshold plot
      pp <- pp +
        ggiraph::geom_vline_interactive(
          xintercept = thresholds[m1],
          tooltip    = paste0(
            m1, ' phenopositivity threshold: ',
            round(x = thresholds[m1], digits = 2)
          ),
          linewidth  = 2.,
          color      = 'darkorchid1',
          alpha      = .6
        ) +
        ggplot2::theme(
          legend.position = 'none',
          plot.margin     = grid::unit(c(.5,.1,.1,.1), 'cm')
        )
    }
    if (m2%in%names(thresholds)) {
      
      ## Add y-axis marker phenopositivity threshold to plot
      pp <- pp +
        ggiraph::geom_hline_interactive(
          yintercept = thresholds[m2],
          tooltip    = paste0(
            m2, ' phenopositivity threshold: ',
            round(x = thresholds[m2], digits = 2)
          ),
          linewidth  = 2.,
          color      = 'darkorchid1',
          alpha      = .6
        ) +
        ggplot2::theme(
          legend.position = 'none',
          plot.margin     = grid::unit(c(.5,.1,.1,.1), 'cm')
        )
    }
  } else if (
    view=='noise'&&!is.null(noise)
  ) { # noise cut-offs to be drawn
    
    if (m1%in%noise$Marker) {
      
      ## Extract and label x-axis marker noise cut-offs per batch
      d1 <- noise[noise$Marker==m1, , drop = FALSE]
      d1$Label <- paste0(
        m1, ' batch ', d1$Batch, ' cut-off: ',
        round(x = d1$Cutoff, digits = 2)
      )
      
      ## Add cut-offs to plot
      pp <- pp +
        ggiraph::geom_vline_interactive(
          data = d1,
          mapping = ggplot2::aes(
            xintercept = .data$Cutoff,
            col        = .data$Batch,
            tooltip    = .data$Label
          )
        ) +
        ggplot2::theme(
          legend.position = 'none',
          plot.margin     = grid::unit(c(.5,.1,.1,.1), 'cm')
        )
    }
    if (m2%in%noise$Marker) {
      
      ## Extract and label y-axis marker noise cut-offs per batch
      d2 <- noise[noise$Marker==m2, , drop = FALSE]
      d2$Label <- paste0(
        m2, ' batch ', d2$Batch, ' cut-off: ',
        round(x = d2$Cutoff, digits = 2)
      )
      
      ## Add cut-offs to plot
      pp <- pp +
        ggiraph::geom_hline_interactive(
          data = d2,
          mapping = ggplot2::aes(
            yintercept = .data$Cutoff,
            col        = .data$Batch,
            tooltip    = .data$Label
          )
        ) +
        ggplot2::theme(
          legend.position  = 'none',
          plot.margin      = grid::unit(c(.5,.1,.1,.1), 'cm')#,
          # panel.background = ggplot2::element_rect(fill = '#f2f2f2')
        )
    }
  }
  
  ## Specify axis limits
  pp <- pp + ggplot2::xlim(xl) + ggplot2::ylim(yl)
  
  pp
}

## Function: plot metacluster-population composition ----

plot_comp_profile <- function(
    prof,   # metacluster profiles
    counts, # metacluster-sample cell counts matrix
    comp,   # metacluster of interest
    n       # cell counts per metacluster
) {
  
  ## Return empty plot if no metacluster profiles given
  if (is.null(prof)) {
    
    return(
      ggplot2::ggplot() +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(
            fill = '#f7f7f7', colour = '#f7f7f7'
          ),
          plot.background = ggplot2::element_rect(
            fill = '#f7f7f7', colour = '#f7f7f7'
          )
        )
    )
  }
  
  ## Resolve base feature-extraction model (left in for future compatibility)
  base_model <- attributes(prof)$BaseModel
  
  ## Get metacluster and total cell counts
  this_n <- n[comp]
  all_n  <- sum(n)
  
  ## Get metacluster title and index  
  tidy_comp <- gsub('^MC', 'Metacluster ', comp)
  comp      <-  gsub('^MC', '', comp)
  
  ## Gather cell counts per population in metacluster
  pops   <- as.factor(colnames(prof))
  cc     <- prof[comp, ]
  d <- data.frame(
    'Population'    = pops,
    'Count'         = cc,
    'Proportion'    = cc/sum(cc),
    'LabelPosition' = NA
  )
  rownames(d) <- NULL
  
  ## Determine content and positions of text labels per population
  d$csum <- rev(cumsum(rev(d$Proportion)))
  pos    <- d$Proportion/2+lead(d$csum, 1)
  pos    <- dplyr::if_else(is.na(pos), d$Proportion/2, pos)
  pos[d$Proportion < 0.01] <- NA
  d$LabelPosition <- pos
  
  ## Get metacluster size as proportion in the aggregate expression data
  comp_size <- this_n/all_n
  
  ## Generate stacked bar plot of metacluster-population composition
  p <-
    ggplot2::ggplot(
      data = d,
      mapping = ggplot2::aes(
        x    = '',
        y    = .data$Proportion,
        fill = .data$Population
      )
    ) +
    ggplot2::scale_fill_manual(
      values = rep(c('#f7feff', '#edfeff'), times = 100)
    ) +
    ggplot2::geom_bar(
      stat = 'identity'
    ) +
    ggplot2::geom_col(
      width = 1, color = 1
    ) +
    ggplot2::theme_void(
    ) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = dplyr::if_else(
          condition = .data$Proportion>=.05,
          true      = paste0(.data$Population, ': ', round(.data$Proportion*100, digits = 1), '%'),
          false     = ''
        )
      ),
      size = 4,
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::ggtitle(
      paste0(
        '<b>', tidy_comp, ' composition</b><br><i>(',
        round(x = comp_size*100, digits = 2),
        '% of cells in aggregate)</i>'
      )
    ) +
    ggplot2::theme(
      legend.position  = 'none',
      plot.title       = ggtext::element_markdown(hjust = 0.5, size = 14),
      panel.background = ggplot2::element_rect(
        fill = '#f7f7f7', colour = '#f7f7f7'
      ),
      plot.background = ggplot2::element_rect(
        fill = '#f7f7f7', colour = '#f7f7f7'
      )
    )
  
  p
}

## Function: plot metacluster signal quarties ----

plot_comp_quartiles <- function(
    comp,                  # metacluster of interest
    quartiles,             # metacluster signal quartiles
    lineage_markers,       # names of analysed lineage markers
    state_markers,         # names of analysed state markers
    noise = NULL,          # noise levels per marker and batch
    thresholds = NULL,     # phenopositivity thresholds per marker
    show_unstained = FALSE # whether to show unstained signal as well
) {
  
  ## Gather signal quartiles per marker for metacluster of interest
  d <- quartiles[
    quartiles$Compartment==comp|
      quartiles$Compartment=='All'|
      quartiles$Compartment=='Unstained', ,
    drop = FALSE
  ]
  
  ## Determine whether unstained signal is available
  unstained_avail <- any(quartiles$Compartment=='Unstained')
  
  ## Label markers by lineage/state
  d$Type <- as.factor(dplyr::if_else(
    d$Marker%in%lineage_markers,
    'Lineage', 'State'
  ))
  
  ## Adapt metacluster signal label for plotting
  lab_comp <- paste0(comp, ' signal')
  levels(d$Compartment)[levels(d$Compartment)==comp] <- lab_comp
  
  if (unstained_avail && show_unstained) { # unstained signal to be plotted
    
    ## Adapt unstained signal labels for plotting
    levels(d$Compartment)[
      levels(d$Compartment)=='All'
    ] <- 'Aggregate signal (stained)'
    levels(d$Compartment)[
      levels(d$Compartment)=='Unstained'
    ] <- 'Aggregate signal (unstained)'
    d$Compartment <-
      relevel(d$Compartment, lab_comp)
    d$Compartment <-
      relevel(d$Compartment, 'Aggregate signal (stained)')
    d$Compartment <-
      relevel(d$Compartment, 'Aggregate signal (unstained)')
    
    ## Specify colour per signal type
    palette <- c('azure4', 'cyan3', 'coral2')
  } else { # unstained signal not to be plotted
    
    ## Exclude unstained signal data
    if (!show_unstained) {
      d <- d[d$Compartment!='Unstained', ]
    }
    
    ## Adapt aggregate signal labels for plotting
    levels(d$Compartment)[levels(d$Compartment)=='All'] <-
      'Aggregate signal'
    d$Compartment <- relevel(d$Compartment, 'Aggregate signal')
    
    ## Specifcy colour per signal type
    palette <- c('cyan3', 'coral2')
  }
  
  ## Compute average noise cut-offs per marker
  if (!is.null(noise)) {
    
    cutoffs <- sapply(
      levels(d$Marker),
      function(m) mean(noise$Cutoff[noise$Marker==m])
    )
    d$Cutoff <- cutoffs[d$Marker]
  }
  
  ## Gather phenopositivity thresholds
  if (!is.null(thresholds)) {
    
    d$Threshold <- thresholds[d$Marker]
  }
  
  ## Resolve whether to plot thresholds and/or noise cut-offs
  plot_thresholds <-
    !is.null(thresholds) &&
    length(thresholds)>1 &&
    sum(!is.na(d$Threshold))>0
  plot_noise_cutoffs <-
    !is.null(noise) &&
    sum(!is.na(d$Cutoff))>0
  
  ## Resolve plot subtitle
  subtitle <- ''
  if (plot_thresholds) {
    subtitle <- 'phenopositivity thresholds in purple'
  }
  if (plot_noise_cutoffs) {
    subtitle <- paste(
      c(subtitle, 'average cut-offs for background signal in grey'),
      collapse = ', '
    )
  }
  
  ## Initialise plot
  p <- ggplot2::ggplot(
    data = d,
    mapping = ggplot2::aes(
      x   = .data$Compartment,
      y   = .data$Value,
      col = .data$Compartment
    ))
  
  if (plot_thresholds) {
    
    ## Add phenopositivity thresholds to the plot
    p <- p +
      ggiraph::geom_hline_interactive(
        data = d[d$Quantile==0.50, , drop = FALSE],
        mapping = ggplot2::aes(
          yintercept = Cutoff,
          tooltip    = paste0(
            .data$Marker, ' average noise cut-off: ',
            round(x = .data$Cutoff, digits = 2)
          )
        ),
        color     = '#595959',
        alpha     = .6,
        linewidth = .8
      )
  }
  if (plot_noise_cutoffs) {
    
    ## Add noise cut-offs to the plot
    p <- p +
      ggiraph::geom_hline_interactive(
        data = d[d$Quantile==0.50, , drop = FALSE],
        mapping = ggplot2::aes(
          yintercept = Threshold,
          tooltip    = paste0(
            .data$Marker, ' phenopositivity threshold: ',
            round(x = .data$Threshold, digits = 2)
          )
        ),
        color     = 'darkorchid1',
        alpha     = .4,
        linewidth = .8
      )
  }
  
  ## Add signal quartiles to the plot
  p <- p +
    ggplot2::scale_colour_manual(
      values = palette
    ) +
    ggiraph::geom_point_interactive(
      mapping = ggplot2::aes(tooltip = .data$Marker)
    ) +
    ggplot2::geom_line(
    ) +
    ggh4x::facet_nested_wrap( # facet by lineage-vs-state and by specific marker
      ~.data$Type+.data$Marker,
      nrow = 1,
      labeller = ggiraph::labeller_interactive(
        mapping = ggplot2::aes(tooltip = .data$Marker)
      )
    ) +
    ggplot2::theme_grey(
    ) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 12),
      legend.title       = ggplot2::element_blank(),
      legend.position    = 'top',
      legend.margin      = ggplot2::margin(c(0, 0, 0, 0)),
      panel.background   = ggplot2::element_rect(fill = '#fcfcf7'),
      panel.border       = ggplot2::element_rect(
        colour = 'black', fill = NA, linewidth = .3),
      axis.ticks.x       = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_blank(),
      text               = ggplot2::element_text(size = 10),
      legend.text        = ggplot2::element_text(size = 10),
      strip.background   = ggplot2::element_rect(
        fill = 'white', colour = 'darkgrey', linewidth = .35
      ),
      panel.grid.major.y = ggplot2::element_line(
        color = 'darkgrey', linewidth = .2, linetype = 2
      ),
      panel.grid.minor.y = ggplot2::element_line(
        color = 'darkgrey', linewidth = .2, linetype = 2
      ),
      axis.title.x       = ggtext::element_markdown(
        face = 'italic', color = '#595959', size = 10
      )
    ) +
    ggplot2::ylab(
      'Signal quartiles'
    ) +
    ggplot2::xlab(
      paste0('(', subtitle, ')')
    )
  p
}

