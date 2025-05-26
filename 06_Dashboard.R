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

## iidx user-level module: 06_Preprocessing.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Opens an interactive dashboard that allows the domain expert to
# explore, filter, interpret and export results of analysis.

message(
  '## iidx: 06_Dashboard.R STARTED ',
  as.character(round(Sys.time()))
)

## Load required packages ----

for (pkg in c(
  'tidyverse', 'ggh4x', 'ggrepel', 'ggExtra', 'ggnewscale', 'grid', 'gridExtra',
  'ggridges', 'flowCore', 'umap', 'RColorBrewer', 'shiny', 'shinyjs',
  'shinythemes', 'shinyWidgets', 'ggiraph', 'scattermore', 'bslib',
  'shinycssloaders', 'shinyBS', 'htmltools', 'bsplus'
)) {
  library(pkg, character.only = TRUE)
}

## Get parameters and auxiliary functions ----

input_files_needed <- FALSE
source('99_SetUp.R')
source('06a_ResultExtraction.R')
source('06b_Plotting.R')
source('06c_InfoText.R')

## Load analysis inputs and results ----

## Analysed lineage and state markers
lineage_markers <- as.vector(channels[idcs_channels_lineage])
state_markers   <- as.vector(channels[idcs_channels_state])

## Features
inputs <- prep_inputs(
  fname_features_postnorm, fname_sample_outliers, annotation, batch_remove
)

## Parameters of fitted statistical models
res_da           <- readRDS(fname_fsom_da)
res_ds_mfi       <- readRDS(fname_fsom_ds_mfi)
res_ds_phenotype <- readRDS(fname_fsom_ds_pheno)

## Unique metaclusters, biological predictors and biological confounders
bfsom_mcs       <- res_da$joint[[1]]$Compartment
mc_n            <- colSums(inputs$counts)
bfsom_preds     <- attributes(res_da)$Predictors
bfsom_confounders <- c('none', get_confounders(res_da, bfsom_preds[1]))

## Metacluster compositions by population
prof <- NULL
if (gating_available) {
  prof <- readRDS(fname_mc_profiles)
}

## Signal quartiles by metacluster
mc_quart <- readRDS(fname_mc_quartiles)

## FlowSOM model and metacluster labels
fsom             <- readRDS(fname_fsom)
metaclustering   <- readRDS(fname_metaclustering)

## Unstained data, noise levels and robustness status per sample
agg_u <- NULL
bfsom_robustness_all <- NULL
bfsom_noise          <- NULL
if (unstained_samples_available) { # unstained FCS files processed
  
  agg_u <- flowCore::read.FCS(fname_agg_unstained) # aggregated unstained signal
  
  bfsom_robustness_all <- readRDS(fname_robustness) # sample robustness status
  bfsom_noise          <- readRDS(fname_noise) # noise levels
  bfsom_noise$Batch    <- factor(
    as.integer(as.character(bfsom_noise$Batch)),
    levels = sort(unique(as.integer(as.character(bfsom_noise$Batch))))
  ) # make sure batch is as factor
  
  ## Use median + 1 MAD as default noise cut-off
  bfsom_noise$Cutoff <- bfsom_noise$Median+bfsom_noise$MAD*1
}

## Set up dashboard frontend ----

## Define auxialiary function to use spinner as loading icon for UI element
f_spinner <- function(ui_element) {
  withSpinner(
    ui_element,
    type    = 8,
    color   = '#3579bc',
    size    = .75,
    hide.ui = FALSE # do not hide previous output while re-calculating
  )
}

## Set a Shiny fluid page layout
ui <- fluidPage(
  
  ## Custom CSS for background colours, slider layout and card appearance
  tags$style(
    type = "text/css",
    "
      .container-fluid {background-color: #fafeff;}
      .irs-max {visibility: hidden !important;}
      .irs-min {visibility: hidden !important;}
      .slider .form-group {
        display:flex; flex-direction:row; width: 100%;
      }
      .slider label {
        margin-right: 2rem; margin-top: 0.6rem; margin-bottom: 0rem;
        align-self: right; text-align: right; flex-basis: 240px;
      }
      .slider .irs {flex-basis: 90%;}
      .card {box-shadow: 0px 0px 0px; border: 1px solid #d4d4d4;}
    "
  ),
  
  ## Message-handler tweak
  tags$script(
    "
      Shiny.addCustomMessageHandler('bfsom_volcano', function(x) {
        Shiny.setInputValue('bfsom_volcano', x);
      });
    "
  ),
  
  ## Use 'flat' design for slider element
  chooseSliderSkin(
    skin  = 'Flat',
    color = c('#3579bc')
  ),
  
  ## Add all scripts used by `shinyjs` to the <head> tag
  useShinyjs(),
  
  ## Use `bsplus` tooltips for interactive help
  use_bs_tooltip(),
  
  ### Section: navigation bar at the top ----
  
  #### Button: user guide ----
  actionButton(
    inputId = 'btn_help',
    label   = '',
    icon    = icon(
      name = 'glyphicon glyphicon-question-sign',
      lib  = 'glyphicon'
    ),
    style = 'position: absolute; top: 5px; right: 45px; z-index:10000;'
  ) %>% 
    bs_embed_tooltip(
      title = 'Display a user guide for this dashboard'
    ),
  
  #### Button: pipeline info ----
  actionButton(
    inputId = 'btn_info',
    label   = '',
    icon    = icon(
      name = 'glyphicon glyphicon-info-sign',
      lib  = 'glyphicon'
    ),
    style = 'position: absolute; top: 5px; right: 130px; z-index:10000;'
  ) %>%
    bs_embed_tooltip(
      title = 'Display info about the differential expression workflow'
    ),
  
  page_navbar(
    bg        = '#fafeff',
    title     = 'Differential expression dashboard',
    fillable  = FALSE,
    underline = TRUE,
    selected  = 'FlowSOM+(G)LMMs', 
    # ^ analysis results are in a single tab, allowing for additional tabs if
    #   the analysis is expanded to include other methodologies in the future
    
    nav_panel(
      title = 'FlowSOM+(G)LMMs',
      style = 'font-size:80%;',
      
      ## Main content (under navigation bar)
      fluidRow(
        
        ### Section: model configuration ----
        
        sidebarPanel(
          
          style = 'background-color: #f7f7f7;',
          width = 3,
          
          #### Dropdown: DE analysis type ----
          selectInput(
            inputId = 'bfsom_select_analysis_type',
            label   = 'Analysis type',
            choices = c(
              'Differential Abundance',
              'Differential State by MFIs',
              'Differential State by Phenopositives'
            ),
            selected = 'Differential Abundance'
          ),
          
          fluidRow(
            column(
              width = 5,
              align = 'left',
              
              #### Dropdown: biological predictor ----
              selectInput(
                inputId  = 'bfsom_select_predictor',
                label    = 'Predictor',
                choices  = bfsom_preds,
                selected = bfsom_preds[1]
              )
            ),
            column(
              width = 2,
              align = 'center',
              
              #### Button: swap predictor and confounder ----
              actionButton(
                inputId = 'bfsom_btn_swap',
                label = '',
                class = 'btn-sm',
                icon  = icon(
                  name = 'glyphicon glyphicon-resize-horizontal',
                  lib  = 'glyphicon'
                ),
                style = paste0(
                  'position: relative; top: 27px; width: 100%; font-size:124%;',
                  ' border-color: #a6a6a6; background-color: #f8feff;'
                )
              )
            ),
            column(
              width = 5,
              align = 'left',
              
              #### Dropdown: biological confounder ----
              selectInput(
                inputId  = 'bfsom_select_confounder',
                label    = 'Covariate',
                choices  = bfsom_confounders,
                selected = 'none'
              )
            )
          ),
          
          ### Section: metacluster/marker configuration ----
          fluidRow(
            column(
              width = 6,
              align = 'left',
              
              ### Dropdown menu: metacluster ----
              selectInput(
                inputId  = 'bfsom_select_mc',
                label    = 'Metacluster',
                choices  = bfsom_mcs,
                selected = bfsom_mcs[1]
              )
            ),
            column(
              width = 6,
              align = 'left',
              
              #### Dropdown menu: state marker for DS ----
              selectInput(
                inputId  = 'bfsom_select_state_marker',
                label    = 'State marker',
                choices  = state_markers,
                selected = state_markers[1]
              )
            )
          ),
          
          ### Section: CSV export ----
          fluidRow(
            column(
              width = 6,
              align = 'center',
              
              #### Button: Export results of selected model ----
              actionButton(
                inputId = 'bfsom_button_export',
                label   = 'Export selected',
                class   = 'btn-sm',
                style   = paste0(
                  'width: 100%; font-size:85%; border-color: #a6a6a6; ',
                  'background-color: #f8feff; margin-top:15px; ',
                  'margin-bottom:-5px;'
                )
              ) %>%
                bs_embed_tooltip(
                  title = paste0(
                    'Export fitted model parameters of the selected test as a ',
                    'CSV file'
                  )
                )
            ),
            column(
              width = 6,
              align = 'center',
              
              #### Button: Export results of all models ----
              actionButton(
                inputId = 'bfsom_button_export_all',
                label   = 'Export all',
                class   = 'btn-sm',
                style   = paste0(
                  'width: 100%; font-size:85%; border-color: #a6a6a6; ',
                  'background-color: #f8feff; margin-top:15px; ',
                  'margin-bottom:-5px;'
                )
              ) %>%
                bs_embed_tooltip(
                  title = paste0(
                    'Export fitted model parameters of all tests for this ',
                    'analysis type in a CSV file'
                  )
                )
            )
          ),
          hr(),
          
          ### Section: robustness filtering configuration ----
          div(
            class = 'slider',
            
            #### Slider: signal filter strength ----
            sliderInput(
              inputId = 'bfsom_robustness_signal_filter_strength',
              label   = 'Signal filter strength:',
              min     = 1L,
              max     = 5L,
              value   = 1L,
              step    = 1L
            ) %>%
              bs_embed_tooltip(
                title = paste0(
                  'How many MADs above median signal of background noise from ',
                  'unstained sample should be considered robust'
                ),
                placement = 'bottom'
              )
          ),
          div(
            class = 'slider',
            
            #### Slider: sample filter strength ----
            sliderInput(
              inputId = 'bfsom_robustness_sample_filter_strength',
              label   = 'Sample filter strength:',
              min     = 0.,
              max     = 1.,
              value   = 0.,
              step    = .05
            ) %>%
              bs_embed_tooltip(
                title = paste0(
                  'Minimum proportion of samples with signal above the noise ',
                  'threshold'
                ),
                placement = 'bottom'
              )
          ),
          hr(),
          
          ### Section: interactive plot configuration ----
          fluidRow(
            column(
              width = 12,
              align = 'left',
              
              div(
                style = 'margin-top: -10px;',
                
                #### Checkbox: filter volcano plot by marker ----
                checkboxInput(
                  inputId = 'bfsom_checkbox_volcano_onlymarker',
                  label   = 'Filter volcano plot by selected state marker',
                  value   = FALSE
                ) %>%
                  bs_embed_tooltip(
                    title = paste0(
                      'Whether to only show differential state hits for the ',
                      'selected state marker in the volcano plot (top-left ',
                      'pane)'
                    )
                  )
              ),
              div(
                style = 'margin-top: -10px;',
                
                #### Checkbox: plot noise cut-offs with MFIs ----
                checkboxInput(
                  inputId = 'bfsom_checkbox_mfi_hit_cutoffs',
                  label   = 'Show noise cut-offs with DS-MFI results',
                  value   = FALSE
                ) %>%
                  bs_embed_tooltip(
                    title = paste0(
                      'Whether to indicate background noise cut-offs per batch',
                      ' when displaying Differential State MFI model results ',
                      '(top-right pane)'
                    )
                  )
              )
            ),
            column(
              width = 4,
              align = 'left',
              
              #### Radio: results view ----
              radioButtons(
                inputId  = 'bfsom_radio_view',
                label    = 'Results view:',
                choices  = c(
                  'By batch',  # colour-code points (feature values) by batch
                  'By cohort', # colour-code points (feature values) by cohort
                  'Binned',    # coarsely bin and colour-code by density
                  'Smooth'     # draw a smooth contour plot
                ),
                selected = 'By batch'
              ) %>%
                bs_embed_tooltip(
                  title = paste0(
                    'Mode of plotting feature values for selected compartment ',
                    '(top-right pane)'
                  )
                )
            ),
            column(
              width = 8,
              align = 'left',
              
              #### Dropdown: scatterplot marker (x-axis) ----
              selectInput(
                inputId = 'bfsom_select_marker1',
                label   = 'Scatterplot marker 1:',
                choices = list(
                  'Lineage' = lineage_markers,
                  'State'   = state_markers
                ),
                selected = lineage_markers[1]
              ),
              
              #### Dropdown: scatterplot marker (y-axis) ----
              selectInput(
                inputId = 'bfsom_select_marker2',
                label   = 'Scatterplot marker 2:',
                choices = list(
                  'Lineage' = lineage_markers,
                  'State'   = state_markers
                ),
                selected = lineage_markers[2]
              )
            ),
            column(
              width = 6,
              align = 'left',
              
              div(
                style = 'margin-top: 0px; margin-bottom: -20px;',
                
                #### Radio: pheno+ vs. noise cut-off in scatterplot ----
                radioButtons(
                  inputId  = 'bfsom_radio_scatter',
                  label    = 'Scatterplot thresholds:',
                  choices  = c(
                    'Phenopositivity',
                    'Noise cut-offs'
                  ),
                  selected = 'Phenopositivity'
                ) %>%
                  bs_embed_tooltip(
                    title = paste0(
                      'Whether to indicate phenopositivity thresholds per ',
                      'marker or background noise cut-offs per marker (if ',
                      'available) in the biaxial expression plot (bottom-left)'
                    )
                  )
              )
            ),
            column(
              width = 6,
              align = 'left',
              
              div(
                style = 'margin-top: 5.2px; margin-bottom: -50px;',
                br(),
                
                #### Checkbox: highlight metacluster in scatterplot ----
                checkboxInput(
                  inputId = 'bfsom_checkbox_show_mc',
                  label   = 'Highlight metacluster',
                  value   = TRUE
                ) %>%
                  bs_embed_tooltip(
                    title = paste0(
                      'Whether to highlight signal from selected metacluster',
                      ' in biaxial expression plot (bottom-left pane)'
                    )
                  ),
                
                div(
                  style = 'margin-top: -50px; margin-bottom: 0px;',
                  br(), br(),
                  
                  #### Checkbox: unstained signal in scatterplot ----
                  checkboxInput(
                    inputId = 'bfsom_checkbox_bg_unstained',
                    label   = 'Show unstained',
                    value   = FALSE
                  ) %>%
                    bs_embed_tooltip(
                      title = paste0(
                        'Whether background in biaxial expression plot ',
                        '(bottom-left pane) should be the background signal ',
                        'from unstained samples'
                      )
                    )
                )
              )
            )
          ),
          hr(),
          
          ### Section: metacluster-population profile ----
          plotOutput(
            outputId = 'bfsom_profile',
            width    = '345px',
            height   = '420px'
          ) %>% f_spinner
        ),
        
        ## Main panel with interactive plots
        mainPanel(
          width = 9,
          
          fluidRow(
            column(
              width = 6,
              
              card(
                card_body(
                  padding = 5
                ),
                
                ### Interactive: volcano plot of leads ----
                girafeOutput(
                  outputId = 'bfsom_volcano',
                  height   = '420px'
                ) %>%
                  f_spinner
              )
            ),
            column(
              width = 6,
              
              card(
                card_body(
                  padding = 5
                ),
                
                ### Interactive: feature values ----
                girafeOutput(
                  outputId = 'bfsom_hit',
                  height   = '420px'
                ) %>%
                  f_spinner
              )
            ),
          ),
          fluidRow(
            column(
              width = 6,
              
              card(
                card_body(
                  padding = 5
                ),
                
                ### Interactive: scatterplot ----
                girafeOutput(
                  outputId = 'bfsom_scatterplot',
                  height   = '390px'
                ) %>%
                  f_spinner
              )
            ),
            column(
              width = 6,
              
              card(
                card_body(
                  padding = 5
                ),
                
                ### Interactive: batch/confounder statistics ----
                girafeOutput(
                  outputId = 'bfsom_batch_stats',
                  height   = '390px'
                ) %>%
                  f_spinner
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              
              card(
                card_body(
                  padding = 5
                ),
                
                ### Interactive: signal quartiles ----
                girafeOutput(
                  outputId = 'bfsom_quartiles',
                  height   = '396px' # (to align with the sidebar)
                ) %>%
                  f_spinner
              )
            )
          )
        )
      )
    )
  )
)

## Set up dashboard backend ----

server <- function(
    input,
    output,
    session
) {
  
  ### Modal: user guide ----
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title = NULL,
      HTML(
        text = paste0(
          '<center><img src="https://raw.githubusercontent.com/davnovak/',
          'davnovak.github.io/main/dashboard_help.png" width="1000px" />',
          '</center>'
        )
      ),
      size      = 'xl',
      easyClose = TRUE,
      footer    = NULL
    ))
  })
  
  ## Modal: info ----
  observeEvent(input$btn_info, {
    showModal(info_modal)
  })
  
  ## Initialise reactive values ----
  react <- reactiveValues()
  
  ## Initialise analysis type choice
  react$bfsom_analysis_type  <- 'Differential Abundance'
  react$bfsom_res            <- res_da
  
  ## Initialise predictors and confounders and default choice
  react$bfsom_preds          <- bfsom_preds
  react$bfsom_pred           <- bfsom_preds[1]
  react$bfsom_confs          <- c('none', bfsom_confounders)
  react$bfsom_conf           <- 'none'
  
  ## Initialise metacluster and state marker choice
  react$bfsom_mc             <- bfsom_mcs[1]
  react$bfsom_state_marker   <- state_markers[1]
  
  ## Default configuration for robustness filtering
  react$bfsom_robustness_       <- 0. # sample filtering disabled by default
  react$bfsom_robustness_signal_filter_strength <-
    1L # 1 MAD from median as default cut-off
  react$bfsom_noise             <- bfsom_noise
  react$bfsom_robustness_rates  <- NULL # robustness rates assigned dynamically
  
  ## Default configuration for volcano plot
  react$bfsom_volcano_onlymarker <- FALSE # do not filter volcano plot by marker
  
  ## Default configuration for feature values plot
  react$bfsom_mfi_hit_cutoffs <- FALSE # do not plot DS-MFI noise cut-offs
  react$bfsom_view <- 'By batch' # colour-code feature values by batch
  
  ## Default configuration for scatterplot
  react$bfsom_marker1 <- as.vector(channels)[idcs_channels_state[1]] # x-axis
  react$bfsom_marker2 <- as.vector(channels)[idcs_channels_state[2]] # y-axis
  react$bfsom_show_mc      <- TRUE    # highlight selected metacluster
  react$bfsom_bg_unstained <- FALSE   # use stained signal for background
  react$bfsom_scatter_th   <- 'pheno' # display phenopositivity thresholds
  
  
  ## Disable unused UI elements ----
  disable('bfsom_btn_swap') # predictor-confounder swap button at beginning
  if (is.null(agg_u)) { # no aggregate unstained data available
    disable('bfsom_checkbox_bg_unstained') # no option to display it
  }
  if (is.null(bfsom_noise)) { # no data for noise cut-offs available
    disable('bfsom_radio_scatter') # no option to display them
  }
  
  ## Observe signal filter strength ----
  observeEvent(
    eventExpr   = input$bfsom_robustness_signal_filter_strength,
    handlerExpr = {
      
      ## Update
      react$bfsom_robustness_signal_filter_strength <-
        input$bfsom_robustness_signal_filter_strength
      
      ## Re-compute signal noise cut-off levels
      average     <- react$bfsom_noise$Median
      deviation   <- react$bfsom_noise$MAD
      n_deviation <- react$bfsom_robustness_signal_filter_strength
      
      react$bfsom_noise$Cutoff <- average+deviation*n_deviation
    }
  )
  
  ## Observe volcano plot filtering ----
  observeEvent(
    eventExpr   = input$bfsom_checkbox_volcano_onlymarker,
    handlerExpr = {
      
      ## Update
      react$bfsom_volcano_onlymarker <- input$bfsom_checkbox_volcano_onlymarker
    }
  )
  
  ## Observe noise cut-off in feature values plot ----
  observeEvent(
    eventExpr   = input$bfsom_checkbox_mfi_hit_cutoffs,
    handlerExpr = {
      
      ## Update
      react$bfsom_mfi_hit_cutoffs <- input$bfsom_checkbox_mfi_hit_cutoffs
    }
  )
  
  ## Observe analysis type ----
  observeEvent(
    eventExpr   = input$bfsom_select_analysis_type,
    handlerExpr = {
      
      ## Note previous setting
      prev <- react$bfsom_analysis_type
      
      ## Update
      react$bfsom_analysis_type <- input$bfsom_select_analysis_type
      
      if (react$bfsom_analysis_type=='Differential Abundance') {
        
        ## Load results
        react$bfsom_res <- res_da
        
        ## Nullify robustness rates
        react$bfsom_robustness_rates <- NULL
        
        ## Disable DS-specific UI elements
        disable('bfsom_robustness_sample_filter_strength')
        disable('bfsom_select_state_marker')
        disable('bfsom_checkbox_volcano_onlymarker')
        disable('bfsom_checkbox_mfi_hit_cutoffs')
      } else if (react$bfsom_analysis_type=='Differential State by MFIs') {
        
        ## Load results
        react$bfsom_res <- res_ds_mfi
        
        ## Resolve robustness rate and enable filtering
        if (!is.null(bfsom_robustness_all)) {
          
          ## Allow to regulate sample filter strength
          enable('bfsom_robustness_sample_filter_strength')
          
          ## Compute robustness rates
          react$bfsom_robustness_rates <-
            robustness_rate_per_test(
              res        = react$bfsom_res,
              robustness = bfsom_robustness_all,
              predictor  = react$bfsom_pred,
              confounder = react$bfsom_conf,
              tidy       = TRUE
            )[[
              as.integer(react$bfsom_robustness_signal_filter_strength)
            ]]
        }
        
        ## Enable DS-specific UI elements
        enable('bfsom_select_state_marker')
        enable('bfsom_checkbox_volcano_onlymarker')
        enable('bfsom_checkbox_mfi_hit_cutoffs')
      } else if (react$bfsom_analysis_type=='Differential State by Phenopositives') {
        
        ## Load results
        react$bfsom_res <- res_ds_phenotype
        
        ## Resolve robustness rate and enable filtering
        if (!is.null(bfsom_robustness_all)) {
          
          ## Allow to regulate sample filter strength
          enable('bfsom_robustness_sample_filter_strength')
          
          ## Compute robustness rates
          react$bfsom_robustness_rates <-
            robustness_rate_per_test(
              res        = react$bfsom_res,
              robustness = bfsom_robustness_all,
              predictor  = react$bfsom_pred,
              confounder = react$bfsom_conf,
              tidy       = TRUE
            )[[
              as.integer(react$bfsom_robustness_signal_filter_strength)
            ]]
        }
        
        ## Enable DS-specific UI elements
        enable('bfsom_select_state_marker')
        enable('bfsom_checkbox_volcano_onlymarker')
        disable('bfsom_checkbox_mfi_hit_cutoffs')
      }
    }
  )
  
  ## Observe unstained signal plotting ----
  observeEvent(
    eventExpr = input$bfsom_checkbox_bg_unstained,
    handlerExpr = {
      
      ## Update
      react$bfsom_bg_unstained <- input$bfsom_checkbox_bg_unstained
    }
  )
  
  ## Observe highlighting metacluster signal ----
  observeEvent(
    eventExpr = input$bfsom_checkbox_show_mc,
    handlerExpr = {
      
      ## Update
      react$bfsom_show_mc <- input$bfsom_checkbox_show_mc
    }
  )
  
  ## Observe confounder choice ----
  observeEvent(
    eventExpr = input$bfsom_select_confounder,
    handlerExpr = {
      
      ## Update
      react$bfsom_conf <- input$bfsom_select_confounder
      
      ## Check possible confounder and predictor choices
      target_confs <- get_confounders(react$bfsom_res, react$bfsom_conf)
      all_preds    <- get_predictors(react$bfsom_res)
      
      ## If confounder and predictor swappable, allow swapping
      if (input$bfsom_select_confounder=='none' ||
          !(
            react$bfsom_conf%in%all_preds &&
            react$bfsom_pred%in%target_confs
          )) {
        disable('bfsom_btn_swap')
      } else {
        enable('bfsom_btn_swap')
      }
    }
  )
  
  ## Observe predictor-confounder swap ----
  observeEvent(
    eventExpr = input$bfsom_btn_swap,
    handlerExpr = {
      
      ## Get values for swap
      tmp_conf <- input$bfsom_select_confounder
      tmp_pred <- input$bfsom_select_predictor
      
      ## Update predictor choice
      updateSelectInput(
        inputId  = 'bfsom_select_predictor',
        selected = tmp_conf
      )
      react$bfsom_pred <- tmp_conf
      
      ## Update possible confounders
      react$bfsom_confs <- c(
        'none',
        get_confounders(react$bfsom_res, react$bfsom_pred)
      )
      
      ## Update confounder choice
      updateSelectInput(
        inputId  = 'bfsom_select_confounder',
        choices  = react$bfsom_confs,
        selected = tmp_pred
      )
      react$bfsom_conf <- tmp_pred
    }
  )
  
  ## Observe sample filter strength ----
  observeEvent(
    eventExpr   = input$bfsom_robustness_sample_filter_strength,
    handlerExpr = {
      
      ## Update
      react$bfsom_robustness_sample_filter_strength <-
        input$bfsom_robustness_sample_filter_strength
    }
  )
  
  ## Observe scatterplot threshold choice ----
  observeEvent(
    eventExpr   = input$bfsom_radio_scatter,
    handlerExpr = {
      
      ## Update
      x <- input$bfsom_radio_scatter
      if (x=='Phenopositivity') {
        x <- 'pheno'
      } else if (x=='Noise cut-offs') {
        x <- 'noise'
      }
      react$bfsom_scatter_th <- x
    }
  )
  
  ## Observe predictor choice ----
  observeEvent(
    eventExpr   = input$bfsom_select_predictor,
    handlerExpr = {
      
      ## Update
      react$bfsom_pred  <- input$bfsom_select_predictor
      
      ## Get possible confounders
      react$bfsom_confs <- c(
        'none',
        get_confounders(react$bfsom_res, react$bfsom_pred)
      )
      
      ## If current confounder choice unavailable, choose 'none'
      if (!react$bfsom_conf%in%react$bfsom_confs) {
        react$bfsom_conf <- 'none'
      }
      
      ## Update possible confounders
      updateSelectInput(
        inputId  = 'bfsom_select_confounder',
        choices  = react$bfsom_confs,
        selected = react$bfsom_conf
      )
    }
  )
  
  ## Observe choice of lead from volcano plot ----
  observeEvent(
    eventExpr   = input$bfsom_volcano_selected,
    handlerExpr = {
      
      ## Get metacluster and (for DS) marker
      comp <- input$bfsom_volcano_selected # selected point label
      if (
        react$bfsom_analysis_type!='Differential Abundance'
      ) {
        mc           <- gsub('[ ].*$', '', comp)
        state_marker <- gsub('^MC[0-9]+[ ]', '', comp)
        state_marker <- gsub(': [0-9].*[$%].*$', '', state_marker)
        
        ## Update state marker choice
        updateSelectInput(
          inputId = 'bfsom_select_state_marker',
          selected = state_marker
        )
      } else {
        
        mc <- comp
      }
      
      ## Update metacluster choice
      updateSelectInput(
        inputId = 'bfsom_select_mc',
        selected = mc
      )
    }
  )
  
  ## Observe metacluster and state marker choice ----
  observe(
    x = {
      
      ## Update
      react$bfsom_mc           <- input$bfsom_select_mc
      react$bfsom_state_marker <- input$bfsom_select_state_marker
    }
  )
  
  ## Observe export of current model stats ----
  observeEvent(
    eventExpr   = input$bfsom_button_export,
    handlerExpr = {
      
      ## Extract parameters of selected model
      st <- get_stats(
        res         = react$bfsom_res,
        predictor   = react$bfsom_pred,
        confounder  = react$bfsom_conf,
        metacluster = react$bfsom_mc,
        marker      = react$bfsom_state_marker,
        robustness  = bfsom_robustness_all,
        wide_format = FALSE
      )
      
      ## Ensure results directory exists
      if (!file.exists(fpath_res06)) {
        dir.create(fpath_res06)
      }
      
      ## Extract name for export CSV file
      name  <- attributes(st)$Name
      fname <- file.path(fpath_res06, paste0(name, '.csv'))
      
      ## Write results to CSV file
      write.csv(x = st, file = fname)
      
      ## Show notification of export
      showNotification(
        ui       = paste0('Stats exported to ', fpath_res06),
        action   = NULL,
        duration = 2,
        type     = 'message'
      )
    }
  )
  
  ## Observe export of all model stats ----
  observeEvent(
    eventExpr   = input$bfsom_button_export_all,
    handlerExpr = {
      
      ## Extract parameters of all models
      st <- get_stats(
        res           = react$bfsom_res,
        predictor     = react$bfsom_pred,
        confounder    = react$bfsom_conf,
        metacluster   = NULL,
        marker        = NULL,
        robustness    = bfsom_robustness_all,
        wide_format   = FALSE,
        state_markers = state_markers
      )
      
      ## Ensure results directory exists
      if (!file.exists(fpath_res06)) {
        dir.create(fpath_res06)
      }
      
      ## Extract name for export CSV file
      name  <- attributes(st)$Name
      fname <- file.path(fpath_res06, paste0(name, '.csv'))
      
      ## Write results to CSV file
      write.csv(x = st, file = fname)
      
      ## Show notification of export
      showNotification(
        ui = paste0('Stats exported to ', fpath_res06),
        action = NULL, duration = 2, type = 'message'
      )
    }
  )
  
  ## Observe metacluster signal quartiles ----
  observe(
    x = {
      
      ## Update interactive plot
      output$bfsom_quartiles <- renderGirafe(
        expr = {
          ggiraph::girafe(
            ggobj = plot_comp_quartiles(
              comp            = react$bfsom_mc,
              quartiles       = mc_quart,
              lineage_markers = lineage_markers,
              state_markers   = state_markers,
              noise           = react$bfsom_noise,
              thresholds      = thresholds
            ),
            width_svg  = 12,
            height_svg = 4
          ) %>%
            ggiraph::girafe_options(
              .,
              ggiraph::opts_hover(
                css = 'fill:#8a0087; stroke:#6e006b; r:5pt;'
              ),
              ggiraph::opts_selection(
                type = 'single', only_shiny = TRUE
              )
            )
        }
      )
    }
  )
  
  ## Observe feature value plot view ----
  observeEvent(
    eventExpr  = input$bfsom_radio_view,
    handlerExpr = {
      
      ## Update
      react$bfsom_view <- input$bfsom_radio_view
    }
  )
  
  ## Observe model & robustness filter choices ----
  observe(
    x = {
      
      ## Check change in analysis type
      bfsom_res <- react$bfsom_res
      
      ## Check change in predictor or confounder
      bfsom_pred <- react$bfsom_pred
      bfsom_conf <- react$bfsom_conf
      
      ## Check changes in robustness filtering
      sample_filter_strength <-
        react$bfsom_robustness_sample_filter_strength
      signal_filter_strength <-
        react$bfsom_robustness_signal_filter_strength
      
      ## Check change in volcano plot filtering by marker
      bfsom_volcano_onlymarker <- react$bfsom_volcano_onlymarker
      marker                   <- react$bfsom_state_marker
      
      ## Recompute robustness rates
      if (
        !is.null(bfsom_robustness_all) &&
        react$bfsom_analysis_type!='Differential Abundance'
      ) {
        
        react$bfsom_robustness_rates <-
          robustness_rate_per_test(
            res        = bfsom_res,
            robustness = bfsom_robustness_all,
            predictor  = bfsom_pred,
            confounder = bfsom_conf,
            tidy       = TRUE
          )[[as.integer(signal_filter_strength)]]
      } else {
        
        react$bfsom_robustness_rates <- NULL
      }
      
      p <- plot_volcano(
        res = bfsom_res,
        marker =
          if (bfsom_volcano_onlymarker) {
            marker
          } else {
            NULL
          },
        predictor         = bfsom_pred,
        confounder        = bfsom_conf,
        robustness_rates  = react$bfsom_robustness_rates,
        robustness_cutoff = react$bfsom_robustness_sample_filter_strength,
        interactive       = TRUE,
        state_markers     = state_markers
      )
      output$bfsom_volcano <- renderGirafe(
        expr = {
          ggiraph::girafe(
            ggobj      = p,
            width_svg  = 6,
            height_svg = 5
          ) %>%
            ggiraph::girafe_options(
              .,
              ggiraph::opts_hover(
                css = 'fill:#8a0087; stroke:#6e006b; r:5pt;'
              ),
              ggiraph::opts_selection(
                type = 'single', only_shiny = TRUE
              )
            )
        }
      )
    }
  )
  
  ## Observe choices for feature values plot ----
  observe(
    x = {
      
      ## Update feature values interactive plot
      output$bfsom_hit <- renderGirafe(
        expr = {
          ggiraph::girafe(
            ggobj = plot_single_hit(
              res                = react$bfsom_res,
              mc                 = react$bfsom_mc,
              marker             = react$bfsom_state_marker,
              noise              = react$bfsom_noise,
              show_noise_cutoffs = react$bfsom_mfi_hit_cutoffs,
              annotation         = annotation,
              data =
                if (
                  react$bfsom_analysis_type=='Differential Abundance'
                ) {
                  inputs$percentages
                } else if (
                  react$bfsom_analysis_type=='Differential State by MFIs'
                ) {
                  inputs$MFIs
                } else if (
                  react$bfsom_analysis_type==
                  'Differential State by Phenopositives'
                ) {
                  inputs$PercPos
                },
              predictor        = react$bfsom_pred,
              confounder       = react$bfsom_conf,
              robustness_rates = react$bfsom_robustness_rates,
              view             = react$bfsom_view,
              n_bins           = 50,
              lim_quantiles    = NULL,
              interactive      = TRUE
            ),
            width_svg  = 6,
            height_svg = 5
          ) %>% ggiraph::girafe_options(
            .,
            ggiraph::opts_hover(
              css = 'fill:#8a0087; stroke:#6e006b; r:5pt;'
            ),
            ggiraph::opts_selection(
              type = 'single', only_shiny = TRUE
            )
          )
        })
    }
  )
  
  ## Observe metacluster choice ----
  observe(
    x = {
      
      ## Update metacluster-population profile plot
      bfsom_mc <- react$bfsom_mc
      output$bfsom_profile <- renderPlot({
        plot_comp_profile(
          prof   = prof,
          counts = inputs$counts,
          comp   = bfsom_mc,
          n      = mc_n
        )
      })
    }
  )
  
  ## Observe x-axis scatterplot marker ----
  observeEvent(
    eventExpr   = input$bfsom_select_marker1,
    handlerExpr = {
      
      ## Update
      react$bfsom_marker1 <- input$bfsom_select_marker1
    }
  )
  
  ## Observe y-axis scatterplot marker ----
  observeEvent(
    eventExpr   = input$bfsom_select_marker2,
    handlerExpr = {
      
      ## Update
      react$bfsom_marker2 <- input$bfsom_select_marker2
    }
  )
  
  ## Observe choices for scatterplot ----
  observe(
    x = {
      
      mc        <- react$bfsom_mc
      m1        <- react$bfsom_marker1
      m2        <- react$bfsom_marker2
      unstained <- react$bfsom_bg_unstained
      bg_only   <- !react$bfsom_show_mc
      th        <- react$bfsom_scatter_th
      
      ## Update interactive scatterplot
      output$bfsom_scatterplot <- renderGirafe(
        expr = {
          ggiraph::girafe(
            ggobj = plot_comp_scatter(
              comp          = mc,
              m1            = m1,
              m2            = m2,
              channels      = channels,
              thresholds    = thresholds,
              noise         =
                if (th=='noise') {
                  react$bfsom_noise
                } else {
                  NULL
                },
              view          = th,
              fsom          = fsom,
              agg_unstained = agg_u,
              bg_unstained  = unstained,
              bg_only       = bg_only
            ),
            width_svg  = 7,
            height_svg = 5
          ) %>%
            ggiraph::girafe_options(
              .,
              ggiraph::opts_hover(
                css = 'fill:#8a0087; stroke:#6e006b; r:5pt;'
              ),
              ggiraph::opts_selection(
                type = 'single', only_shiny = TRUE
              )
            )
        }
      )
    }
  )
  
  ## Observe choices for batch/confounder stats ----
  observe(
    x = {
      
      res          <- react$bfsom_res
      model        <- react$bfsom_pred
      mc           <- react$bfsom_mc
      state_marker <- react$bfsom_state_marker
      conf         <- react$bfsom_conf
      
      ## Update interactive batch/confounder stats plot
      output$bfsom_batch_stats <- renderGirafe(
        expr = {
          ggiraph::girafe(
            ggobj = plot_batch_stats(
              res        = res,
              mc         = mc,
              marker     = state_marker,
              annotation = annotation,
              predictor  = model,
              confounder =
                if (conf=='none') {
                  NULL
                } else {
                  conf
                }
            ),
            width_svg  = 7,
            height_svg = 5
          ) %>%
            ggiraph::girafe_options(
              .,
              ggiraph::opts_hover(
                css = 'fill:#8a0087; stroke:#6e006b; r:5pt;'
              ),
              ggiraph::opts_selection(
                type = 'single', only_shiny = TRUE
              )
            )
        }
      )
    }
  )
}

## Launch dashboard in browser (port 3030) ----

runApp(
  appDir = list(
    'ui'     = ui,
    'server' = server
  ),
  launch.browser = TRUE,
  port = 3030
)

message(
  '## iidx: 06_Dashboard.R FINISHED ',
  as.character(round(Sys.time()))
)


