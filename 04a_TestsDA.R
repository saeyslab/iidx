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

## iidx internal module: 04a_TestsDA.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines functions for differential abundance testing using edgeR
# with extensions for post-hoc explainability (for a single experimental
# design).


## Function: fit a DA model ----

fit_da_model <- function(
    counts,             # cell counts per metacluster and sample
    samples,            # sample names
    annotation,         # sample-level annotation
    predictor,          # biological predictor to be modelled
    confounders = NULL, # confounders: Batch and/or 1 biological variable
    famstr = FALSE,     # whether annotation$FamilyID should be used to account
    # for siblings using fixed intercepts
    interaction = FALSE # whether to also test potential interaction between
    # predictor and biological confounder
) {
  
  ## Resolve specification of biological confounders
  wconf <- !is.null(confounders)
  if (wconf) {
    if (length(confounders[!confounders%in%'Batch'])>1) {
      stop('Only one biological confounder at a time is currently allowed')  
    }
  } else {
    if (interaction) {
      stop(
        'Cannot model interaction term in the absence ',
        'of a biological confounder'
      )
    }
  }
  
  ## Check additional covariate to add (for post-hoc testing)
  ecov <- Sys.getenv('IIDX_EXTRA_COVARIATE')
  if (ecov=='') {
    ecov <- NULL
  }
  
  ## Set up experimental design
  experiment <- prep_experiment(
    files          = samples,
    annotation     = annotation,
    fixed_effects  = 
      if (famstr) {
        c(predictor, confounders, ecov, 'FamilyID')
      } else {
        c(predictor, confounders, ecov)
      },
    random_effects = c(), # no random effects allowed in edgeR model
    interactions   = FALSE
  )
  design     <- experiment$Design
  na_samples <- experiment[['NA']]
  
  ## Exclude samples with missing value for predictor or chosen covariate
  mask <- colnames(counts)%in%na_samples
  counts  <- counts[, !mask, drop = FALSE]
  samples <- samples[!mask]
  
  ## Establish metaclusters as compartments for testing
  comps  <- rownames(counts)
  
  ## Calculate trimmed mean of M-values normalisation factors
  nf <- edgeR::calcNormFactors(counts, method = 'TMM')
  
  ## Represent abundances as a DGEList object
  y <- edgeR::DGEList(counts, norm.factors = nf)
  
  ## Estimate negative-binomial dispersions
  y_disp <- edgeR::estimateDisp(y, design, robust = TRUE)
  
  ## Fit negative-binomial model per compartment
  suppressWarnings({
    fit <- edgeR::glmFit(y_disp, design)
  })
  
  ## Conduct likelihood ratio tests for predictor
  lrt    <- edgeR::glmLRT(fit, coef = 2)
  
  ## Extract fitted parameters for predictor
  logfcs <- lrt$table$logFC                   # log fold changes
  fcs    <- sign(logfcs)*(2^abs(logfcs))      # fold changes
  p      <- lrt$table$PValue                  # raw p-values
  p_adj  <- stats::p.adjust(p, method = 'BH') # adjusted p-values
  names(logfcs) <-
    names(fcs) <-
    names(p) <-
    names(p_adj) <- comps
  
  ## Initialise fitted params for confounder and interaction
  logfcs_conf  <- fcs_conf <- p_conf <- p_adj_conf <- NA
  inter_logfcs <- inter_fcs <- inter_p <- inter_p_adj <- NA
  bio_conf    <- confounders[confounders!='Batch'] # do not consider batch here
  cont_conf   <- NULL # continuous confounder?
  conf_spec   <- !is.null(bio_conf) && length(bio_conf)>0 # confounder specified
  if (conf_spec) {
    
    ## Conduct likelihood ratio tests and extract fitted params for confounder
    lrt_conf    <- edgeR::glmLRT(fit, coef = ncol(fit$design))
    
    logfcs_conf <- lrt_conf$table$logFC                     # log fold changes
    fcs_conf    <- sign(logfcs_conf) * (2^abs(logfcs_conf)) # fold changes
    p_conf      <- lrt_conf$table$PValue                    # raw p-values
    p_adj_conf  <- stats::p.adjust(p_conf, method = 'BH')   # adjusted p-values
    names(logfcs_conf) <-
      names(fcs_conf) <-
      names(p_conf) <-
      names(p_adj_conf) <- comps
    
    cont_conf <- is.numeric(experiment$Data[, bio_conf])
    
    if (interaction) {
      
      ## Train models with interaction term
      experiment <- prep_experiment(
        files          = samples,
        annotation     = annotation,
        fixed_effects  = 
          if (famstr) {
            c(predictor, confounders, 'FamilyID')
          } else {
            c(predictor, confounders)
          },
        random_effects = c(), # no random effects allowed in edgeR model
        interactions = TRUE
      )
      design     <- experiment$Design
      y_disp <- edgeR::estimateDisp(y, design, robust = TRUE)
      suppressWarnings({
        fit <- edgeR::glmFit(y_disp, design)
      })
      
      ## Conduct likelihood ratio tests for interaction term
      lrt    <- edgeR::glmLRT(fit, coef = ncol(fit))
      
      ## Extract fitted parameters for interaction term
      inter_logfcs <- lrt$table$logFC                   
      inter_fcs    <- sign(inter_logfcs)*(2^abs(inter_logfcs))      
      inter_p      <- lrt$table$PValue                  
      inter_p_adj  <- stats::p.adjust(inter_p, method = 'BH') 
      names(inter_logfcs) <-
        names(inter_fcs) <-
        names(inter_p) <-
        names(inter_p_adj) <- comps
    }
  }
  
  ## Gather fitted and post-hoc stats
  res <- data.frame(
    'Compartment'        = comps,
    'logFC'              = logfcs,
    'FC'                 = fcs,
    'PValue'             = p,
    'AdjPVal'            = p_adj,
    
    'logFCConfounder'    = logfcs_conf,
    'FCConfounder'       = fcs_conf,
    'PValueConfounder'   = p_conf,
    'AdjPValConfounder'  = p_adj_conf,
    
    'logFCInteraction'   = inter_logfcs,
    'FCInteraction'      = inter_fcs,
    'PValueInteraction'  = inter_p,
    'AdjPValInteraction' = inter_p_adj
  )
  
  ## Add metadata
  attributes(res)$Continuous <- is.numeric(experiment$Data[, predictor])
  attributes(res)$ContinuousConfounder <- cont_conf
  attributes(res)$na_annotation <-
    na_samples # samples excluded due to missing annotation
  attributes(res)$confounder_specified <- conf_spec
  attributes(res)$IIDX_EXTRA_CONFOUNDER <- ecov
  
  res
}
