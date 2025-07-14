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

## iidx internal module: 04b_FullExperiments.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines functions for differential expression experiments using
# differential abundance, differential state via MFIs and differential state via
# phenopositivity (for all requested experimental designs).

source('04a_TestsDA.R')
source('04a_TestsDS-MFI.R')
source('04a_TestsDS-Pheno.R')

## Function: perform all DA tests ----

test_da <- function(
    counts,               # cell counts per metacluster and sample
    annotation,           # sample-level annotation
    predictors   = NULL,  # biological predictors
    confounders  = c(),   # biological confounders
    interactions = FALSE, # whether to also test for predictor-confounder
                          # interactions
    verbose      = TRUE   # whether to show progress
) {
  
  ## Save input sample names before any filtering
  input_samples   <- rownames(counts)
  n_input_samples <- length(input_samples)
  
  ## Gather inputs
  samples         <- rownames(counts) # will be filtered later
  cc              <- t(counts) # (edgeR requires sample names in columns)
  
  ## Resolve predictors, confounders, batches & family structure
  batches <- unique(annotation$Batch[annotation$FileName%in%samples])
  nba     <- length(batches)
  npred   <- length(predictors)
  nconf   <- length(confounders)
  wconf   <- nconf>0
  famstr  <- 'FamilyID'%in%colnames(annotation)&&
    length(unique(annotation$FamilyID[!is.na(annotation$FamilyID)]))>1
  
  if (verbose) {
    
    message(
      'Batch-adjusted',
      if (famstr) { ', sibling-adjusted '} else { ' ' },
      'edgeR models for ',
      npred, ' predictors in ',
      nrow(cc), ' compartments with up to ',
      nconf, ' potential ', ifelse(nconf==1, 'confounder', 'confounders'),
      ' will be fitted'
    )
  }
  
  ## Initialise joint & batch-wise models without confounders
  res_joint_conf0 <-
    vector(mode = 'list', length = npred) # all batches, without confounders
  names(res_joint_conf0) <- predictors
  res_joint_conf1 <- NULL
  res_batch_conf0 <-
    vector(mode = 'list', length = nba) # by batch, without confounders
  names(res_batch_conf0) <- batches
  res_batch_conf1 <- NULL
  
  ## Initialise excluded samples for joint & batch-wise tests w/out confounders
  na_annotation_joint_conf0 <- vector(mode = 'list', length = npred)
  na_annotation_joint_conf1 <- NULL
  na_annotation_batch_conf0 <- vector(mode = 'list', length = nba)
  na_annotation_batch_conf1 <- NULL
  names(na_annotation_joint_conf0) <- predictors
  names(na_annotation_batch_conf0) <- paste0('Batch', batches)
  
  ## Initialise joint & batch-wise models with confounders
  if (wconf) {
    
    res_joint_conf1 <-
      vector(mode = 'list', length = npred) # all batches, with confounders
    res_batch_conf1 <-
      vector(mode = 'list', length = nba) # by batch, with confounders
    names(res_joint_conf1) <- predictors
    names(res_batch_conf1) <- paste0('Batch', batches)
    
    na_annotation_joint_conf1 <-
      vector(mode = 'list', length = npred)
    na_annotation_batch_conf1 <-
      vector(mode = 'list', length = nba)
    names(na_annotation_joint_conf1) <- predictors
    names(na_annotation_batch_conf1) <- paste0('Batch', batches)
  }
  
  ## Fit joint models (multivariate models with batch encoding as fixed effects)
  if (verbose) {
    message('Fitting joint models')
  }
  pb <- utils::txtProgressBar(
    min = 0, max = npred, style = 3
  ) # progress bar
  for (idx_pred in seq_along(predictors)) {
    
    predictor <- predictors[idx_pred]
    data <- cc
    
    ## Fit joint model without biological covariate term
    fit <- fit_da_model(
      counts      = data,
      samples     = samples,
      annotation  = annotation,
      predictor   = predictor,
      confounders = c('Batch'),
      famstr      = famstr,
      interaction = FALSE
    )
    
    ## Store fitted params and samples excluded due to missing labels
    res_joint_conf0[[idx_pred]] <- fit
    na_annotation_joint_conf0[[idx_pred]] <-
      attributes(fit)$na_annotation
    attributes(fit)$na_annotation <- NULL
    
    if (wconf) {
      
      ## Determine biological confounders for given predictor
      this_confs <- confounders[confounders!=predictor]
      this_nconf <- length(this_confs)
      if (this_nconf>0) {
        
        ## Initialise models for this predictor per confounder
        res_joint_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        names(res_joint_conf1[[idx_pred]]) <- this_confs
        
        na_annotation_joint_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        names(na_annotation_joint_conf1[[idx_pred]]) <- this_confs
        
        for (idx_conf in seq_along(this_confs)) {
          
          confounder <- this_confs[idx_conf]
          
          ## Fit joint model with biological covariate term
          fit <- fit_da_model(
            counts      = data,
            samples     = samples,
            annotation  = annotation,
            predictor   = predictor,
            confounders = c('Batch', confounder),
            famstr      = famstr,
            interaction = interactions
          )
          
          ## Store fitted params and samples excluded due to missing labels
          res_joint_conf1[[idx_pred]][[idx_conf]] <- fit
          na_annotation_joint_conf1[[idx_pred]][[idx_conf]] <-
            attributes(fit)$na_annotation
          attributes(fit)$na_annotation <- NULL
        }
      }
    }
    
    utils::setTxtProgressBar(pb, idx_pred)
  }
  close(pb)
  
  ## Fit batch-wise models (input limited to samples in one batch)
  if (verbose) {
    message('Fitting single-batch models for ', nba, ' batches')
  }
  pb <- utils::txtProgressBar(
    min = 0, max = npred*nba, style = 3
  ) # progress bar
  counter <- 0
  for (idx_ba in seq_along(batches)) {
    
    batch <- batches[idx_ba]
    
    ## Initialise models for this predictor per batch w/out confounders
    res_batch_conf0[[idx_ba]] <-
      vector(mode = 'list', length = npred)
    names(res_batch_conf0[[idx_ba]]) <- predictors
    
    ## Initialise excluded samples for tests per batch w/out confounders
    na_annotation_batch_conf0[[idx_ba]] <-
      vector(mode = 'list', length = npred)
    names(na_annotation_batch_conf0[[idx_ba]]) <- predictors
    
    if (wconf) {
      
      ## Initialise models for this predictor per batch w/ confounders
      res_batch_conf1[[idx_ba]] <-
        vector(mode = 'list', length = npred)
      names(res_batch_conf1[[idx_ba]]) <- predictors
      
      ## Initialise excluded samples for tests per batch w/ confounders
      na_annotation_batch_conf1[[idx_ba]] <-
        vector(mode = 'list', length = npred)
      names(na_annotation_batch_conf1[[idx_ba]]) <- predictors
    }
    
    for (idx_pred in seq_along(predictors)) {
      
      predictor   <- predictors[idx_pred]
      batch_annot <- annotation[
        annotation$Batch==batch &
          annotation$FileName%in%colnames(cc), ,
        drop = FALSE
      ]
      batch_samples <- batch_annot$FileName
      
      data <- cc[, batch_samples, drop = FALSE]
      
      ## Fit batch-restricted model w/out biological covariate term
      fit <- fit_da_model(
        counts      = data,
        samples     = batch_samples, # input samples restricted to 1 batch
        annotation  = batch_annot,
        predictor   = predictor,
        confounders = NULL,
        famstr      = famstr,
        interaction = FALSE
      )
      
      ## Store fitted params and samples excluded due to missing labels
      res_batch_conf0[[idx_ba]][[idx_pred]] <- fit
      na_annotation_batch_conf0[[idx_ba]][[idx_pred]] <-
        attributes(fit)$na_annotation
      attributes(fit)$na_annotation <- NULL
      
      if (wconf) {
        
        ## Determine biological confounders for this predictor
        this_confs <- confounders[confounders!=predictor]
        this_nconf <- length(this_confs)
        if (this_nconf>0) {
          
          ## Initialise models for this predictor per confounder
          res_batch_conf1[[idx_ba]][[idx_pred]] <-
            vector(mode = 'list', length = this_nconf)
          names(res_batch_conf1[[idx_ba]][[idx_pred]]) <- this_confs
          
          na_annotation_batch_conf1[[idx_ba]][[idx_pred]] <-
            vector(mode = 'list', length = this_nconf)
          names(na_annotation_batch_conf1[[idx_ba]][[idx_pred]]) <- this_confs
          
          for (idx_conf in seq_along(this_confs)) {
            
            confounder <- this_confs[idx_conf]
            
            ## Fit batch-restricted model w/ biological covariate term
            fit <- fit_da_model(
              counts      = data,
              samples     = batch_samples, # input samples restricted to 1 batch
              annotation  = batch_annot,
              predictor   = predictor,
              confounders = confounder,
              famstr      = famstr,
              interaction = FALSE # no interaction term in batch-restricted
                                  # setting
            )
            
            ## Store fitted params and samples excluded due to missing labels
            res_batch_conf1[[idx_ba]][[idx_pred]][[idx_conf]] <- fit
            na_annotation_batch_conf1[[idx_ba]][[idx_pred]][[idx_conf]] <-
              attributes(fit)$na_annotation
            attributes(fit)$na_annotation <- NULL
          }
        }
      }
      counter <- counter + 1
      utils::setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  
  ## Determine which predictors are continuous
  cont <- sapply(predictors, function(effect) is.numeric(annotation[, effect]))
  names(cont) <- predictors
  
  ## Gather results for models with and without confounder, joint and by batch
  res <- list(
    'joint'                           = res_joint_conf0,
    'batch'                           = res_batch_conf0,
    'na_annotation'                   = na_annotation_joint_conf0,
    'na_annotation_batch'             = na_annotation_batch_conf0,
    'confounders_joint'               = res_joint_conf1,
    'confounders_batch'               = res_batch_conf1,
    'confounders_na_annotation'       = na_annotation_joint_conf1,
    'confounders_na_annotation_batch' = na_annotation_batch_conf1
  )
  
  ## Add metadata
  attributes(res)$AnalysisType   <- 'Differential Abundance'
  attributes(res)$BaseModel      <- 'FlowSOM'
  attributes(res)$Predictors     <- predictors
  attributes(res)$Continuous     <- cont
  attributes(res)$Annotation     <- annotation
  attributes(res)$InputSamples   <- input_samples
  attributes(res)$NInputSamples  <- n_input_samples
  attributes(res)$FamilyAdjusted <- famstr
  
  res
}

## Function: perform all DS tests ----

test_ds <- function(
    annotation,            # sample-level annotation
    phenopos      = NULL,  # phenopositivity rates per compartment per sample
    mfi           = NULL,  # MFI values per compartment per sample
    counts        = NULL,  # counts per compartment per sample
    predictors    = NULL,  # biological predictors
    confounders   = c(),   # biological confounders
    interactions  = FALSE, # whether to also test for predictor-confounder
                           # interactions
    state_markers = NULL,  # all state markers
    parallel      = FALSE, # whether to use multi-threading
    verbose       = TRUE   # whether to show progress
) {
  
  ## Resolve type of outcome
  if (!is.null(phenopos) && !is.null(mfi)) {
    stop('Either `phenopos` or `mfi` must be specified, not both')
  }
  if (!is.null(mfi) && is.null(counts)) {
    stop('If `mfi` is specified, `counts` must be specified as well')
  }
  outcome <- ifelse(!is.null(phenopos), 'phenopos', 'mfi')
  
  ## Gather outcome values
  if (outcome=='phenopos') {
    
    ## Save input sample names before any filtering
    input_samples   <- rownames(phenopos)
    n_input_samples <- length(input_samples)
    
    str_model       <- 'Beta-GLMM'
    samples         <- rownames(phenopos) # will get filtered later
    perc_pos        <- t(phenopos) # phenopos rates with samples in columns
    
    ## Make sure only state markers used
    if (!is.null(state_markers)) {
      
      markers  <- gsub('^MC[0-9]+[ ]', '', rownames(perc_pos))  
      perc_pos <- perc_pos[markers%in%state_markers, , drop = FALSE]
    }
    
    comps           <- rownames(perc_pos) # compartments
  } else if (outcome=='mfi') {
    
    ## Save input sample names before any filtering
    input_samples   <- rownames(mfi)
    n_input_samples <- length(input_samples)
    
    str_model       <- 'LMM'
    samples         <- rownames(mfi) # will get filtered later
    mfi_signal      <- t(mfi) # median signal values with samples in columns
    
    ## Make sure only state markers used
    if (!is.null(state_markers)) {
      
      markers    <- gsub('^MC[0-9]+[ ]', '', rownames(mfi_signal))  
      mfi_signal <- mfi_signal[markers%in%state_markers, , drop = FALSE]
    }
    
    comps           <- rownames(mfi_signal) # compartments
  }
  
  ## Determine sample weights
  weights <- rowSums(counts)
  
  ## Resolve predictors, their covariates, batches & family structure
  nconf <- length(confounders)
  wconf <- nconf>0
  npred <- length(predictors)
  ncomp <- length(comps)
  batches <- unique(
    as.integer(annotation$Batch[annotation$FileName%in%samples])
  )
  nba   <- length(batches)
  famstr  <- 'FamilyID'%in%colnames(annotation)&&
    length(unique(annotation$FamilyID[!is.na(annotation$FamilyID)]))>1
  if (verbose) {
    
    message(
      'Batch-adjusted',
      if (famstr) { ', sibling-adjusted '} else { ' ' },
      str_model,
      ' models for ',
      npred, ' predictors in ',
      ncomp, ' compartments with up to ',
      nconf, ' potential ', ifelse(nconf==1, 'confounder', 'confounders'),
      ' will be fitted'
    )
  }
  
  ## Initialise fits, intercepts per batch and R^2 per batch w/out confounders
  res_conf0   <- vector(mode = 'list', length = npred)
  inter_conf0 <- vector(mode = 'list', length = npred)
  rsq_conf0   <- vector(mode = 'list', length = npred)
  names(res_conf0) <-
    names(inter_conf0) <-
    names(rsq_conf0) <- predictors
  
  ## Initialise names of excluded samples for tests w/out confounders
  na_annotation_conf0 <- vector(mode = 'list', length = npred)
  na_outcome_conf0    <- vector(mode = 'list', length = npred)
  names(na_annotation_conf0) <-
    names(na_outcome_conf0) <- predictors
  
  if (wconf) { # biological confounder specified
    
    ## Initialise fits, intercepts per batch and R^2 per batch w/ confounders
    res_conf1   <- vector(mode = 'list', length = npred)
    inter_conf1 <- vector(mode = 'list', length = npred)
    rsq_conf1   <- vector(mode = 'list', length = npred)
    names(res_conf1) <-
      names(inter_conf1) <-
      names(rsq_conf1) <- predictors
    
    ## Initialise names of excluded samples for tests w/out confounders
    na_annotation_conf1 <- vector(mode = 'list', length = npred)
    na_outcome_conf1    <- vector(mode = 'list', length = npred)
    names(na_annotation_conf1) <-
      names(na_outcome_conf1) <- predictors
  }
  
  ## Iterate over predictors
  for (idx_pred in seq_along(predictors)) {
    
    predictor <- predictors[idx_pred]
    if (verbose) {
      message('- predictor: ', predictor, ' (', idx_pred, '/', npred, ')')
    }
    
    ## Fit models w/out confounder
    this_res <-
      if (outcome=='phenopos') {
        
        fit_ds_pheno_model(
          phenopos    = perc_pos,
          weights     = weights,
          samples     = samples,
          annotation  = annotation,
          predictor   = predictor,
          confounder  = NULL,
          famstr      = famstr,
          interaction = FALSE,
          parallel    = parallel,
          verbose     = verbose
        )
      } else if (outcome=='mfi') {
        
        fit_ds_mfi_model(
          mfi         = mfi_signal,
          weights     = weights,
          samples     = samples,
          annotation  = annotation,
          predictor   = predictor,
          confounder  = NULL,
          famstr      = famstr,
          interaction = FALSE,
          parallel    = parallel,
          verbose     = verbose
        )
      }
    res_conf0[[predictor]]   <- this_res$main
    inter_conf0[[predictor]] <- this_res$random_intercepts
    rsq_conf0[[predictor]]   <- this_res$batch_r_squared
    
    na_outcome_conf0[[predictor]]    <- this_res$na_outcome
    na_annotation_conf0[[predictor]] <- this_res$na_annotation
    
    if (wconf) {
      
      ## Determine biological confounders for this predictor
      this_confs <- confounders[confounders!=predictor]
      this_nconf <- length(this_confs)
      if (this_nconf>0) {
        
        ### Initialise fits, intercepts and R^2 per confounder
        res_conf1[[idx_pred]]   <- vector(mode = 'list', length = this_nconf)
        inter_conf1[[idx_pred]] <- vector(mode = 'list', length = this_nconf)
        rsq_conf1[[idx_pred]]   <- vector(mode = 'list', length = this_nconf)
        names(res_conf1[[idx_pred]]) <-
          names(inter_conf1[[idx_pred]]) <-
          names(rsq_conf1[[idx_pred]]) <- this_confs
        
        ## Initialise names of excluded samples for each test per confounder
        na_outcome_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        na_annotation_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        names(na_outcome_conf1[[idx_pred]]) <-
          names(na_annotation_conf1[[idx_pred]]) <- this_confs
        
        for (idx_conf in seq_along(this_confs)) {
          
          confounder <- this_confs[idx_conf]
          
          ### Fit models with biological confounder
          this_res <-
            if (outcome=='phenopos') {
              
              fit_ds_pheno_model(
                phenopos    = perc_pos,
                samples     = samples,
                weights     = weights,
                annotation  = annotation,
                predictor   = predictor,
                confounder  = confounder,
                famstr      = famstr,
                interaction = interactions,
                parallel    = parallel,
                verbose     = verbose
              )
            } else if (outcome=='mfi') {
              
              fit_ds_mfi_model(
                mfi         = mfi_signal,
                weights     = weights,
                samples     = samples,
                annotation  = annotation,
                predictor   = predictor,
                confounder  = confounder,
                famstr      = famstr,
                interaction = interactions,
                parallel    = parallel,
                verbose     = verbose
              )
            }
          res_conf1[[predictor]][[confounder]]   <- this_res$main
          inter_conf1[[predictor]][[confounder]] <- this_res$random_intercepts
          rsq_conf1[[predictor]][[confounder]]   <- this_res$batch_r_squared
          na_outcome_conf1[[predictor]][[confounder]] <-
            this_res$na_outcome
          na_annotation_conf1[[predictor]][[confounder]] <-
            this_res$na_annotation
        }
      }
    }
  }
  
  ## Determine which predictors are continuous
  cont <- sapply(
    predictors,
    function(effect) is.numeric(annotation[, effect])
  )
  names(cont) <- predictors
  
  ## Gather results for models with and without confounder
  res <- list(
    'main'                          = res_conf0,
    'random_intercepts'             = inter_conf0,
    'batch_r_squared'               = rsq_conf0,
    'na_annotation'                 = na_annotation_conf0,
    'na_outcome'                    = na_outcome_conf0,
    'confounders_main'              = res_conf1,
    'confounders_random_intercepts' = inter_conf1,
    'confounders_batch_r_squared'   = rsq_conf1,
    'confounders_na_annotation'     = na_annotation_conf1,
    'confounders_na_outcome'        = na_outcome_conf1
  )
  
  ## Add metadata
  attributes(res)$AnalysisType <- 
    paste0(
      'Differential State (',
      ifelse(outcome=='phenopos', 'Phenopositivity', 'MFI'), ')'
    )
  attributes(res)$BaseModel   <- 'FlowSOM'
  attributes(res)$Predictors  <- predictors
  attributes(res)$Confounders <- confounders
  attributes(res)$Continuous  <- cont
  attributes(res)$Annotation  <- annotation
  attributes(res)$InputSamples   <- input_samples
  attributes(res)$NInputSamples  <- n_input_samples
  attributes(res)$FamilyAdjusted <- famstr
  
  res
}
