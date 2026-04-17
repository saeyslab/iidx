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
    counts,              # cell counts per metacluster and sample
    annotation,          # sample-level annotation
    predictors   = NULL, # biological predictors
    confounders  = c(),  # biological confounders
    interactions = FALSE,# whether to also test for predictor-confounder
    batch_aware  = TRUE, # whether to model Batch as random intercept
    parallel     = FALSE,# whether to use multi-threading
    verbose      = TRUE  # whether to show progress
) {

  ## Save input sample names before any filtering
  input_samples   <- rownames(counts)
  n_input_samples <- length(input_samples)

  ## Gather inputs
  samples <- rownames(counts) # will be filtered later
  cc      <- t(counts)        # samples in columns (matches fit_da_model input)

  ## Resolve predictors, confounders & family structure
  npred  <- length(predictors)
  nconf  <- length(confounders)
  wconf  <- nconf>0
  if (!wconf) {
    interactions <- FALSE
  }
  famstr <- 'FamilyID'%in%colnames(annotation)&&
    length(unique(annotation$FamilyID[!is.na(annotation$FamilyID)]))>1

  if (verbose) {
    message(
      'Batch-adjusted',
      if (famstr) { ', sibling-adjusted ' } else { ' ' },
      'NB-GLMM models for ',
      npred, ' predictors in ',
      nrow(cc), ' compartments with up to ',
      nconf, ' potential ', ifelse(nconf==1, 'confounder', 'confounders'),
      ' will be fitted'
    )
  }

  ## Initialise joint model results & random-effect diagnostics (no confounder)
  res_joint_conf0 <- vector(mode = 'list', length = npred)
  inter_conf0     <- vector(mode = 'list', length = npred) # batch RI
  rsq_conf0       <- vector(mode = 'list', length = npred) # batch R^2
  names(res_joint_conf0) <-
    names(inter_conf0) <-
    names(rsq_conf0) <- predictors

  na_annotation_joint_conf0 <- vector(mode = 'list', length = npred)
  names(na_annotation_joint_conf0) <- predictors

  ## Initialise joint model results & diagnostics (with confounder)
  res_joint_conf1 <- inter_conf1 <- rsq_conf1 <- NULL
  na_annotation_joint_conf1 <- NULL
  if (wconf) {
    res_joint_conf1 <- vector(mode = 'list', length = npred)
    inter_conf1     <- vector(mode = 'list', length = npred)
    rsq_conf1       <- vector(mode = 'list', length = npred)
    names(res_joint_conf1) <-
      names(inter_conf1) <-
      names(rsq_conf1) <- predictors

    na_annotation_joint_conf1 <- vector(mode = 'list', length = npred)
    names(na_annotation_joint_conf1) <- predictors
  }

  ## Fit joint NB-GLMM models (Batch & FamilyID as random intercepts)
  if (verbose) {
    message('Fitting joint NB-GLMM models')
  }
  pb <- utils::txtProgressBar(min = 0, max = npred, style = 3)
  for (idx_pred in seq_along(predictors)) {

    predictor <- predictors[idx_pred]

    ## Fit joint model without biological confounder
    fit_result <- fit_da_model(
      counts      = cc,
      samples     = samples,
      annotation  = annotation,
      predictor   = predictor,
      confounder  = NULL,
      batch_aware = batch_aware,
      famstr      = famstr,
      interaction = FALSE,
      parallel    = parallel,
      verbose     = FALSE
    )
    res_joint_conf0[[idx_pred]]        <- fit_result$main
    inter_conf0[[idx_pred]]            <- fit_result$random_intercepts
    rsq_conf0[[idx_pred]]              <- fit_result$batch_r_squared
    na_annotation_joint_conf0[[idx_pred]] <- fit_result$na_annotation

    if (wconf) {

      ## Determine biological confounders for this predictor
      this_confs <- confounders[confounders!=predictor]
      this_nconf <- length(this_confs)
      if (this_nconf>0) {

        ## Initialise per-confounder storage
        res_joint_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        inter_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        rsq_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        names(res_joint_conf1[[idx_pred]]) <-
          names(inter_conf1[[idx_pred]]) <-
          names(rsq_conf1[[idx_pred]]) <- this_confs

        na_annotation_joint_conf1[[idx_pred]] <-
          vector(mode = 'list', length = this_nconf)
        names(na_annotation_joint_conf1[[idx_pred]]) <- this_confs

        for (idx_conf in seq_along(this_confs)) {

          confounder <- this_confs[idx_conf]

          ## Fit joint model with biological confounder
          fit_result <- fit_da_model(
            counts      = cc,
            samples     = samples,
            annotation  = annotation,
            predictor   = predictor,
            confounder  = confounder,
            batch_aware = batch_aware,
            famstr      = famstr,
            interaction = interactions,
            parallel    = parallel,
            verbose     = FALSE
          )
          res_joint_conf1[[idx_pred]][[idx_conf]] <- fit_result$main
          inter_conf1[[idx_pred]][[idx_conf]]     <- fit_result$random_intercepts
          rsq_conf1[[idx_pred]][[idx_conf]]       <- fit_result$batch_r_squared
          na_annotation_joint_conf1[[idx_pred]][[idx_conf]] <-
            fit_result$na_annotation
        }
      }
    }

    utils::setTxtProgressBar(pb, idx_pred)
  }
  close(pb)

  ## Determine which predictors are continuous
  cont <- sapply(predictors, function(effect) is.numeric(annotation[, effect]))
  names(cont) <- predictors

  ## Gather results: joint fits + batch random-effect diagnostics
  res <- list(
    'joint'                         = res_joint_conf0,
    'random_intercepts'             = inter_conf0,
    'batch_r_squared'               = rsq_conf0,
    'na_annotation'                 = na_annotation_joint_conf0,
    'confounders_joint'             = res_joint_conf1,
    'confounders_random_intercepts' = inter_conf1,
    'confounders_batch_r_squared'   = rsq_conf1,
    'confounders_na_annotation'     = na_annotation_joint_conf1
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
    batch_aware   = TRUE,
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
  if (!wconf) {
    interactions <- FALSE
  }
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
  
  res_conf1   <- NULL
  inter_conf1 <- NULL
  rsq_conf1   <- NULL
  na_annotation_conf1 <- NULL
  na_outcome_conf1    <- NULL
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
          batch_aware = batch_aware,
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
          batch_aware = batch_aware,
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
                batch_aware = batch_aware,
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
                batch_aware = batch_aware,
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
