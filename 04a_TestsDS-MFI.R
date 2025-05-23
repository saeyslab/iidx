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

## iidx internal module: 04a_TestsDS-MFI.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines functions for differential state testing via MFIs using
# limma with extensions for post-hoc explainability (for a single experimental
# design).


## Function: fit DS-MFI for single compartment ----

ds_mfi_singlefit <- function(
    comps,      # all compartments (metacluster-marker combinations)
    idx_comp,   # index of compartment to test
    mfi,        # sample-wise MFI values per compartment
    weights,    # weights per sample
    experiment, # experiment design from `prep_experiment`
    batches,    # batch per sample
    nbatches,   # number of unique batches
    wconf       # whether confounder is specified
) {
  
  ## Gather inputs for model
  comp <- comps[idx_comp]
  y    <- mfi[idx_comp, , drop = TRUE]
  
  ## Use samples with non-missing outcome values for test
  mask       <- is.na(y)
  na_outcome <- names(y)[mask]
  y          <- y[!mask]
  weights    <- weights[!mask]
  samples    <- names(y)
  
  ## Apply scaling to the sampling weights
  weights <- weights/median(weights)
  
  ## Gather inputs required to fit the model
  dat <-
    experiment$Data[match(samples, experiment$Data$sample_id), , drop = FALSE]
  d <- cbind(
    'y' = y,
    'w' = weights,
    dat
  )
  
  ## Fit model
  suppressMessages({
    suppressWarnings({
      fit <- lmerTest::lmer(
        formula = experiment$Formula,
        data    = d,
        weights = d$w
      )
    })
  })
  
  ## Estimate 95% confidence intervals for batch random intercepts (predictor)
  fi <-
    lme4::fixef(fit)['(Intercept)'] # fixed intercept
  re <-
    as.data.frame(lme4::ranef(fit)) # all random-intercept errors
  re <- re[re$grpvar=='Batch', ] # batch-intercept errors
  z_crit <-
    stats::qnorm(.975) # critical value for CI (for two-tailed Gaussian)
  intercepts <-
    stats::coef(fit)$Batch[, '(Intercept)'] # batch intercepts
  re$ci_min <-
    intercepts-z_crit*re$condsd # CI lower bounds
  re$ci_max <-
    intercepts+z_crit*re$condsd # CI upper bounds
  
  ## Calculate goodness-of-fit estimate per batch
  sr <-
    residuals(fit)**2 # squared residuals
  tot <-
    (d$y-mean(d$y, na.rm = TRUE))**2 # total squares
  rss <-
    sum(sr, na.rm = TRUE) # residual sum of squares (RSS)
  tss <-
    sum(tot, na.rm = TRUE) # total sum of squares (TSS)
  rss_batch <-
    sapply(batches, function(b) sum(sr[d$Batch==b], na.rm = TRUE)) # batch RSS
  names(rss_batch) <- batches
  tss_batch        <-
    sapply(batches, function(b) sum(tot[d$Batch==b], na.rm = TRUE)) # batch TSS
  names(tss_batch) <- batches
  rsq <-
    1-rss/tss # overall R^2
  rsq_batch <- 
    sapply(batches, function(b) 1-rss_batch[[b]]/tss_batch[[b]]) # batch R^2
  names(rsq_batch) <- batches
  
  ## Extract fitted parameters for predictor
  pval  <- summary(fit)$coefficients[, 'Pr(>|t|)'][-1] # p-value
  coeff <- unlist(stats::coef(fit)$Batch[1, -1]) # effect
  
  ## Gather fitted parameters
  res <- data.frame(
    'Compartment'       = comp,
    'logFC'             = sign(coeff[1])*log2(1+abs(coeff[1])),
    'change'            = coeff[1],
    'PValue'            = pval[1],
    'AdjPVal'           = NA,
    'Rsq'               = rsq,
    'logFCConfounder'   = NA,
    'changeConfounder'  = NA,
    'PValueConfounder'  = NA,
    'AdjPValConfounder' = NA
  )
  if (wconf) {
    res['PValueConfounder'] <- pval[2]
    res['logFCConfounder']  <- sign(coeff[2])*log2(1+abs(coeff[2]))
    res['changeConfounder'] <- coeff[2]
  }
  rownames(res) <- NULL
  
  ## Gather random intercept parameters per batch
  random_intercepts <- data.frame(
    'Compartment'           = rep(comp, times = nbatches),
    'Batch'                 = batches,
    'Intercept'             = intercepts,
    'ConfidenceIntervalMin' = re$ci_min,
    'ConfidenceIntervalMax' = re$ci_max
  )
  rownames(random_intercepts) <- NULL
  
  ## Gather goodness-of-fit estimates per batch
  batch_rsq <- data.frame(
    'Compartment' = rep(comp, times = nbatches),
    'Batch'       = batches,
    'Rsq'         = rsq_batch
  )
  rownames(batch_rsq) <- NULL
  
  ## Gather names of sampled excluded due to missing outcome values
  na_outcome        <- list(na_outcome)
  names(na_outcome) <- comp
  
  list(
    'res'               = res,
    'random_intercepts' = random_intercepts,
    'batch_rsq'         = batch_rsq,
    'na_outcome'        = na_outcome
  )
}

## Function: fit DS-MFI for all compartments ----

fit_ds_mfi_model <- function(
    mfi,                # sample-wise MFI values per compartment
    weights,            # weights per sample
    samples,            # sample names
    annotation,         # sample-level annotation
    predictor,          # biological predictor to be modelled
    confounder = NULL,  # biological confounder to be modelled
    famstr     = FALSE, # whether annotation$FamilyID should be used to account
    # for siblings using fixed intercepts
    parallel   = FALSE, # whether to use multi-threading
    verbose    = TRUE   # whether to show progress
) {
  
  ## Only allow one biological covariate/confounder
  wconf <- !is.null(confounder)&&length(confounder)>0
  if (wconf && length(confounder)!=1) {
    stop('One biological confounder at a time is currently allowed')
  }
  
  ## Set up experiment design
  experiment <- prep_experiment(
    samples,
    annotation,
    fixed_effects = c(predictor, confounder),
    random_effects =
      if (famstr) {
        c('Batch', 'FamilyID')
      } else {
        'Batch'
      }, # (batch and maybe family ID modelled via random intercepts)
    force_ls       = TRUE
  )
  na_annotation <- experiment[['NA']]
  
  ## Exclude samples with missing predictor/covariate values
  mask       <- colnames(mfi)%in%na_annotation
  mfi        <- mfi[, !mask, drop = FALSE]
  weights    <- weights[!mask]
  
  ## Establish metacluster-marker combinations as compartments for testing
  comps  <- rownames(mfi)
  ncomps <- length(comps)
  
  ## Gather sample batches
  batches  <- levels(experiment$Data$Batch)
  nbatches <- length(batches)
  
  ## Fit models per compartment
  if (parallel) { # multi-threading enabled
    
    ## Set up parallel processing
    cores <- detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoSNOW(cl)
    
    if (verbose) {
      pb <- utils::txtProgressBar(
        min = 0, max = ncomps, style = 3
      ) # progress bar
      opts <- list(
        'progress' = function(i) utils::setTxtProgressBar(pb, i)
      )
    } else {
      opts <- list()
    }
    out <- foreach(
      idx_comp = seq_len(ncomps),
      .combine = function(res1, res2) {
        list(
          'res' = rbind(res1$res, res2$res),
          'random_intercepts' =
            rbind(res1$random_intercepts, res2$random_intercepts),
          'batch_rsq' = rbind(res1$batch_rsq, res2$batch_rsq),
          'na_outcome' = c(res1$na_outcome, res2$na_outcome)
        )
      },
      .inorder      = TRUE, # (set to FALSE to get more speed-up)
      .options.snow = opts,
      .export       = c('ds_mfi_singlefit'), # required function
      .packages     = c('lme4', 'lmerTest')  # required packages
    ) %dopar% {
      ds_mfi_singlefit(
        comps, idx_comp, mfi, weights, experiment, batches, nbatches, wconf
      )
    }
    if (verbose) {
      close(pb)
    }
    stopCluster(cl)
    
    ## Extract fitted parameters
    res <- out$res
    
    ## Apply multiple testing correction to predictor p-values
    res$AdjPVal <- stats::p.adjust(res$PValue, method = 'BH')
    
    ## Apply multiple testing correction to confounder p-values
    if (wconf) {
      res$AdjPValConfounder <-
        stats::p.adjust(res$PValueConfounder, method = 'BH')
    }
    
    ## Extract random intercepts and their confidence interval estimates
    random_intercepts <- out$random_intercepts
    
    ## Extract goodness-of-fit estimates
    batch_rsq <- out$batch_rsq
    
    rownames(res) <-
      rownames(random_intercepts) <-
      rownames(batch_rsq) <- NULL
    
    ## Extract samples excluded based on missing outcome values
    na_outcome <- out$na_outcome
  } else { # multi-threading disabled
    
    if (verbose) {
      pb <- utils::txtProgressBar(
        min = 0, max = ncomps, style = 3
      ) # progress bar
    }
    
    out <- lapply(
      seq_len(ncomps),
      function(idx_comp) {
        if (verbose) {
          utils::setTxtProgressBar(pb, idx_comp)
        }
        ds_mfi_singlefit(
          comps, idx_comp, mfi, weights, experiment, batches, nbatches, wconf
        )
      }
    )
    if (verbose) {
      close(pb)
    }
    
    ## Extract fitted parameters
    res <- do.call(rbind, lapply(out, function(x) x$res))
    
    ## Apply multiple testing correction to predictor p-values
    res$AdjPVal <- stats::p.adjust(res$PValue, method = 'BH')
    
    ## Apply multiple testing correction to confounder p-values
    if (wconf) {
      res$AdjPValConfounder <-
        stats::p.adjust(res$PValueConfounder, method = 'BH')
    }
    
    ## Extract random intercepts and their confidence interval estimates
    random_intercepts <-
      do.call(rbind, lapply(out, function(x) x$random_intercepts))
    
    ## Extract goodness-of-fit estimates
    batch_rsq <- do.call(rbind, lapply(out, function(x) x$batch_rsq))
    
    rownames(res) <-
      rownames(random_intercepts) <-
      rownames(batch_rsq) <- NULL
    
    ## Extract samples excluded based on missing outcome values
    na_outcome <- out$na_outcome
  }
  
  ## Gather fitted and post-hoc stats
  list(
    'main'              = res,
    'random_intercepts' = random_intercepts,
    'batch_r_squared'   = batch_rsq,
    'na_annotation'     = na_annotation,
    # ^ samples excluded due to missing annotation (for all compartments)
    'na_outcome'        = na_outcome
    # ^ samples excluded due to missing outcome values (per compartment)
  )
}
