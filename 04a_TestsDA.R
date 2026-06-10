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
# Description: Defines functions for differential abundance testing using
# NB-GLMM (glmmTMB) with extensions for post-hoc explainability (for a
# single experimental design). FamilyID (family structure) is modelled as a
# random intercept, as is Batch. TMM normalisation factors are retained as an
# offset on the log scale, as in the previous edgeR implementation.


## Function: fit NB-GLMM DA model for single metacluster ----

da_singlefit <- function(
    comps,       # all compartments (metaclusters)
    idx_comp,    # index of compartment to test
    counts,      # count matrix (MCs x samples)
    nf,          # TMM normalisation factors per sample
    experiment,  # experiment design from prep_experiment
    predictor,   # biological predictor column name
    confounder,  # biological confounder column name (or NULL)
    batch_aware, # whether to model Batch as random intercept
    famstr,      # whether to model FamilyID as random intercept
    batches,     # unique batch levels (character)
    nbatches,    # number of unique batches
    wconf,       # whether confounder is specified
    interaction  # whether to test predictor:confounder interaction
) {
  # fun: fit NB-GLMM for one metacluster, extract LRT results, batch
  #      random intercepts w/ 95% CI & batch-wise R^2 goodness-of-fit

  comp    <- comps[idx_comp]
  count_g <- as.integer(counts[comp, ])
  names(count_g) <- colnames(counts)

  ## Use samples with non-missing counts
  mask    <- is.na(count_g)
  count_g <- count_g[!mask]
  nf_g    <- nf[!mask]
  samples <- names(count_g)

  ## Gather annotation for samples present in this fit
  dat <- experiment$Data[
    match(samples, as.character(experiment$Data$sample_id)), , drop = FALSE
  ]
  d <- data.frame('y' = count_g, 'nf' = nf_g, dat)

  ## Ensure FamilyID is unordered factor (numeric IDs need explicit coercion)
  if ('FamilyID'%in%colnames(d)) {
    d$`FamilyID` <- factor(d$`FamilyID`, ordered = FALSE)
  }

  ## Build random-effects part of formula
  re <- c(
    if (batch_aware) 'Batch',
    if (famstr) 'FamilyID'
  )
  rhs_re <- if (length(re)>0) {
    paste(paste0('(1|', re, ')'), collapse = ' + ')
  } else {
    NULL
  }

  ## Build fixed-effects part of formula
  fe <- c(predictor, if (wconf) confounder)

  ## Construct full-model formula
  rhs_full <- paste(
    c(fe, rhs_re, 'offset(log(nf))'),
    collapse = ' + '
  )
  formula_full <- as.formula(paste0('y ~ ', rhs_full))

  ## Construct reduced formula (predictor dropped) for predictor LRT
  fe_nopred <- fe[fe!=predictor]
  rhs_nopred <- paste(
    c(if (length(fe_nopred)>0) fe_nopred, rhs_re, 'offset(log(nf))'),
    collapse = ' + '
  )
  formula_nopred <- as.formula(paste0('y ~ ', rhs_nopred))

  ## Construct reduced formula (confounder dropped) for confounder LRT
  if (wconf) {
    fe_noconf  <- fe[fe!=confounder]
    rhs_noconf <- paste(
      c(if (length(fe_noconf)>0) fe_noconf, rhs_re, 'offset(log(nf))'),
      collapse = ' + '
    )
    formula_noconf <- as.formula(paste0('y ~ ', rhs_noconf))
  }

  ## Fit full NB-GLMM
  suppressMessages(suppressWarnings({
    fit_full <- glmmTMB::glmmTMB(
      formula = formula_full,
      data    = d,
      family  = glmmTMB::nbinom2()
    )
  }))

  ## Fit reduced model (no predictor) & compute predictor LRT
  suppressMessages(suppressWarnings({
    fit_nopred <- glmmTMB::glmmTMB(
      formula = formula_nopred,
      data    = d,
      family  = glmmTMB::nbinom2()
    )
  }))
  lrt_pred <- stats::anova(fit_nopred, fit_full)
  p        <- lrt_pred$`Pr(>Chisq)`[2]
  if (is.na(p)) { p <- 1. }

  ## Extract log2FC for predictor
  coef_pred <- glmmTMB::fixef(fit_full)$cond[2]
  logfc     <- coef_pred/log(2)
  fc        <- sign(logfc)*(2^abs(logfc))

  ## Initialise confounder & interaction params
  logfc_conf <- fc_conf <- p_conf <- NA
  logfc_inter <- fc_inter <- p_inter <- NA

  if (wconf) {

    ## Fit reduced model (no confounder) & compute confounder LRT
    suppressMessages(suppressWarnings({
      fit_noconf <- glmmTMB::glmmTMB(
        formula = formula_noconf,
        data    = d,
        family  = glmmTMB::nbinom2()
      )
    }))
    lrt_conf <- stats::anova(fit_noconf, fit_full)
    p_conf   <- lrt_conf$`Pr(>Chisq)`[2]
    if (is.na(p_conf)) { p_conf <- 1. }

    ## Extract log2FC for confounder
    coef_conf  <- glmmTMB::fixef(fit_full)$cond[confounder]
    logfc_conf <- coef_conf/log(2)
    fc_conf    <- sign(logfc_conf)*(2^abs(logfc_conf))

    if (interaction) {

      ## Construct interaction formula (full model + predictor:confounder term)
      inter_term  <- paste0(predictor, ':', confounder)
      rhs_inter   <- paste(
        c(fe, inter_term, rhs_re, 'offset(log(nf))'),
        collapse = ' + '
      )
      formula_inter <- as.formula(paste0('y ~ ', rhs_inter))

      ## Fit interaction model & compute interaction LRT vs. main-effects model
      suppressMessages(suppressWarnings({
        fit_inter <- glmmTMB::glmmTMB(
          formula = formula_inter,
          data    = d,
          family  = glmmTMB::nbinom2()
        )
      }))
      lrt_inter <- stats::anova(fit_full, fit_inter)
      p_inter   <- lrt_inter$`Pr(>Chisq)`[2]
      if (is.na(p_inter)) { p_inter <- 1. }

      ## Extract log2FC for interaction term
      coef_inter  <- tail(glmmTMB::fixef(fit_inter)$cond, 1)
      logfc_inter <- coef_inter/log(2)
      fc_inter    <- sign(logfc_inter)*(2^abs(logfc_inter))
    }
  }

  ## Gather per-compartment fitted parameters
  res <- data.frame(
    'Compartment'        = comp,
    'logFC'              = logfc,
    'FC'                 = fc,
    'PValue'             = p,
    'AdjPVal'            = NA,
    'Rsq'                = NA,
    'logFCConfounder'    = logfc_conf,
    'FCConfounder'       = fc_conf,
    'PValueConfounder'   = p_conf,
    'AdjPValConfounder'  = NA,
    'logFCInteraction'   = logfc_inter,
    'FCInteraction'      = fc_inter,
    'PValueInteraction'  = p_inter,
    'AdjPValInteraction' = NA
  )

  random_intercepts <- batch_rsq <- NA
  rsq <- NA

  if (batch_aware) {

    ## Extract full batch intercepts (β₀ + û_j) on log scale
    ri <- stats::coef(fit_full)$cond$Batch[, '(Intercept)']

    ## Extract conditional SD of batch random intercepts for CI construction
    re_df  <- as.data.frame(lme4::ranef(fit_full))
    re_df  <- re_df[re_df$grpvar=='Batch', ]
    z_crit <- stats::qnorm(.975) # critical value for 2-tailed 95% CI
    ci_min <- ri-z_crit*re_df$condsd
    ci_max <- ri+z_crit*re_df$condsd

    ## Transform intercepts & CI bounds to count scale via exp()
    t_ri     <- exp(ri)
    t_ci_min <- exp(ci_min)
    t_ci_max <- exp(ci_max)

    ## Compute overall & batch-wise goodness-of-fit (R^2) from joint model
    y_pred <- predict(fit_full, newdata = d, type = 'response') # fitted counts
    sr     <- (d$y-y_pred)**2              # squared residuals
    tot    <- (d$y-mean(d$y, na.rm = TRUE))**2 # total squares
    rss    <- sum(sr, na.rm = TRUE)        # residual sum of squares
    tss    <- sum(tot, na.rm = TRUE)       # total sum of squares
    rss_batch <-
      sapply(
        batches,
        function(b) sum(sr[as.character(d$Batch)==b], na.rm = TRUE)
      ) # batch-wise RSS
    names(rss_batch) <- batches
    tss_batch <-
      sapply(
        batches,
        function(b) sum(tot[as.character(d$Batch)==b], na.rm = TRUE)
      ) # batch-wise TSS
    names(tss_batch) <- batches
    rsq       <- 1-rss/tss # overall R^2
    rsq_batch <-
      sapply(batches, function(b) 1-rss_batch[[b]]/tss_batch[[b]]) # batch R^2
    names(rsq_batch) <- batches

    res$Rsq <- rsq

    ## Gather batch random intercept parameters
    random_intercepts <- data.frame(
      'Compartment'           = rep(comp, times = nbatches),
      'Batch'                 = batches,
      'Intercept'             = t_ri,
      'ConfidenceIntervalMin' = t_ci_min,
      'ConfidenceIntervalMax' = t_ci_max
    )

    ## Gather batch-wise goodness-of-fit estimates
    batch_rsq <- data.frame(
      'Compartment' = rep(comp, times = nbatches),
      'Batch'       = batches,
      'Rsq'         = rsq_batch
    )
  }

  list(
    'res'               = res,
    'random_intercepts' = random_intercepts,
    'batch_rsq'         = batch_rsq
  )
}


## Function: fit NB-GLMM DA models across all metaclusters ----

fit_da_model <- function(
    counts,             # cell counts per metacluster and sample (MCs x samples)
    samples,            # sample names
    annotation,         # sample-level annotation
    predictor,          # biological predictor to be modelled
    confounder  = NULL, # biological confounder (at most 1)
    batch_aware = TRUE, # whether to model Batch as random intercept
    famstr      = FALSE,# whether to model FamilyID as random intercept
    interaction = FALSE,# whether to test predictor:confounder interaction
    parallel    = FALSE,# whether to use multi-threading
    verbose     = TRUE  # whether to show progress
) {

  ## Resolve confounder
  wconf <- !is.null(confounder)&&length(confounder)>0
  if (wconf&&length(confounder)!=1) {
    stop('One biological confounder at a time is currently allowed')
  }

  ## Check if an extra covariate should be taken into account
  extra_covar <- Sys.getenv('IIDX_EXTRA_COVARIATE')
  if (extra_covar=='') {
    extra_covar <- NULL
  }

  ## Set up experiment design
  experiment <- prep_experiment(
    files          = samples,
    annotation     = annotation,
    fixed_effects  = c(predictor, confounder, extra_covar),
    random_effects = c(
      if (batch_aware) 'Batch',
      if (famstr) 'FamilyID'
    ),
    force_ls = TRUE
  )
  na_annotation <- experiment[['NA']]

  ## Exclude samples with missing predictor or covariate values
  mask   <- colnames(counts)%in%na_annotation
  counts <- counts[, !mask, drop = FALSE]

  ## Establish metaclusters as compartments for testing
  comps  <- rownames(counts)
  ncomps <- length(comps)

  ## Calculate effective library sizes usign TMM
  nf <- edgeR::calcNormFactors(counts, method = 'TMM')*colSums(counts)

  ## Gather batch levels
  batches  <- if (batch_aware) levels(experiment$Data$Batch) else character(0)
  nbatches <- length(batches)

  ## Fit models per metacluster
  if (parallel) {

    ## Set up parallel processing
    cores <- detectCores()
    cl    <- makeCluster(cores[1]-1)
    registerDoSNOW(cl)

    if (verbose) {
      pb   <- utils::txtProgressBar(min = 0, max = ncomps, style = 3)
      opts <- list('progress' = function(i) utils::setTxtProgressBar(pb, i))
    } else {
      opts <- list()
    }

    out <- foreach(
      idx_comp = seq_len(ncomps),
      .combine = function(r1, r2) {
        list(
          'res'               = rbind(r1$res, r2$res),
          'random_intercepts' = rbind(r1$random_intercepts, r2$random_intercepts),
          'batch_rsq'         = rbind(r1$batch_rsq, r2$batch_rsq)
        )
      },
      .inorder      = TRUE,
      .options.snow = opts,
      .export       = c('da_singlefit'),
      .packages     = c('glmmTMB', 'lme4', 'edgeR')
    ) %dopar% {
      da_singlefit(
        comps, idx_comp, counts, nf, experiment,
        predictor, confounder, batch_aware, famstr,
        batches, nbatches, wconf, interaction
      )
    }
    if (verbose) { close(pb) }
    stopCluster(cl)

    res       <- out$res
    rand_int  <- out$random_intercepts
    batch_rsq <- out$batch_rsq

  } else { # sequential

    if (verbose) {
      pb <- utils::txtProgressBar(min = 0, max = ncomps, style = 3)
    }

    out <- lapply(
      seq_len(ncomps),
      function(idx_comp) {
        if (verbose) { utils::setTxtProgressBar(pb, idx_comp) }
        da_singlefit(
          comps, idx_comp, counts, nf, experiment,
          predictor, confounder, batch_aware, famstr,
          batches, nbatches, wconf, interaction
        )
      }
    )
    if (verbose) { close(pb) }

    res       <- do.call(rbind, lapply(out, function(x) x$res))
    rand_int  <- do.call(rbind, lapply(out, function(x) x$random_intercepts))
    batch_rsq <- do.call(rbind, lapply(out, function(x) x$batch_rsq))
  }

  rownames(res) <- NULL
  if (is.data.frame(rand_int))  { rownames(rand_int)  <- NULL }
  if (is.data.frame(batch_rsq)) { rownames(batch_rsq) <- NULL }

  ## Apply BH correction for multiple testing across metaclusters
  res$AdjPVal <- stats::p.adjust(res$PValue, method = 'BH')
  if (wconf) {
    res$AdjPValConfounder <-
      stats::p.adjust(res$PValueConfounder, method = 'BH')
    if (interaction) {
      res$AdjPValInteraction <-
        stats::p.adjust(res$PValueInteraction, method = 'BH')
    }
  }

  ## Add metadata
  attributes(res)$Continuous <- is.numeric(experiment$Data[, predictor])
  attributes(res)$ContinuousConfounder <-
    if (wconf) is.numeric(experiment$Data[, confounder]) else NULL
  attributes(res)$na_annotation       <- na_annotation
  attributes(res)$confounder_specified <- wconf

  list(
    'main'              = res,
    'random_intercepts' = rand_int,
    'batch_r_squared'   = batch_rsq,
    'na_annotation'     = na_annotation
  )
}
