
input_files_needed <- FALSE
source('99_SetUp.R')
source('04b_FullExperiments.R')
library(tidyverse)

Sys.setenv(IIDX_NCORES = '3')

## Full w/ & w/out inter ----

inputs <- prep_inputs(fname_features_postnorm, fname_sample_outliers, annotation, batch_remove)
Sys.setenv(IIDX_EXTRA_COVARIATE = '')

message('B DA regular')
res_da <- test_da(counts = inputs$counts, annotation = annotation, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Full_Regular.RDS'); rm(res_da); gc(); gc(); gc()
message('B DA interactions')
res_da <- test_da(counts = inputs$counts, annotation = annotation, predictors = predictors, confounders = confounders, interactions = TRUE, batch_aware = FALSE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Full_Interactions.RDS'); rm(res_da); gc(); gc(); gc()

message('B DS-MFI regular')
res_ds_mfi <- test_ds(annotation = annotation, mfi = inputs$MFIs, counts = inputs$counts, predictors  = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE, parallel = TRUE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Regular.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
message('B DS-MFI interactions')
res_ds_mfi <- test_ds(annotation = annotation, mfi = inputs$MFIs, counts = inputs$counts, predictors = predictors, confounders = confounders, interactions = TRUE, batch_aware = FALSE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE, parallel = TRUE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Interactions.RDS'); rm(res_ds_mfi); gc(); gc(); gc()

message('B DS-Pheno regular')
res_ds_pheno <- test_ds(annotation = annotation, phenopos = inputs$PercPos, counts = inputs$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], parallel = TRUE, verbose = TRUE)
saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Regular.RDS'); rm(res_ds_pheno); gc(); gc(); gc()
message('B DS-Pheno interactions')
res_ds_pheno <- test_ds(annotation = annotation, phenopos = inputs$PercPos, counts = inputs$counts, predictors = predictors, confounders = confounders, interactions = TRUE, batch_aware = FALSE, state_markers = as.vector(channels)[idcs_channels_state], verbose = TRUE, parallel = TRUE)
saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Interactions.RDS'); rm(res_ds_pheno); gc(); gc(); gc()

## By cohort ----

inputs <- prep_inputs(fname_features_postnorm, fname_sample_outliers, annotation, batch_remove)

mp <- match(rownames(inputs$counts), annotation$FileName)
coh <- annotation$Cohort[mp]
fid <- annotation$FamilyID[mp]
mask_vrc <- coh=='VRC'
mask_twins <- coh=='UK' & (fid%in%unique(fid[duplicated(fid)]))
mask_twin1 <- mask_twins & duplicated(fid)
mask_twin2 <- mask_twins & !duplicated(fid)

inputs_vrc             <- inputs
inputs_vrc$counts      <- inputs$counts[     mask_vrc, ]
inputs_vrc$percentages <- inputs$percentages[mask_vrc, ]
inputs_vrc$MFIs        <- inputs$MFIs[       mask_vrc, ]
inputs_vrc$PercPos     <- inputs$PercPos[    mask_vrc, ]

inputs_uk             <- inputs
inputs_uk$counts      <- inputs$counts[     mask_twins, ]
inputs_uk$percentages <- inputs$percentages[mask_twins, ]
inputs_uk$MFIs        <- inputs$MFIs[       mask_twins, ]
inputs_uk$PercPos     <- inputs$PercPos[    mask_twins, ]

inputs_twin1             <- inputs
inputs_twin1$counts      <- inputs$counts[     mask_twin1, ]
inputs_twin1$percentages <- inputs$percentages[mask_twin1, ]
inputs_twin1$MFIs        <- inputs$MFIs[       mask_twin1, ]
inputs_twin1$PercPos     <- inputs$PercPos[    mask_twin1, ]

inputs_twin2             <- inputs
inputs_twin2$counts      <- inputs$counts[     mask_twin2, ]
inputs_twin2$percentages <- inputs$percentages[mask_twin2, ]
inputs_twin2$MFIs        <- inputs$MFIs[       mask_twin2, ]
inputs_twin2$PercPos     <- inputs$PercPos[    mask_twin2, ]

inputs_vrc_twin1             <- inputs
inputs_vrc_twin1$counts      <- inputs$counts[     mask_vrc|mask_twin1, ]
inputs_vrc_twin1$percentages <- inputs$percentages[mask_vrc|mask_twin1, ]
inputs_vrc_twin1$MFIs        <- inputs$MFIs[       mask_vrc|mask_twin1, ]
inputs_vrc_twin1$PercPos     <- inputs$PercPos[    mask_vrc|mask_twin1, ]

inputs_vrc_twin2             <- inputs
inputs_vrc_twin2$counts      <- inputs$counts[     mask_vrc|mask_twin2, ]
inputs_vrc_twin2$percentages <- inputs$percentages[mask_vrc|mask_twin2, ]
inputs_vrc_twin2$MFIs        <- inputs$MFIs[       mask_vrc|mask_twin2, ]
inputs_vrc_twin2$PercPos     <- inputs$PercPos[    mask_vrc|mask_twin2, ]

ann_woutfam <- annotation[, colnames(annotation)[colnames(annotation)!='FamilyID']]

message('B DA cohort: VRC')
res_da <- test_da(counts = inputs_vrc$counts, annotation = ann_woutfam, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Cohort_VRC.RDS'); rm(res_da); gc(); gc(); gc()
message('B DA cohort: UK')
res_da <- test_da(counts = inputs_uk$counts, annotation = annotation, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Cohort_UK.RDS'); rm(res_da); gc(); gc(); gc()
message('B DA cohort: Twin1')
res_da <- test_da(counts = inputs_twin1$counts, annotation = ann_woutfam, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Cohort_Twin1.RDS'); rm(res_da); gc(); gc(); gc()
message('B DA cohort: Twin2')
res_da <- test_da(counts = inputs_twin2$counts, annotation = ann_woutfam, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Cohort_Twin2.RDS'); rm(res_da); gc(); gc(); gc()
message('B DA cohort: VRC+Twin1')
res_da <- test_da(counts = inputs_vrc_twin1$counts, annotation = ann_woutfam, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Cohort_VRC+Twin1.RDS'); rm(res_da); gc(); gc(); gc()
message('B DA cohort: VRC+Twin2')
res_da <- test_da(counts = inputs_vrc_twin2$counts, annotation = ann_woutfam, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_Cohort_VRC+Twin2.RDS'); rm(res_da); gc(); gc(); gc()

message('B DS-MFI cohort: VRC')
res_ds_mfi <- test_ds(annotation = ann_woutfam, mfi = inputs_vrc$MFIs, counts = inputs_vrc$counts, predictors= predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Cohort_VRC.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
message('B DS-MFI cohort: UK')
res_ds_mfi <- test_ds(annotation = annotation, mfi = inputs_uk$MFIs, counts = inputs_uk$counts, predictors= predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Cohort_UK.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
message('B DS-MFI cohort: Twin1')
res_ds_mfi <- test_ds(annotation = ann_woutfam, mfi = inputs_twin1$MFIs, counts = inputs_twin1$counts, predictors= predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Cohort_Twin1.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
message('B DS-MFI cohort: Twin2')
res_ds_mfi <- test_ds(annotation = ann_woutfam, mfi = inputs_twin2$MFIs, counts = inputs_twin2$counts, predictors= predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Cohort_Twin2.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
message('B DS-MFI cohort: VRC+Twin1')
res_ds_mfi <- test_ds(annotation = ann_woutfam, mfi = inputs_vrc_twin1$MFIs, counts = inputs_vrc_twin1$counts, predictors= predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Cohort_VRC+Twin1.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
message('B DS-MFI cohort: VRC+Twin2')
res_ds_mfi <- test_ds(annotation = ann_woutfam, mfi = inputs_vrc_twin2$MFIs, counts = inputs_vrc_twin2$counts, predictors= predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_Cohort_VRC+Twin2.RDS'); rm(res_ds_mfi); gc(); gc(); gc()

# message('B DS-Pheno cohort: VRC')
# res_ds_pheno <- test_ds(annotation = ann_woutfam, phenopos = inputs_vrc$PercPos, counts = inputs_vrc$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Cohort_VRC.RDS'); rm(res_ds_pheno); gc(); gc(); gc()
# message('B DS-Pheno cohort: UK')
# res_ds_pheno <- test_ds(annotation = ann_woutfam, phenopos = inputs_uk$PercPos, counts = inputs_uk$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Cohort_UK.RDS'); rm(res_ds_pheno); gc(); gc(); gc()
# message('B DS-Pheno cohort: Twin1')
# res_ds_pheno <- test_ds(annotation = ann_woutfam, phenopos = inputs_twin1$PercPos, counts = inputs_twin1$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Cohort_Twin1.RDS'); rm(res_ds_pheno); gc(); gc(); gc()
# message('B DS-Pheno cohort: Twin2')
# res_ds_pheno <- test_ds(annotation = ann_woutfam, phenopos = inputs_twin2$PercPos, counts = inputs_twin2$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Cohort_Twin2.RDS'); rm(res_ds_pheno); gc(); gc(); gc()
# message('B DS-Pheno cohort: VRC+Twin1')
# res_ds_pheno <- test_ds(annotation = ann_woutfam, phenopos = inputs_vrc_twin1$PercPos, counts = inputs_vrc_twin1$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Cohort_VRC+Twin1.RDS'); rm(res_ds_pheno); gc(); gc(); gc()
# message('B DS-Pheno cohort: VRC+Twin2')
# res_ds_pheno <- test_ds(annotation = ann_woutfam, phenopos = inputs_vrc_twin2$PercPos, counts = inputs_vrc_twin2$counts, predictors = predictors, confounders = confounders, interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_Cohort_VRC+Twin2.RDS'); rm(res_ds_pheno); gc(); gc(); gc()

## CMV ----


message('B DA CMV')
Sys.setenv(IIDX_EXTRA_COVARIATE = 'Sex'); res_da <- test_da(counts = inputs$counts, annotation = annotation, predictors = c('CMV', 'Age'), confounders = c('CMV', 'Age'), interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_CMV_Age.RDS'); rm(res_da); gc(); gc(); gc()
Sys.setenv(IIDX_EXTRA_COVARIATE = 'Age'); res_da <- test_da(counts = inputs$counts, annotation = annotation, predictors = c('CMV', 'Sex'), confounders = c('CMV', 'Sex'), interactions = FALSE, batch_aware = TRUE, verbose = FALSE, parallel = TRUE)
saveRDS(res_da, '../NEW RESULTS/B_DA_CMV_Sex.RDS'); rm(res_da); gc(); gc(); gc()
message('B DS-MFI CMV')
Sys.setenv(IIDX_EXTRA_COVARIATE = 'Sex'); res_ds_mfi <- test_ds(annotation = annotation, mfi = inputs$MFIs, counts = inputs$counts, predictors = c('CMV', 'Age'), confounders = c('CMV', 'Age'), interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_CMV_Age.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
Sys.setenv(IIDX_EXTRA_COVARIATE = 'Age'); res_ds_mfi <- test_ds(annotation = annotation, mfi = inputs$MFIs, counts = inputs$counts, predictors = c('CMV', 'Sex'), confounders = c('CMV', 'Sex'), interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
saveRDS(res_ds_mfi, '../NEW RESULTS/B_DS-MFI_CMV_Sex.RDS'); rm(res_ds_mfi); gc(); gc(); gc()
# message('B DS-Pheno CMV')
# res_ds_pheno <- test_ds(annotation = annotation, phenopos = inputs$PercPos, counts = inputs$counts, predictors = c('CMV', 'Age'), confounders = c('CMV', 'Age'), interactions = FALSE, batch_aware = TRUE, state_markers = as.vector(channels)[idcs_channels_state], verbose = FALSE)
# saveRDS(res_ds_pheno, '../NEW RESULTS/B_DS-Pheno_CMV.RDS'); rm(res_ds_pheno); gc(); gc(); gc()

message('Done')









