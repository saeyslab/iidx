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

## iidx internal module: 06c_InfoText.R
## github.com/saeyslab/iidx

# Maintainer: David Novak (davidnovak9000@gmail.com)
# Description: Defines a modal dialog window within the dashboard with iidx
# workflow information.

## Define info text modal ----

info_modal <- shiny::modalDialog(
  title     = NULL,
  size      = 'l',
  easyClose = TRUE,
  footer    = NULL,
  shiny::HTML("
<h3>iidx <b>differential expression workflow</b></h3>
<hr style='margin-top:1em; margin-bottom: 1em;' />
  
This dashboard presents the results of a differential expression analysis of a 
large, multi-batch cytometry dataset.<br><br>

We draw from and expand on the <b>
<a href='https://www.nature.com/articles/s41596-021-00550-0' target='_blank' rel='noopener noreferrer'><i>FlowSOM</i></a>
clustering protocol</b> and the <b>
<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6517415/' target='_blank' rel='noopener noreferrer'><i>diffcyt</i></a> 
differential expression (DE) analysis framework</b> to provide an advanced 
workflow for computational and biological scientists to collaborate in data 
exploration. Our solution focuses on 
<ul>
  <li>mitigating and explaining <b>technical batch effects</b></li>
  <li><b>disentangling interactions</b> between biological predictors of phenotype</li>
  <li>allowing domain experts to <b>explore results interactively</b></li>
</ul>

Key features of our workflow are summarised below.<br><br>

<h4><b>Marker function and stability</b></h4>

We distinguish between 
<ul>
  <li><b>lineage markers</b>: assumed stable, define cell identities</li>
  <li><b>state markers</b>: assumed transient, define cell states</li>
</ul>

Bimodally expressed markers are assigned <b>phenopositivity thresholds</b> for 
distinction between negative and positive populations.<br><br>

If serious batch effect is detected in lineage markers, we correct it via <b>
<a href='https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.23904' target='_blank' rel='noopener noreferrer'><i>CytoNorm</i></a> 
signal normalisation</b>, to ensure stable cell phenotyping across batches. 
Any <b>remaining batch effect in state markers is tackled in statistical 
modelling</b> downstream.<br><br>

<h4><b>Feature extraction and outlier detection</b></h4>

We apply <i>FlowSOM</i> metaclustering to a large aggregate of expression data 
and map each donor sample to the fitted model to extract phenotypic features:
<ul>
  <li><b>abundances</b>: sizes and proportions of each metacluster per donor</li>
  <li><b>average state marker expression</b>: median fluorescence/signal intensities (MFIs) of per state marker per metacluster per donor</li>
  <li><b>phenopositivity rates</b>: proportions of phenopositive cells per bimodally expressed state marker per metacluster per donor</li>
</ul>

Samples with extreme outlier values for metacluster abundances are flagged for 
removal as outliers.<br><br>

<h4><b>Differential expression modelling</b></h4>

In line with the <i>diffcyt</i> framework, we model changes in marker 
expression in terms of <b>differential abundance (DA)</b> and <b>differential 
state (DS)</b>. Both DA and DS seek to model phenotypic changes based on 
sample-level attributes (age, sex, disease status, ...). <b>DA describes 
relative over- or under-representation of cell compartments</b> (FlowSOM 
metaclusters). <b>DS describes shifts in expression of state markers</b> within 
these compartments.<br><br>

The exploratory DE analysis facilitates eventual post-hoc testing of hypotheses 
for combined changes in multiple compartments or markers at a single-donor or 
single-cell level.<br><br>

We describe our DA and DS methodologies below.<br><br>

<h5>&mdash; <b>Differential abundance</b></h5>

<p style='margin-left: 20px;'>
The presented DA approach uses an adaptation of <b>
<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/' target='_blank' rel='noopener noreferrer'><i>edgeR</i></a>
negative-binomial generalised linear models (GLMs)</b> with 
<a href='https://projecteuclid.org/journals/annals-of-applied-statistics/volume-10/issue-2/Robust-hyperparameter-estimation-protects-against-hypervariable-genes-and-improves-power/10.1214/16-AOAS920.full' target='_blank' rel='noopener noreferrer'>empirical Bayes estimation of dispersion</a> 
and
<a href='https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25'>TMM</a> 
normalisation of sample sizes, due to the robustness of this approach. 
<b>Predictor, covariate/confounder (if used) and sample batch are taken as fixed 
effects</b> in a multivariate model fitted to data from all batches. Parameters 
of models fitted per single batch are extracted to monitor potential 
inconsistencies in effect size and significance across batches (due to 
remaining technical batch effect or sampling bias).<br><br>
</p>

<h5>&mdash; <b>Differential state</b></h5>

<p style='margin-left: 20px;'>
To test DS, we consider changes in both state marker MFI (<b>DS-MFI</b>) and their 
respective phenopositivity rates (<b>DS-Pheno</b>). Since batch effect in state 
marker signal is assumed to be more persistent than that in lineage markers, 
<b>sample batch is modelled via random intercepts</b> in both analyses. This 
corrects for spurious effects driven by technical artifacts.<br><br>

DS-MFI uses <b>
<a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402510/' target='_blank' rel='noopener noreferrer'><i>limma</i></a>
linear mixed models (LMMs)</b> and DS-Pheno uses <b>Beta-generalised linear 
mixed models (Beta-GLMMs)</b>. To interpret the stability of observed effects,
random intercepts per batch (transformed to the original feature space) are 
displayed along with corresponding estimates of 95% confidence intervals. 
Additionally, batch-level R<sup>2</sup>-like goodness-of-fit estimates are 
computed and presented to show the influence of individual batches and reveal 
potential inconsistencies.<br><br>

Just as DA, DS-MFI and DS-Pheno allow for the <b>inclusion of a confounder 
variable to disentangle combined effects</b> using a multivariate model.<br>
</p><br>

<h4><b>Filtering by robustness w.r.t. background signal</b></h4>

If FCS data is available for unstained samples (<i>i.e.</i>, without 
antibody-conjugated fluorophore or heavy metal tags), we determine the expected 
levels of background signal due to autofluorescence or imperfect instrument 
calibration, per sample batch. This makes it possible to <b>filter differential 
state results by robustness with respect to background noise</b>. Each hit (<i>
i.e.</i>, test result for specific metacluster, marker, predictor and, 
optionally, confounder) is filtered by the proportion of samples with MFI that 
clears a threshold of signal considered to exceed the background noise. Such 
samples are considered to have biologically meaningful signal.<br><br>

The stringency of this filtering can be adjusted in the dashboard in two ways:
<ul>
  <li><b>signal filter strength</b>: the number of median abolute deviations (MADs) above background signal median considered as the robustness threshold for marker signal</li>
  <li><b>sample filter strength</b>: the minimum proportion of samples in a compartment that need to clear the signal threshold for the hit to be considered robust</li>
</ul>
  ")
)
