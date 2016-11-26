#Enhancing gene regulatory network inference through data integration with markov random fields
<br />
Home of the **GRACE** (**G**ene **R**egulatory network inference **AC**curacy **E**nhancement) algorithm <br />

<br />
<br />
*Summary* <br />
A gene regulatory network links transcription factors to their target genes and represents a map of transcriptional regulation. Much progress has been made in deciphering gene regulatory networks computationally. However, gene regulatory network inference for most eukaryotic organisms remain challenging. To improve the accuracy of gene regulatory network inference and facilitate candidate selection for experimentation, we developed an algorithm called GRACE (Gene Regulatory network inference ACcuracy Enhancement). GRACE exploits biological *a priori* and heterogeneous data integration to generate high-confidence network predictions for eukaryotic organisms using Markov Random Fields in a semi-supervised fashion. GRACE uses a novel optimization scheme to integrate regulatory evidence and biological relevance. It is particularly suited for model learning with sparse regulatory gold standard data. We show GRACE's potential to produce high confidence regulatory networks compared to state of the art approaches using *Drosophila melanogaster* and *Arabidopsis thaliana* data. In an *A. thaliana* developmental gene regulatory network, GRACE recovers cell cycle related regulatory mechanisms and further hypothesizes several novel regulatory links, including a putative control mechanism of vascular structure formation due to modifications in cell proliferation.
<br />

![Alt text](/inference_novel.jpg?raw=true "GRACE workflow")

*General usage* <br />
Step 1) Initial gene regulatory network inference (Figure 1 A) is based on fast random forest regression followed by DNA binding prediction map filtering (GRACE_initial_GRN_inference.R) <br/>
<br/>
*Specific parameters and datasets* <br/>
th.grnCutoff <- 0.999 // initial cutoff for random forest regression on gene expression network (should result in not more than 20000 - 30000 links) <br/>
n.cpus <- 2 // number of available cpus for parallel random forest regression <br/>
m.expression <- ... // load expression matrix for grn inference - matrix - rows (genes) x cols (conditions) <br/>
v.tfs <- ... // character vector of transcription factor genes <br/>
df.dna_binding <- // load DNA binding set (dataframe with 2 columns) <br/>
<br/>
Step 2) Configure GRACE_load_datasets.R to set paths to the location of all needed datasets.<br/>
<br/>
<br/>
Step 3) Configure GRACE_pipeline_template.R, as described within the script, to run GRACE on the initial network using cofunctional network information and regulatory as well as cofunctional evidence data as training sets (Figure 1 B-C)<br/>
<br/>
*Specific parameters and datasets* <br/>
n.cpus <- 2 <br/>
beta <- 1 # equal emphasis on precision (atrm) and recall (bp coreg pairs) - introducing (novel) combined f-measure (physical + functional evidence, singularity handling, minimum size) <br/>
b.normalize_precision <- TRUE <br/>
b.jaccard = TRUE <br/>
n.sets <- 20 <br/>
lambda.gridSearch <- c(0.01,seq(0.5,2.5,0.5)) <br/>
th.percentage.hyperparameters = 0.99 # grid search <br/>
max.call =  200 # simulated annealing <br/>
n.models <- 100 # 100 models <br/>
n.sample_size <- 0.632 # as in traditional bootrapping approaches <br/>
<br/>
*Tutorials*<br />
In addition, we have prepared two standalone tutorials (for A. thaliana and D. melanogaster) in order to reproduce results <br/>
GRACE_tutorial_athaliana.R - represents a less optimized version of GRACE - can be used to recompute the predictions as given in the Supplement <br/> 

GRACE_tutorial_dmelanogaster.R - based on the current version of GRACE - recompute the published predictions as well comparision with

Reference for the datasets used within the D. melanogaster model:
Marbach D, Roy S, Ay F, Meyer PE, Candeias R, Kahveci T, Bristow CA, Kellis M. Predictive regulatory models in Drosophila melanogaster by integrative inference of transcriptional networks. Genome Res. 2012 Jul;22(7):1334-49.
<br />

By default, the precomputed GRACE models are loaded for further processing, otherwise they can be recomputed. <br/>

<br />
For help or questions please contact: <br />
mbanf.research(at)gmail.com

<br />

