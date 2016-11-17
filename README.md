#Enhancing gene regulatory network inference through data integration with markov random fields
<br />
Home of the **GRACE** (**G**ene **R**egulatory network inference **AC**curacy **E**nhancement) algorithm <br />

<br />
<br />
*Summary* <br />
A gene regulatory network links transcription factors to their target genes and represents a map of transcriptional regulation. Much progress has been made in deciphering gene regulatory networks computationally. However, gene regulatory network inference for most eukaryotic organisms remain challenging. To improve the accuracy of gene regulatory network inference and facilitate candidate selection for experimentation, we developed an algorithm called GRACE (Gene Regulatory network inference ACcuracy Enhancement). GRACE exploits biological *a priori* and heterogeneous data integration to generate high-confidence network predictions for eukaryotic organisms using Markov Random Fields in a semi-supervised fashion. GRACE uses a novel optimization scheme to integrate regulatory evidence and biological relevance. It is particularly suited for model learning with sparse regulatory gold standard data. We show GRACE's potential to produce high confidence regulatory networks compared to state of the art approaches using *Drosophila melanogaster* and *Arabidopsis thaliana* data. In an *A. thaliana* developmental gene regulatory network, GRACE recovers cell cycle related regulatory mechanisms and further hypothesizes several novel regulatory links, including a putative control mechanism of vascular structure formation due to modifications in cell proliferation.
<br />

![Alt text](/inference_novel.jpg?raw=true "GRACE workflow" style="width=10px;"}

*General usage* <br />
Initial gene regulatory network inference (figure 1 A) is based on fast random forest regression followed by DNA binding prediction map filtering (see below). To construct the initial network follow the steps in GRACE_initial_GRN_inference.R. 

Second, configure GRACE_load_datasets.R to set paths to the location of all needed datasets.

Subsequently, the GRACE_pipeline_template.R script can be configured, as described within the script, to run GRACE on the initial network using cofunctional network information and regulatory as well as cofunctional evidence data as training sets. 

*Tutorials*<br />
In addition, we have prepared two standalone tutorials (for A. thaliana and D. melanogaster) in order to reproduce results, (GRACE_tutorial_athaliana.R and GRACE_tutorial_dmelanogaster.R). By default, the precomputed GRACE models are loaded for further processing, otherwise they can be recomputed.

<br />
*DNA binding constraint random forest regression based gene regulatory network inference*

![Alt text](/inference_randomForest.jpg?raw=true "random forest regression based network inference")



<br />
For help or question please contact: <br />
michael.banf(at)gmx.net

<br />

