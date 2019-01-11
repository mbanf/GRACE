# Enhancing gene regulatory network inference through data integration with markov random fields


[![GRACE](http://img.youtube.com/vi/KOxbO0EI4MA/0.jpg)](https://www.youtube.com/watch?v=sQ02p4CU6UM "GRACE")


## Notes
---------- Important: GRACE is currently down for maintenence - will be up again shortly  --------


<br>

[Contact the author](mailto:michael@junkdna.ai) for help or questions. 


## About
A gene regulatory network links transcription factors to their target genes and represents a map of transcriptional regulation. Much progress has been made in deciphering gene regulatory networks computationally. However, gene regulatory network inference for most eukaryotic organisms remain challenging. To improve the accuracy of gene regulatory network inference and facilitate candidate selection for experimentation, we developed an algorithm called GRACE (Gene Regulatory network inference ACcuracy Enhancement). GRACE exploits biological *a priori* and heterogeneous data integration to generate high-confidence network predictions for eukaryotic organisms using Markov Random Fields in a semi-supervised fashion. GRACE uses a novel optimization scheme to integrate regulatory evidence and biological relevance. It is particularly suited for model learning with sparse regulatory gold standard data. We show GRACE's potential to produce high confidence regulatory networks compared to state of the art approaches using *Drosophila melanogaster* and *Arabidopsis thaliana* data. In an *A. thaliana* developmental gene regulatory network, GRACE recovers cell cycle related regulatory mechanisms and further hypothesizes several novel regulatory links, including a putative control mechanism of vascular structure formation due to modifications in cell proliferation.
<br />

![Alt text](/inference_novel.jpg?raw=true "GRACE workflow")

*General usage* <br />
Step 1) Initial gene regulatory network inference (Figure 1 A) is based on fast random forest regression followed by DNA binding prediction map filtering (GRACE_initial_GRN_inference.R) <br/>
<br/>
Step 2) Configure GRACE_load_datasets.R to set paths to the location of all needed datasets.<br/>
<br/>
Step 3) Configure GRACE_pipeline_template.R, as described within the script, to run GRACE on the initial network using cofunctional network information and regulatory as well as cofunctional evidence data as training sets (Figure 1 B-C)<br/>
<br/>
<br/>
*Tutorials*<br />

Get the original code and datasets from onedrive: [GRACE original code and datasets](https://1drv.ms/u/s!Avm82Xhe9EZj411IGorBFt8zpwKp). 

In addition, we have prepared two standalone tutorials (for A. thaliana and D. melanogaster) in order to reproduce results <br/>
GRACE_tutorial_athaliana.R - represents a less optimized version of GRACE - can be used to recompute the prioritized predictions (unpack A. thaliana dataset in GRACE folder) <br/> 

GRACE_tutorial_dmelanogaster.R - based on the current version of GRACE - recompute the prioritized predictions (unpack D. melanogaster dataset in GRACE folder) <br/>

By default, the precomputed GRACE models are loaded for further processing, otherwise they can be recomputed. <br/>

# References



GRACE algorithm and results<br />
Banf M, and Rhee S. Enhancing gene regulatory network inference through data integration with markov random fields. Scientific Reports,7:41174, Nature publishing group. (http://www.nature.com/articles/srep41174)
<br />
Datasets used within the D. melanogaster model<br />
Marbach D, Roy S, Ay F, Meyer PE, Candeias R, Kahveci T, Bristow CA, Kellis M. Predictive regulatory models in Drosophila melanogaster by integrative inference of transcriptional networks. Genome Res. 2012 Jul;22(7):1334-49.

<br />

