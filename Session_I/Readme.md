# Materials for session 1
Contact: Sambhawa Priya (*priya030 at umn.edu*)

This contains script and data for Session 1 tutorial on applying machine learning to microbiome data.

### Directory structure ### 
Please create a folder named _Session_1_tutorial_ at a relevant location on your computer.  
Place the Rscript _MLMB_session_1_tutorial.R_ in _Session_1_tutorial_.  
Next, create a directory within _Session_1_tutorial_ called _input_ and place the input files described below (metadata and otu table) in the _input_ directory.  

### Input data (available under folder "input") ###

This data was originally published by [Singh, Pallavi et al. “Intestinal microbial communities associated with acute enteric infections and disease recovery.” Microbiome vol. 3 45. 22 Sep. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4579588/). It was complied as part of _MicrobiomeHD_ as described in [Duvallet et al. 2017](https://doi.org/10.1038/s41467-017-01973-8), and further analyzed and described by [Zhou and Gallins et al. 2019](https://www.frontiersin.org/article/10.3389/fgene.2019.00579). 

- **edd_singh.metadata.txt**: The metadata file containing description of samples in the dataset.
- **edd_singh_otu_table**: The OTU table containing microbiome abundance data (16S amplicon sequencing).

### R code ###

**MLMB_session_1_tutorial.R**: R script to implement machine learning pipeline on microbiome data.
   
