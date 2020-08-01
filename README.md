## Metatranscriptome survey, station ALOHA & coastal California 

Data analysis for metatranscriptome survey from station ALOHA and coastal California  

* Hu et al. (_in prep_). Distinct transcriptional signatures of the protistan community at the SPOT and ALOHA ocean time-series stations characterize different depths and locations.

## Links to available data
* [Zenodo](https://zenodo.org/record/3954884#.XyXWNi-ZNTY)
* [SRA][https://github.com/shu251/CA-ALOHA-metatranscriptome-analysis/blob/master/download-raw-metaT-fastq.sh)

## Contents of repo
* ```compile-raw-data/``` source code to compile raw count, taxonomic assignment, and Kegg identities. Generates raw count table, normalized data sets, and associated  R objects, these are available through linked Zenodo.
* ```environmental-parameters/``` plot vertical profiles from station ALOHA, the SPOT station, Catalina Island, and the Port of LA
* ```metatranscriptome-analysis-CA-ALOHA-SKH.ipynb``` - Jupyter notebook with complete R script to reproduce all figures and analyses. Flat .R script and statis .html script also included.
* _.RData_ - R data files to run in provided R code
* _Custom_KO_list_30-04-2020.txt_ - text file with list of genes and pathways of interest for study
* _taxonomic-assign-reference.txt_ - text file listing taxonomic assignments and curation of taxonomic placement

https://github.com/shu251/CA-ALOHA-metatranscriptome-analysis/blob/master/metatranscriptome-analysis-CA-ALOHA-SKH.html




