# Genomic prediction within and across maize landrace derived populations using haplotypes

## Description
See [DataDescription](https://rawcdn.githack.com/yan-cheng-lin/HBGP/646f3b62edd013a573a3d17bfbaa9471568badf7/Rscripts/00-01_DataDescription.html) for the description of data used in the analysis and [TablesAndFigures](https://rawcdn.githack.com/yan-cheng-lin/HBGP/646f3b62edd013a573a3d17bfbaa9471568badf7/Rscripts/03-01_TablesAndFigures.html) for results description and tables & figures reproduction of the paper. All the haplotype construction and genomic prediction in the paper can be reproduced via scripts in [Rscripts](https://github.com/yan-cheng-lin/HBGP/tree/main/Rscripts) and with the [Data](https://github.com/yan-cheng-lin/HBGP/tree/main/Data) \(please contact authors for some large files\).  

**00-01_DataDescription.Rmd:** Rmarkdown for Description of the genotypic and phenotypic data.  
**00-02_CVSamples.R:** Sampling for cross-validation of genomic predcition.  
**01-01_FixedHB_construction.R:** RScript for *FixedHB* haplotype construction.  
**01-02_HaploView_construction.R:** RScript for *HaploView* haplotype construction. Should be implemented with shellscript **01-02_Run_HaploView.sh**.  
**01-03_HaploBlocker_construction.R:** RScript for *HaploBlocker* haplotype construction.  
**02-02_CrossValidation_WeightedGBLUP.R** RScript for demonstration of weighting the haplotype length in genomic relationship matrix.  
**02-01_CrossValidation_GBLUP.R:** RScript for haplotype-based genomic prediction for all three scenario.  
**03-01_TablesAndFigures.Rmd:** Rmarkdown for source code of the tables and figures.  

## Authors
Yan-Cheng Lin

## Contact
yancheng.lin@tum.de
