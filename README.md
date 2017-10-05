# agedbrain
This repository contains code and files for reproducing the analysis presented in "Neuropathological and transcriptomic characteristics of the aged brain," [LINK/INFO HERE] the platform paper for the Aging, Dementia, and TBI Study. 

To reproduce the analysis, do the following:
1) Download everything from the Aging, Dementia, and TBI Study download page (http://aging.brain-map.org/download/index) and save it in the same folder
2) Download this github repository files.  Unzip "dataFiles.zip" and save to a folder.  Save "TbT_normalization.RData" in the same folder
3) Save the six .r files from this github repository in a different folder
4) Update the top several lines of "Code01_ReadInDataAndFormat.r" with the appropriate folder names (from above)
5) Install the following R packages: sva, WGCNA, gridExtra, ggplot2, gplots, baySeq, edgeR, DESeq, NBPSeq, ROC, and mclust
6) Run the five numbered code files in order

Some useful websites:
1) Allen Institute data portal: http://brain-map.org/
2) The Aging, Dementia, and TBI Study main page: http://aging.brain-map.org/
3) Raw data (FASTQ/bam files) for the RNA-Seq data presented in this study: https://www.niagads.org/datasets/ng00059
4) Link to the open-access eLife article describing this resource: [Not yet available]

If you find any errors in this code, please let me know so I can correct them.  I have validated this code in R version 3.2.5 (2016-04-14), and don't expect other errors, but you never know.
