#################################################################################################
##  These scripts take data from the the Aging, Dementia, and TBI website along with supplementary 
##  materials and uses it to perform all analyses for the aging and dementia manuscript 
#################################################################################################

# R version 3.2.5 (2016-04-14) -- "Very, Very Secure Dishes" was used for the analysis, but the code should
# be compatible with most of the recent versions of R.

#################################################################################################
print("***** -------------------------------------------")
print("***** Code #1: Read in all data for the analysis.")
print("***** -------------------------------------------")

# THESE MUST BE SPECIFIED IN ADVANCE

# LOCATION OF MAIN FOLDER WITH ALL OF THE FILES
mainFolder    = "\\\\allen/programs/celltypes/workgroups/hct/TBI_Analysis/manuscript/dataAnalysis/"

# LOCATION OF ALL OUTPUT FILES
outputFolder  = paste(mainFolder,"outputFiles_NN/",sep="")     

# LOCATION OF ALL FILES DOWNLOADED FROM THE WEBSITE (http://aging.brain-map.org/download/index)
inputFolder   = paste(mainFolder,"downloads/",sep="")       

# LOCATION OF THE DATA FILES UNZIPPED FROM "dataFiles.zip" (DOWNLOADED FROM GITHUB)
extraFolder   = paste(mainFolder,"additionalFiles/",sep="")

# LOCATION OF THE SCRIPTS DOWNLOADED FROM GITHUB (including this one)
scriptsFolder = paste(mainFolder,"code/",sep="")

#################################################################################################
print("The following libraries are required to run this code: SVA")

library(sva)
library(WGCNA)
library(gridExtra)
library(ggplot2)
library(gplots)
library(baySeq)
library(edgeR)
library(DESeq)
library(NBPSeq)
library(ROC)
library(mclust)

#################################################################################################
print("Read in, update, and properly format sample information and gene expression data")

suppressWarnings(dir.create(outputFolder))
setwd(outputFolder)

print("-----Read in the gene expression data.")
datExprN   = read.csv(paste(inputFolder,"fpkm_table_normalized.csv",sep=""),row.names=1)   # RIN/batch-normalized FPKM
datExprU   = read.csv(paste(inputFolder,"fpkm_table_unnormalized.csv",sep=""),row.names=1) # Unnormalized FPKM
geneInfo   = read.csv(paste(inputFolder,"rows-genes.csv",sep=""),row.names=1)

print("-----Read in the meta-data and QC info about each sample.")
sampleInfo = read.csv(paste(inputFolder,"columns-samples.csv",sep=""),row.names=1)
colnames(datExprN) <- colnames(datExprU) <- rownames(sampleInfo)
extraInfo  = read.csv(paste(extraFolder,"tissueMetrics.csv",sep=""),row.names=1)  # Batch and RIN information
extraInfo  = extraInfo[rownames(sampleInfo),]
sampleInfo$RIN     = extraInfo$RIN
sampleInfo$Batch   = extraInfo$Batch
rm(extraInfo)

print("-----Re-order samples by brain region and then donor.")
ord        = order(sampleInfo$structure_acronym,sampleInfo$donor_name)
sampleInfo = sampleInfo[ord,]
datExprN   = datExprN[,ord]
datExprU   = datExprU[,ord]

print("-----Read in (and re-order) pathology info (i.e., Luminex and IHC) and donor information.")
metricInfo = read.csv(paste(inputFolder,"ProteinAndPathologyQuantifications.csv",sep=""))
metricInfo = metricInfo[order(metricInfo$structure_acronym,metricInfo$donor_name),]
donorInfo  = read.csv(paste(inputFolder,"DonorInformation.csv",sep=""),row.names=1)
donorInfo  = donorInfo[as.character(metricInfo$donor_id),]
rownames(donorInfo) <- rownames(metricInfo) <- samples <- rownames(sampleInfo)
donorInfo$donor_id = metricInfo$donor_id
donorInfo$structure_acronym = metricInfo$structure_acronym

print("-----Convert age ranges to numeric values corresponding to the median.")
ageNum     = as.character(donorInfo$age)
ageNum[ageNum=="90-94"] = "92"
ageNum[ageNum=="95-99"] = "97"
ageNum[ageNum=="100+"]  = "101"
donorInfo$ageNumeric    = as.numeric(as.character(ageNum))

print("-----Update TBI information for four donors who were initially thought not to have had a TBI.")
# Here is the specific changes that are made from the original file on the website based on new information
# tbiInfo = donorInfo[,c(13,6,9,7,1,20)]
# tbiInfo[donorInfo$name=="H14.09.072",1] = "Y"
# tbiInfo[donorInfo$name=="H14.09.072",2] = 50
# tbiInfo[donorInfo$name=="H14.09.072",3] = 1
# tbiInfo[donorInfo$name=="H14.09.072",4] = "10 sec - 1 min"
# tbiInfo[donorInfo$name=="H14.09.018",1] = "Y"
# tbiInfo[donorInfo$name=="H14.09.018",2] = 75
# tbiInfo[donorInfo$name=="H14.09.018",3] = 1
# tbiInfo[donorInfo$name=="H14.09.018",4] = "3-5 min"
# tbiInfo[donorInfo$name=="H14.09.038",1] = "Y"
# tbiInfo[donorInfo$name=="H14.09.038",2] = 64
# tbiInfo[donorInfo$name=="H14.09.038",3] = 1
# tbiInfo[donorInfo$name=="H14.09.038",4] = "10 sec - 1 min"
# tbiInfo[donorInfo$name=="H14.09.034",1] = "Y"
# tbiInfo[donorInfo$name=="H14.09.034",2] = 85
# tbiInfo[donorInfo$name=="H14.09.034",3] = 1
# tbiInfo[donorInfo$name=="H14.09.034",4] = "3-5 min"
# write.csv(tbiInfo,paste(extraFolder,"updatedTBIMetrics.csv",sep=""))  # Add name to first column and reorder

tbiInfo = read.csv(paste(extraFolder,"updatedTBIMetrics.csv",sep=""),row.names=1)
tbiInfo = tbiInfo[samples,]
donorInfo[,colnames(tbiInfo)] = tbiInfo
rm(tbiInfo)

print("-----Since one key piece of information is in dispute, it is removed from the analysis.")
print("Remove a particular sample (H14.09.011) that clusters as male, but is actually a female.")
kpSamp     = donorInfo$name!="H14.09.011"
sampleInfo = sampleInfo[kpSamp,]
datExprN   = datExprN[,kpSamp]
datExprU   = datExprU[,kpSamp]
donorInfo  = donorInfo[kpSamp,]
metricInfo = metricInfo[kpSamp,]

print("-----Re-level a few of the variables that are sorted alphabetically.")
l = c("< 10 sec","10 sec - 1 min","1-2 min","3-5 min","6-9 min","10 min - 1 hr","> 1 hr","Unknown or N/A")
donorInfo$longest_loc_duration = factor(as.character(donorInfo$longest_loc_duration),levels=l)
donorInfo$act_demented = factor(as.character(donorInfo$act_demented),levels=c("No Dementia","Dementia"))
l = c("No Dementia","Alzheimer's Disease Type","Vascular","Multiple Etiologies","Other or Unknown Cause","Other Medical")
donorInfo$dsm_iv_clinical_diagnosis = factor(as.character(donorInfo$dsm_iv_clinical_diagnosis),levels=l)
l = c("No dementia","Possible Alzheimer's Disease","Probable Alzheimer's Disease","Dementia, type unknown")
donorInfo$nincds_arda_diagnosis = factor(as.character(donorInfo$nincds_arda_diagnosis),levels=l)
donorInfo$apo_e4_allele = factor(as.character(donorInfo$apo_e4_allele),levels=c("N","Y","N/A"))

print("-----Save the names of some commonly used variables for later.")
region        = as.character(sampleInfo$structure_acronym)
regions       = c("HIP","TCx","PCx","FWM")
region        = factor(region,levels=regions)
batch         = as.character(sampleInfo$Batch)
sampleRIN     = as.numeric(as.character(sampleInfo$RIN))
TBIorControl  = as.character(donorInfo$ever_tbi_w_loc)
DementedYesNo = substr(as.character(donorInfo$act_demented),1,1)
DementedYesNo[DementedYesNo=="D"] = "Y"


print("-----Create separate tables for containing various combinations of data metrics for convenience.")
siCol = c("Batch", "hemisphere")
diCol = c("sex", "ageNumeric", "education_years",
  "act_demented", "nincds_arda_diagnosis", "dsm_iv_clinical_diagnosis",
  "braak", "cerad", "nia_reagan", "apo_e4_allele",
  "ever_tbi_w_loc", "num_tbi_w_loc", "longest_loc_duration", "age_at_first_tbi")
miCol=c("ihc_at8", "ihc_at8_ffpe", "ihc_tau2_ffpe", "ihc_a_beta", "ihc_a_beta_ffpe",
  "ihc_a_syn", "ihc_ptdp_43_ffpe", "ihc_iba1_ffpe", "ihc_gfap_ffpe",
  "ptau_ng_per_mg", "tau_ng_per_mg", "ab40_pg_per_mg", "ab42_pg_per_mg", "a_syn_pg_per_mg",
  "vegf_pg_per_mg", "tnf_a_pg_per_mg", "rantes_pg_per_mg", "bdnf_pg_per_mg", "mcp_1_pg_per_mg", 
  "il_10_pg_per_mg", "il_6_pg_per_mg", "il_7_pg_per_mg", "ifn_g_pg_per_mg", "mip_1a_pg_per_mg")

comparisonInfo = cbind(sampleInfo[,siCol],donorInfo[,diCol],metricInfo[,miCol])
comparisonInfo$ADorControl = comparisonInfo$dsm_iv_clinical_diagnosis
comparisonInfo$ADorControl[!is.element(comparisonInfo$dsm_iv_clinical_diagnosis,c("No Dementia","Alzheimer's Disease Type"))] = NA
comparisonInfo$ADorControl = droplevels(comparisonInfo$ADorControl)

abetaVals  = c("act_demented","dsm_iv_clinical_diagnosis","cerad","ihc_a_beta","ihc_a_beta_ffpe","ab42_pg_per_mg","ab40_pg_per_mg",
               "apo_e4_allele","braak","ihc_at8","ihc_at8_ffpe","ihc_tau2_ffpe","ptau_ng_per_mg","tau_ng_per_mg","ihc_a_syn","ageNumeric")
abetaInfo  = comparisonInfo[,c(abetaVals)]
abetaInfo[!is.element(abetaInfo[,"dsm_iv_clinical_diagnosis"],c("No Dementia","Alzheimer's Disease Type")),"dsm_iv_clinical_diagnosis"] = NA
abetaInfo[,"dsm_iv_clinical_diagnosis"] = droplevels(abetaInfo[,"dsm_iv_clinical_diagnosis"])
abetaInfo$cerad03 = abetaInfo$cerad
abetaInfo$cerad03[is.element(abetaInfo$cerad,c(1,2))] = NA
abetaInfo$cerad03 = factor(abetaInfo$cerad03)
abetaInfo$braakCombine = rep("3-4",length(abetaInfo$braak))
abetaInfo$braakCombine[abetaInfo$braak<3] = "0-2"
abetaInfo$braakCombine[abetaInfo$braak>4] = "5-6"
abetaInfo$braakCombine = factor(abetaInfo$braakCombine)


## Region colors
regCols   = c("#377eb8","#984ea3","#4daf4a","#e41a1c")
names(regCols) = regions
regionColors   = regCols[as.character(sampleInfo$structure_acronym)]


#################################################################################################
print("Only consider genes here at least 10% of samples from at least 1 brain region have FPKM > 2.")

findFromGroups <- function(datExpr,groupVector,fn="mean"){
  groups   = names(table(groupVector))
  fn       = match.fun(fn)
  datMeans = matrix(0,nrow=dim(datExpr)[2],ncol=length(groups))
  for (i in 1:length(groups)){
    datIn = datExpr[groupVector==groups[i],]
    if (is.null(dim(datIn)[1])) { datMeans[,i] = as.numeric(datIn)
    } else { datMeans[,i] = as.numeric(apply(datIn,2,fn)) }
  };    colnames(datMeans)  = groups;
  rownames(datMeans) = colnames(datExpr)
  return(datMeans)
}

zeroCount  = findFromGroups(t(datExprU>0),region,sum)
pass0Count = apply(zeroCount,1,min)>1
# This avoids some errors at calculating standard deviation later in the code and only fails 1 gene.

q90ByReg   = findFromGroups(t(datExprU),region,function(x) return(quantile(x,probs=0.9)))
isPresent  = (apply(q90ByReg,1,max) >= 2)&pass0Count

# Alternative strategy: Only consider genes with average FPKM > 1 in at least one region and at least 2 samples with FPKM > 0.
# meanExprU  = findFromGroups(t(datExprU),region)
# isPresent2 = (apply(meanExprU,1,max)>=1)&pass0Count




