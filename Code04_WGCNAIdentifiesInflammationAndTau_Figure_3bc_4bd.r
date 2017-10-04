#################################################################################################
##  These scripts take data from the the Aging, Dementia, and TBI website along with supplementary 
##  materials and uses it to perform all analyses for the TBI manuscript 
#################################################################################################

# R version 3.2.5 (2016-04-14) -- "Very, Very Secure Dishes" was used for the analysis, but the code should
# be compatible with most of the recent versions of R.

#################################################################################################
print("***** -----------------------------------------------------------------------------------------------------")
print("***** Code #4: WGCNA on RNA-Seq data identifies gene co-expression modules related to inflammation and tau.")
print("***** -----------------------------------------------------------------------------------------------------")

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################
## BEGIN FUNCTIONS

runRegionalWGCNA <- function(datExpr,kpSubset = 1:dim(datExpr)[2], reg = "All", minKMEtoStay = 0.4, threshM = 0.1, 
 minModSize = 20, minCoreKME = 0.75, N = round(dim(datExpr)[1]/2), power = 14, deepSplit = 2, maxModCount = 20, 
 maxBlockSize = 37500, matchModules = NA, omitLabels = NULL, ignoreMods = NULL, genes = rownames(datExpr), 
 networkType = "signed", ctMarkers=NA){

 print("Rank the genes based on variability.")
 rownames(datExpr) = genes
 N = round(min(c(dim(datExpr)[1],N)) )
 datExprRo  = datExpr[,kpSubset]
 meanR0     = rowMeans(datExprRo)
 varR0      = apply(datExprRo,1,var)
 topN       = names(-sort(-varR0))[1:N]
 datExprRun = t(datExprRo[topN,])

 
 print(paste("Run automated WGCNA using the",N,"most variable genes."))
 varR0      = apply(datExprRo,1,var)
 topN       = names(-sort(-varR0))[1:N]
 datExprRun = t(datExprRo[topN,])
 datExprRun = as.data.frame(datExprRun)
 ignoreMods = c(paste(reg,"M00",sep="_"),ignoreMods)

 # THIS IS THE MAIN FUNCTION CALL FOR THE WGCNA
 blockRun = blockwiseModules(datExprRun, checkMissingData = TRUE, maxBlockSize = maxBlockSize,
    power = power, networkType = networkType, deepSplit = deepSplit, minModuleSize = minModSize, 
    minCoreKMESize = minModSize/3, minKMEtoStay = minKMEtoStay, mergeCutHeight = threshM, 
    numericLabels = TRUE, verbose = 1)
 unmergedLabels = blockRun$colors
 names(unmergedLabels) <- colnames(datExprRun)


 print("Merge modules if more than a specified number are in the data set.")
 mergedLabels = unmergedLabels
 modCount = length(setdiff(unique(mergedLabels),c(ignoreMods,"0")))
 while(modCount>maxModCount){
  mergedLabels = as.numeric(mergedLabels)
  mergedMEs = moduleEigengenes(datExprRun, colors=mergedLabels, verbose=0,excludeGrey=TRUE)$eigengenes
  corMEdiss = 1-cor(mergedMEs);  diag(corMEdiss)=1
  mLab = rownames(which(corMEdiss==min(corMEdiss),arr.ind=TRUE))
  mLab = substr(mLab,3,nchar(mLab))
  mergedLabels[mergedLabels == mLab[1]] = mLab[2]
  modCount = length(setdiff(unique(mergedLabels),c(ignoreMods,"0")))
  print(paste("===Merged",mLab[1],"and",mLab[2],"- Distance =",signif(min(corMEdiss),3),"-",modCount,"modules remain."))
 } 


 print("Assign each gene to a final module based on kME and output some representative plots and data.")
 colorMergedDeleted = mergedLabels
 names(colorMergedDeleted) = topN
 colorMergedDeleted = paste("M",colorMergedDeleted,sep="");    
 a = nchar(colorMergedDeleted)==2
 colorMergedDeleted[a] = paste("M0",substr(colorMergedDeleted[a],2,2),sep="")
 colorMergedDeleted = paste(reg,colorMergedDeleted,sep="_")

 MEs    = moduleEigengenes(datExprRun,colorMergedDeleted,excludeGrey=TRUE,verbose=0,grey=paste(reg,"M00",sep="_"))$eigengenes
 colnames(MEs) = substr(colnames(MEs),3,nchar(colnames(MEs)))
 rownames(MEs) = colnames(datExprRo)
 kME    = cor(t(datExprRo),MEs,use="p")
 kME[is.na(kME)] = 0
 mods   = c(colnames(kME),paste(reg,"M00",sep="_"))
 maxKme = apply(kME,1,max)
 moduleColors = apply(kME,1,function(x,y) return(y[which.max(x)[1]]),mods)
 moduleColors[maxKme<minKMEtoStay] = paste(reg,"M00",sep="_")


 if(is.na(matchModules[1])){
  print("Relabel modules and order based on the (percent neuronal) - (percent glial) content of module.")
  bl  = as.character(unique(ctMarkers$CellType))
  cnt = table(ctMarkers[,1])[bl]
  comparisons = matrix(FALSE,nrow=length(colorMergedDeleted),ncol=length(bl))
  rownames(comparisons) = names(unmergedLabels)
  colnames(comparisons) = bl
  for (b in bl)
    comparisons[is.element(rownames(comparisons),rownames(ctMarkers)[ctMarkers[,1]==b]),b] = TRUE
  colorAssignedF = factor(colorMergedDeleted,levels = names(table(colorMergedDeleted)))
  names(colorAssignedF) = names(unmergedLabels)
  out = table(colorAssignedF)
  for (b in bl) out = cbind(out, table(colorAssignedF[comparisons[,b]]))
  colnames(out) = c("TotalGenes",bl)

  out2 = (-100/out[,1])*(out[,3]/cnt[2]-out[,2]/cnt[1]-out[,4]/cnt[3]-out[,5]/cnt[4])
  out2 = out2[setdiff(names(out2),ignoreMods)]
  ord      = order(out2)
  labels   = c(colnames(MEs),paste(reg,"M00",sep="_"))
  names(labels) = c(colnames(MEs)[ord],paste(reg,"M00",sep="_"))

  colorMergedDeleted = labels[colorMergedDeleted]
  seqs    = paste("M",c(paste("0",0:9,sep=""),as.character(10:99)),sep="")
  labels  = seqs[1:(length(unique(colorMergedDeleted)))]
  labels  = paste(reg,labels,sep="_")
  names(labels)      = sort(unique(colorMergedDeleted))
  moduleColors       = labels[moduleColors]
  colorMergedDeleted = labels[colorMergedDeleted]
  
 } else {
  print("Reassign module colors to roughly match module ordering from other specified module labels.")
  ovTable = overlapTable(moduleColors,matchModules[names(moduleColors)],ignore=c(ignoreMods,"M00"))
  countTable = ovTable$countTable;  
  pTable  = -log10(ovTable$pTable+10^(-299))+countTable/1000;
  pTable[grep("_M",rownames(pTable)),grep("_C",colnames(pTable))] = 0
  pTable[grep("_C",rownames(pTable)),grep("_M",colnames(pTable))] = 0
  ord     = order(apply(pTable,1,which.max))
  seqs    = paste("M",c(paste("0",0:9,sep=""),as.character(10:99)),sep="")
  labels  = seqs[1:(length(rownames(pTable))+1)]
  labels  = paste(reg,labels,sep="_")
  names(labels)      = c(paste(reg,"M00",sep="_"),rownames(pTable)[ord])
  colorMergedDeleted = labels[colorMergedDeleted]
  moduleColorsOld2   = moduleColors
  moduleColors       = labels[moduleColors]
  moduleColorsOld3   = moduleColors
 } 

 
 print("Finalize and save the modules.")
 MEs    = moduleEigengenes(datExprRun,colorMergedDeleted,excludeGrey=TRUE,verbose=0,grey=paste(reg,"M00",sep="_"))$eigengenes
 colnames(MEs) = substr(colnames(MEs),3,nchar(colnames(MEs)))
 rownames(MEs) = colnames(datExprRo)
 kME    = cor(t(datExprRo),MEs,use="p")
 kME[is.na(kME)] = 0
 mods   = colnames(kME)
 maxKme = apply(kME,1,max)
 moduleColors = apply(kME,1,function(x,y) return(y[which.max(x)[1]]),mods)
 moduleColors[maxKme<minKMEtoStay] = paste(reg,"M00",sep="_")
 kMEtable = data.frame(gene=rownames(kME),kME=maxKme,module=moduleColors)
 kMEtable = kMEtable[order(factor(moduleColors,levels=mods),-maxKme),]
 kME      = kME[rownames(kMEtable),mods]

 names(mergedLabels) <- names(colorMergedDeleted) <- names(unmergedLabels)
 write.csv(kMEtable,paste("WGCNA_",reg,"_MEcorrelationTable_kME.csv",sep=""),row.names=FALSE)
 save(MEs,mergedLabels,colorMergedDeleted,blockRun,moduleColors,file=paste("WGCNA_",reg,"_Robject.RData",sep=""))

 
 print("Output dendrograms of the data.")
 pdf(paste("WGCNA_",reg,"_dendrogram.pdf",sep=""),height=8,width=15)
 mName  = c("Intial Modules","Final Assignments")
 cs     = c("grey",standardColors())
 mColor = cbind(labels2colors(mergedLabels),labels2colors(moduleColors[topN],colorSeq=cs))
 plotDendroAndColors(blockRun$dendrograms[[1]],mColor,mName,dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
 dev.off()

 return(list(MEs=MEs,mergedLabels=mergedLabels,colorMergedDeleted=colorMergedDeleted,blockRun=blockRun,
        moduleColors=moduleColors,kMEtable=kMEtable))
}


## END FUNCTIONS
######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

print("=====Perform WGCNA with the goal of finding marker genes and co-expressed groups of genes.")

print("Read in the cell type markers for module ordering and cell type comparison.")
barresLabMarkers = read.csv(paste(extraFolder,"barreslab_rnaseq_DEX_human.csv",sep="/"),row.names=1)
neMarkers        = read.csv(paste(extraFolder,"neuroExpresso_human.csv",sep="/"))

print("Use temporal cortex as baseline region for matching module order.")
wgcna      = list()
r          = "TCx"
datExpr    = datExprRg
genes      = rownames(datExpr)
wgcnaFolder= paste(outputFolder,"WGCNA_",r,"/",sep="")
dir.create(wgcnaFolder)
setwd(wgcnaFolder)
wgcna[[r]] = runRegionalWGCNA(datExpr,kpSubset = region==r, reg = r, genes = genes, ctMarkers=barresLabMarkers) # DEFAULT maxModCount = 20.

print("Run WGCNA on other three regions")
matchModules = wgcna[[r]][["kMEtable"]]$module
names(matchModules) = wgcna[[r]][["kMEtable"]]$gene
matchModules = matchModules[genes]
for (r in c("HIP","PCx","FWM")){
   wgcnaFolder= paste(outputFolder,"WGCNA_",r,"/",sep="")
   dir.create(wgcnaFolder)
   setwd(wgcnaFolder)
   wgcna[[r]] = runRegionalWGCNA(datExpr,kpSubset = region==r, reg = r, genes = genes, matchModules = matchModules, ctMarkers=barresLabMarkers)
}
setwd(outputFolder)

print("For each pair of regions, determine the module overlap.")
modAssigns = NULL
for (r in regions)
 modAssigns = cbind(modAssigns,wgcna[[r]][["moduleColors"]])
colnames(modAssigns) = regions
rownames(modAssigns) = genes
ignoreMods = paste(regions,"M00",sep="_")

write.csv(modAssigns,"SupplementaryTable_3_WGCNA_module_assignments.csv")

pdf("WGCNA_moduleComparisonBetweenRegions2.pdf", height=9,width=15)  
for (r1 in regions) for (r2 in regions){
 ovTable = overlapTable(modAssigns[,r1],modAssigns[,r2],ignore=ignoreMods)
 pTable  = -log10(ovTable$pTable+10^(-299));    
 pTable2 = pmin(pTable,100);     countTable = ovTable$countTable;   
 xLabels = colnames(pTable);     yLabels    = rownames(pTable)

 par(mfrow=c(1,1));   par(cex = 1.0);   par(mar=c(8, 12.4, 2.7, 1)+0.3);
 labeledHeatmap(Matrix = pTable2, xLabels = paste(" ",xLabels), yLabels = paste(" ",yLabels), 
  colorLabels = TRUE, ySymbols = yLabels, main = " ", xSymbols = xLabels, textMatrix = countTable, 
  colors = greenWhiteRed(100)[50:100], setStdMargins = FALSE, cex.text = 1.0, cex.lab = 1.0);
}
dev.off()


#################################################################################################
print("Identify marker genes for cell type using Zhang ... Barres data for comparison.")

phyper2 <- function (total, group1, group2, overlap, verySig=TRUE ,lt=TRUE){
  # This function is the same is phyper, just allows for more sensible input values
  q = overlap
  if(q==0) return(1)
  m = group1
  n = total-group1
  k = group2
  prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
  if (verySig) return(-prob)
  return(1-prob)
}

cellTypeNames    = as.character(unique(barresLabMarkers[,1]))
cellTypeMarkers <- mods <- NULL
for (r in regions) for (m in sort(unique(modAssigns[,r]))){
 cmp = NULL
 for (b in cellTypeNames){
  inM = rownames(modAssigns)[modAssigns[,r]==m]
  inB = rownames(barresLabMarkers)[barresLabMarkers[,1]==b]
  cmp = c(cmp,phyper2(dim(datExprRg)[1],length(inM),length(inB),length(intersect(inM,inB))))
 }
 mods = c(mods,m)
 cellTypeMarkers = rbind(cellTypeMarkers,cmp)
}
rownames(cellTypeMarkers) = mods
colnames(cellTypeMarkers) = cellTypeNames

cellTypeMarkers = pmin(cellTypeMarkers*prod(dim(cellTypeMarkers)),1)

write.csv(cellTypeMarkers,"SupplementaryTable_4_moduleEnrichmentForCellTypes.csv")

#################################################################################################
print("Identify marker genes for cell type using neuroExpresso data for comparison.")

neMarkers2         = neMarkers[is.element(neMarkers[,1],rownames(modAssigns)),]
cellTypeNENames    = as.character(unique(neMarkers2[,2]))
cellTypeNEMarkers <- mods <- NULL
for (r in regions) for (m in sort(unique(modAssigns[,r]))){
 cmp = NULL
 reg = ifelse(r=="HIP","Hippocampus","Cortex")
 for (b in cellTypeNENames){
  inM = rownames(modAssigns)[modAssigns[,r]==m]
  inB = as.character(neMarkers2[,1][(neMarkers2[,2]==b)&(neMarkers2[,3]==reg)])
  cmp = c(cmp,phyper2(dim(datExprRg)[1],length(inM),length(inB),length(intersect(inM,inB))))
 }
 mods = c(mods,m)
 cellTypeNEMarkers = rbind(cellTypeNEMarkers,cmp)
}
rownames(cellTypeNEMarkers) = mods
colnames(cellTypeNEMarkers) = cellTypeNENames

cellTypeNEMarkers = pmin(cellTypeNEMarkers*437,1) # prod(dim(cellTypeNEMarkers)),1)  # Some clusters are too small or only in hippocampus

write.csv(cellTypeNEMarkers,"SupplementaryTable_4b_moduleEnrichmentForCellTypeNEs_NE.csv")



#################################################################################################
print("Find all MEs associated with the metrics using SVA, and show correlations for signficant results.")
print("-----Identify metrics for comparison with module eigengenes using SVA")
print("-----Include select pathology metrics for Abeta and Tau, inflammatory markers, and demographic info.")
cn      = c("Age","act_demented2","dsm_iv_hasAD","casi_irt",
"braak","braakCombine","ihc_at8","ihc_at8_ffpe","ihc_tau2_ffpe","ptau_ng_per_mg","tau_ng_per_mg",
"cerad","ihc_a_beta","ihc_a_beta_ffpe","ab42_pg_per_mg")
an = c("tnf_a_pg_per_mg","il_6_pg_per_mg","mip_1a_pg_per_mg","mcp_1_pg_per_mg","rantes_pg_per_mg",
"il_7_pg_per_mg","ihc_iba1_ffpe","ihc_gfap_ffpe")
dn = c("sex","apo_e4_allele","education_years")
svaInfo2 = cbind(abetaInfo[,cn],metricInfo[,an],donorInfo[,dn])

print("-----Cap outlier values at 3 SD above the mean (after excluding outlier samples).")
for (i in 1:dim(svaInfo2)[2]) if(class(svaInfo2[,i])=="numeric") for (r in regions) {
  kp = region==r
  s3 = mean(svaInfo2[kp,i],na.rm=TRUE)+3*sd(svaInfo2[kp,i],na.rm=TRUE)
  tmp = svaInfo2[,i];
  tmp[is.na(tmp)] = -10000
  kp2 = (region==r)&(tmp<s3)
  s33 = mean(svaInfo2[kp2,i],na.rm=TRUE)+3*sd(svaInfo2[kp2,i],na.rm=TRUE)
  kp3 = (region==r)&(tmp>=s33)
  svaInfo2[kp3,i] = s33
}

print("-----Run the SVA analysis and plot significant modules.")
svaMEs   = list() 
coln     = colnames(svaInfo2)
numTests = length(coln)*
 (length(colnames(wgcna[["HIP"]][["MEs"]]))+length(colnames(wgcna[["TCx"]][["MEs"]]))+
 length(colnames(wgcna[["PCx"]][["MEs"]]))+length(colnames(wgcna[["FWM"]][["MEs"]])))
for (r in regions){
 svaTmp = matrix(0,nrow=length(colnames(wgcna[[r]][["MEs"]])),ncol=length(coln))
 rownames(svaTmp) = colnames(wgcna[[r]][["MEs"]]) 
 colnames(svaTmp) = coln
 for (cn in coln){
  kp  = (region==r)
  kp2 = (!is.na(svaInfo2[kp,cn]))
  siTmp = cbind(rep(1,length(svaInfo2[,cn])),svaInfo2[,cn])
  if(class(siTmp[,1])=="factor")  siTmp = droplevels(siTmp)
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprN)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,][kp2,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,][kp2,])
  pValuesBatch    = f.pvalue(as.matrix(t(wgcna[[r]][["MEs"]][kp2,])),modBatch,mod0Batch)
  qValuesBatch    = pValuesBatch*numTests
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,cn]     = qValuesBatch
 }
 svaMEs[[r]] = svaTmp
}

## Output a table of the p-values
svaOut = NULL
for (r in regions) svaOut = rbind(svaOut,svaMEs[[r]])
write.csv(svaOut,"SupplementaryTable_5_WGCNA_moduleVsMetricCorrelations.csv")

## Make the plot
pThresh = 0.05
pdf("Figure_3b_4b_CorrelationOfMetricsWithSignficantMEs.pdf")
for (r in regions) for (cn in coln) for (mn in colnames(wgcna[[r]][["MEs"]])) if(svaMEs[[r]][mn,cn]<pThresh) {
 kp   = (region==r)
 px   = as.numeric(svaInfo2[kp,cn])
 gy   = wgcna[[r]][["MEs"]][,mn]
 xlab = paste(cn,"quant. - p =",signif(svaMEs[[r]][mn,cn],3))
 verboseScatterplot(px,gy,col=regCols[r],main=mn,ylab="ME expression",xlab=xlab,pch=19)
}
dev.off()


#################################################################################################
print("-----Run the SVA analysis on all genes to determine how many genes are related to each metric.")

datTmp   = list()
for (r in regions) datTmp[[r]] = t(datExprRg[,region==r])
svaGNs   = list() 
coln     = colnames(svaInfo2)
numTests = length(coln)*
 (length(colnames(datTmp[["HIP"]]))+length(colnames(datTmp[["TCx"]]))+
 length(colnames(datTmp[["PCx"]]))+length(colnames(datTmp[["FWM"]])))
for (r in regions){
 svaTmp = matrix(0,nrow=length(colnames(datTmp[[r]])),ncol=length(coln))
 rownames(svaTmp) = colnames(datTmp[[r]]) 
 colnames(svaTmp) = coln
 for (cn in coln){
  kp  = (region==r)
  kp2 = (!is.na(svaInfo2[kp,cn]))
  siTmp = cbind(rep(1,length(svaInfo2[,cn])),svaInfo2[,cn])
  if(class(siTmp[,1])=="factor")  siTmp = droplevels(siTmp)
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprN)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,][kp2,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,][kp2,])
  pValuesBatch    = f.pvalue(as.matrix(t(datTmp[[r]][kp2,])),modBatch,mod0Batch)
  qValuesBatch    = pValuesBatch*numTests
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,cn]     = qValuesBatch
 }
 svaGNs[[r]] = svaTmp
}

DEXcount = NULL
for (r in regions)
 DEXcount = cbind(DEXcount,colSums(svaGNs[[r]]<0.05))
colnames(DEXcount) = regions
print(DEXcount[rowSums(DEXcount)>0,])
# Very few genes are significantly correlated with anything, likely due to the strict multiple correction criteria

# OLD ANALYSIS (DISREGARD)
#                  HIP TCx PCx FWM
# ptau_ng_per_mg     0   0   0   1
# il_6_pg_per_mg     0  50  31  42
# mip_1a_pg_per_mg  30   0   0   0
# mcp_1_pg_per_mg    0  18   1   2
# sex               39  38  39  34

# CURRENT ANALYSIS
#                  HIP TCx PCx FWM
# ihc_at8            0   0   1   0
# ihc_at8_ffpe       0   0   0   1
# ptau_ng_per_mg     0   0   0   1
# il_6_pg_per_mg     0  21  13   8
# mip_1a_pg_per_mg  14   7  12   0
# mcp_1_pg_per_mg    4  13   2   1
# sex               39  38  39  34

#################################################################################################
print("Determine whether inflammation is regional/global based on gene and protein markers.")

donorsU = sort(unique(donorInfo$name))
inflamG  = matrix(NA,nrow=length(donorsU),ncol=4)
rownames(inflamG) = donorsU
infModules  = c("HIP_M18","TCx_M18","PCx_M18","FWM_M13")
infProteins = c("mip_1a_pg_per_mg",rep("il_6_pg_per_mg",3))
colnames(inflamG) <- names(infModules) <- names(infProteins) <-regions
inflamP = inflamG

for (r in regions) {
 kp      = (region==r)
 dn      = as.character(donorInfo$name[kp])
 inflamG[dn,r] = wgcna[[r]]$MEs[,infModules[r]]
 inflamP[dn,r] = svaInfo2[kp,infProteins[r]]
}

# Output p-values for the poorest correlations of the pairs for inclusion in the text
print(cor.test(inflamP[,1],inflamP[,4],use="P"))
print(cor.test(inflamG[,1],inflamG[,4],use="P"))

plotCorText = function(x,y,...){
 cr = as.numeric(signif(cor(x,y,use="p"),2))
 text(mean(range(x,na.rm=TRUE)),mean(range(y,na.rm=TRUE)),cr,cex=2)
}

pdf("Figure_3c_RegionalCorrelationOfInflammatoryGenesAndProteins.pdf")
pairs(inflamG,main="Genes",pch=19,cex=0.75,upper.panel=plotCorText,lower.panel=points)
pairs(inflamP,main="Proteins",pch=19,cex=0.75,upper.panel=points,lower.panel=plotCorText)
dev.off()

print("Output the values marking gene and protein inflammation across regions.")
inflamPG = cbind(inflamP,inflamG)
colnames(inflamPG) = paste(colnames(inflamPG),c(rep("(protein)",4),rep("(gene)",4)))
write.csv(signif(inflamPG,4),"zNotAFigure_inflamPG.csv")


#################################################################################################
print("Compare gene expresson signatures of tau in this data set with gene expression signatures in Blalock et al 2008.")

datTmp  = read.csv(paste(extraFolder,"Blalock2004_geneExpression.csv",sep=""),row.names=1)
sampTmp = read.csv(paste(extraFolder,"Blalock2004_sampleInfo.csv",sep=""),row.names=1) 
if(substr(colnames(datTmp),1,1)[1]=="X")  rownames(sampTmp) = paste("X",rownames(sampTmp),sep="")
kpS     = intersect(rownames(sampTmp),colnames(datTmp))
kpR     = region=="HIP"
kpG     = intersect(rownames(datExprRg),rownames(datTmp))
corCP   = apply(2^datExprRg[kpG,kpR]-1,1,cor,as.numeric(as.character(metricInfo[kpR,"ihc_at8"])))
corO    = apply(as.matrix(datTmp[kpG,kpS]),1,cor,as.numeric(as.character(sampTmp[kpS,"NFT"])),use="p")
plotCol = rep("grey",length(corO))
names(plotCol) = names(corO)
kme = wgcna[["HIP"]]$kMEtable
plotCol[is.element(names(corO),kme[kme[,"module"]=="HIP_M16","gene"])] = "black"

pdf("Figure_4d_comparisonBetweenCurrentDataAndBlalock2004.pdf",width=8,height=8) 
verboseScatterplot(corCP,corO,xlab="IHC AT8 in hippocampus",ylab="NFT quant in Blalock 2004",
  main="",pch=19,cex=0.5+0.5*(plotCol=="black"),col=plotCol)
dev.off()


