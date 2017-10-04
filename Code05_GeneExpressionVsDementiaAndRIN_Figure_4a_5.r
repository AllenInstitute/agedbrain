#################################################################################################
##  These scripts take data from the the Aging, Dementia, and TBI website along with supplementary 
##  materials and uses it to perform all analyses for the TBI manuscript 
#################################################################################################

# R version 3.2.5 (2016-04-14) -- "Very, Very Secure Dishes" was used for the analysis, but the code should
# be compatible with most of the recent versions of R.

#################################################################################################
print("***** ---------------------------------------------------------------------------------------------------")
print("***** Code #5: Produce figures and plots related to gene expression vs. dementia and RIN for Figures 4-5.")
print("***** ---------------------------------------------------------------------------------------------------")


###################################################################################################
print("Compare RIN to dementia status and AD vs. control.")
print("-----(Note that the jitter in these plots is arbitrary and may not match the manuscript.)")

cnsTmp = c("act_demented","ADorControl")
sn     = "RIN"
pdf("Figure_5a_RIN_comparison_plots_Dementia.pdf",height=4,width=8)
grobb = list()
for (cn in cnsTmp) for (r in regions) {
 comr = cbind(comparisonInfo[region==r,cnsTmp],sampleInfo[region==r,sn])
 comr[,1] = gsub("No Dementia","1No",comr[,1])
 comr[,2] = gsub("No Dementia","1No",comr[,2])
 comr[,1] = gsub("Dementia","Dem.",comr[,1])
 comr[,2] = gsub("Alzheimer's Disease Type","AD",comr[,2])
 colnames(comr) = c(cnsTmp,sn)
 kp  = (!is.na(comr[,cn]))&(!is.na(comr[,sn]))
 pval= getAnovaPvalforApply(comr[kp,sn],comr[kp,cn])
 ddf = data.frame(NUMS = comr[kp,sn], GRP = as.character(comr[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("p =",signif(pval,3))) + labs(x=cn,y=paste(r,sn)) + ylim(2.2,9) + 
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=4,nrow=1))
dev.off()

pdf("zNotAFigure_RIN_histograms_Dementia.pdf",width=11,height=5)
par(mfrow=c(1,2))
hist(sampleRIN[DementedYesNo=="N"],xlim=c(2,9),ylim=c(0,50),main="",xlab="RIN, no dementia")
abline(v=1:10,lty="dotted")
hist(sampleRIN[DementedYesNo=="Y"],xlim=c(2,9),ylim=c(0,50),main="",xlab="RIN, dementia")
abline(v=1:10,lty="dotted")
dev.off()

#################################################################################################
print("Identify genes significantly associated with dementia, AD, or RIN using SVA.")

coln    = c("act_demented","ADorControl","RIN")
svaROI  = list() 
corMet  = "BH"
comparisonInfo$RIN = sampleInfo$RIN
for (cn in coln) {
 svaTmp = datExprRg[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprRg)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(as.matrix(datExprRg[,kp]),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaROI[[cn]] = svaTmp
}

svaROIu = list() 
for (cn in coln) {
 svaTmp = datExprTb2[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprTb2)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(as.matrix(datExprTb2[,kp]),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet)
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaROIu[[cn]] = svaTmp
}


###################################################################################################
print("Find and plot correlation between RIN and gene expression (of TbT normalized, but NOT RIN-corrected data).")

datExpr = datExprTb2
rinCor  = apply(datExpr,1,cor,sampleInfo$RIN)
rinCors = NULL
for (r in regions)
  rinCors = cbind(rinCors,apply(datExpr[,region==r],1,cor,sampleInfo$RIN[region==r]))
colnames(rinCors) = regions

pThresh = 0.05
cThresh = 0.5
isP     = (svaROIu[["RIN"]]<pThresh)&(rinCors>cThresh);  cnCntsP = colSums(isP)
isN     = (svaROIu[["RIN"]]<pThresh)&(rinCors< -cThresh);  cnCntsN = colSums(isN)

pdf("Figure_5b_histogramsOfRINcorrelations.pdf",height=5,width=6)
for (r in regions){
  hist(rinCors[,r],col="grey",xlab="RIN correlation",main=r,breaks=50,xlim=c(-1,1),ylim=c(0,900))
  abline(v=c(cThresh,-cThresh),col="red",lty="dashed"); text(-0.8,900,cnCntsN[r]); text(0.8,900,cnCntsP[r])
}
dev.off()


###################################################################################################
print("Find and plot log2 correlations between control/AD, control/dementia both with and without normalization.")

fcADn <- fcDemn <- fcADu <- fcDemu <- NULL
for (r in regions){
 tmp    = findFromGroups(t(datExprTb2[,region==r]),comparisonInfo$act_demented[region==r],function(x) return(mean(x,na.rm=TRUE)))
 fcADu  = cbind(fcADu,tmp[,2]-tmp[,1])
 tmp    = findFromGroups(t(datExprTb2[,region==r]),comparisonInfo$ADorControl[region==r],function(x) return(mean(x,na.rm=TRUE)))
 fcDemu = cbind(fcDemu,tmp[,2]-tmp[,1])
 tmp    = findFromGroups(t(datExprRg[,region==r]),comparisonInfo$act_demented[region==r],function(x) return(mean(x,na.rm=TRUE)))
 fcADn  = cbind(fcADn,tmp[,2]-tmp[,1])
 tmp    = findFromGroups(t(datExprRg[,region==r]),comparisonInfo$ADorControl[region==r],function(x) return(mean(x,na.rm=TRUE)))
 fcDemn = cbind(fcDemn,tmp[,2]-tmp[,1])
}
colnames(fcADn) <- colnames(fcADu) <- colnames(fcDemn) <- colnames(fcDemu) <- regions 

fThresh = log2(1.3)
isPdn   = (svaROI[["act_demented"]]<pThresh)&(fcDemn>fThresh);    cnCntsPdn = colSums(isPdn)
isNdn   = (svaROI[["act_demented"]]<pThresh)&(fcDemn< -fThresh);  cnCntsNdn = colSums(isNdn)
isPan   = (svaROI[["ADorControl"]]<pThresh)&(fcADn>fThresh);      cnCntsPan = colSums(isPan)
isNan   = (svaROI[["ADorControl"]]<pThresh)&(fcADn< -fThresh);    cnCntsNan = colSums(isNan)
isPdu   = (svaROIu[["act_demented"]]<pThresh)&(fcDemu>fThresh);   cnCntsPdu = colSums(isPdu)
isNdu   = (svaROIu[["act_demented"]]<pThresh)&(fcDemu< -fThresh); cnCntsNdu = colSums(isNdu)
isPau   = (svaROIu[["ADorControl"]]<pThresh)&(fcADu>fThresh);     cnCntsPau = colSums(isPau)
isNau   = (svaROIu[["ADorControl"]]<pThresh)&(fcADu< -fThresh);   cnCntsNau = colSums(isNau)

thVal = 0.45
ylim  = c(0,2500)
xlim  = c(-thVal,thVal)
pdf("Figure_4a_5c_histogramsOfDementiaFCs.pdf",height=5,width=6)
for (r in regions){
  hist(pmin(thVal,pmax(fcDemu[,r],-thVal)),col="grey",xlab="log2FC(Dementia), uncorrected",main=r,breaks=50,xlim=xlim,ylim=ylim)
  abline(v=c(fThresh,-fThresh),col="red",lty="dashed"); text(-fThresh*0.8,ylim[2],cnCntsNdu[r]); text(fThresh*0.8,ylim[2],cnCntsPdu[r])
  hist(pmin(thVal,pmax(fcADu[,r],-thVal)),col="grey",xlab="log2FC(AD), uncorrected",main=r,breaks=50,xlim=xlim,ylim=ylim)
  abline(v=c(fThresh,-fThresh),col="red",lty="dashed"); text(-fThresh*0.8,ylim[2],cnCntsNau[r]); text(fThresh*0.8,ylim[2],cnCntsPau[r])
  hist(pmin(thVal,pmax(fcDemn[,r],-thVal)),col="grey",xlab="log2FC(Dementia), RIN corrected",main=r,breaks=50,xlim=xlim,ylim=ylim)
  abline(v=c(fThresh,-fThresh),col="red",lty="dashed"); text(-fThresh*0.8,ylim[2],cnCntsNdn[r]); text(fThresh*0.8,ylim[2],cnCntsPdn[r])
  hist(pmin(thVal,pmax(fcADn[,r],-thVal)),col="grey",xlab="log2FC(AD), RIN corrected",main=r,breaks=50,xlim=xlim,ylim=ylim)
  abline(v=c(fThresh,-fThresh),col="red",lty="dashed"); text(-fThresh*0.8,ylim[2],cnCntsNan[r]); text(fThresh*0.8,ylim[2],cnCntsPan[r])
}
dev.off()


###################################################################################################
print("Output tables of all of these correlations and statistics and whether or not they are significant in each region.")

outU  = cbind(rinCors,fcDemu,fcADu,svaROIu[["RIN"]],svaROIu[["act_demented"]],svaROIu[["ADorControl"]],isP-isN,isPdu-isNdu,isPau-isNau)
colnames(outU) = c(paste("cor_RIN:",regions),paste("log2FC_Dementia:",regions),paste("log2FC_AD:",regions),
 paste("svaPval_RIN:",regions),paste("svaPval_Dementia:",regions),paste("svaPval_AD:",regions),
 paste("sig?_RIN:",regions),paste("sig?_Dementia:",regions),paste("sig?_AD:",regions))
 
outN  = cbind(fcDemn,fcADn,svaROI[["act_demented"]],svaROI[["ADorControl"]],isPdn-isNdn,isPan-isNan)
colnames(outN) = c(paste("log2FC_Dementia:",regions),paste("log2FC_AD:",regions),
 paste("svaPval_Dementia:",regions),paste("svaPval_AD:",regions),
 paste("sig?_Dementia:",regions),paste("sig?_AD:",regions))
  
write.csv(outU,"SupplementaryTable_6and7_GeneExpressionVsRINandDementia_uncorrected.csv")
write.csv(outN,"SupplementaryTable_7_GeneExpressionVsDementia_RINcorrected.csv")

  
###################################################################################################
print("Plot gene ranks between RIN correlation and AD vs. control fold change after accounting for RIN.")

pdf("zNotAFigure_scatterplotRINvsDementiaRankedExpression.pdf",height=7,width=7)
for (r in regions){
 xd = rank(fcDemn[,r])
 xa = rank(fcADn[,r])
 y  = rank(rinCors[,r])
 verboseScatterplot(xd,y,xlab="Rank Dementia vs. control",ylab="Rank RIN corr.",main=r,pch=19,cex=0.3)
 df <- data.frame(x = xd, y = y,
    d = densCols(xd, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
 p <- ggplot(df) +
    geom_point(aes(xd, y, col = d), size = 1) +
    scale_color_identity() +
    theme_bw()
 print(p)
 verboseScatterplot(xa,y,xlab="Rank AD vs. control",ylab="Rank RIN corr.",main=r,pch=19,cex=0.3)
 df <- data.frame(x = xa, y = y,
    d = densCols(xa, y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
 p <- ggplot(df) +
    geom_point(aes(xa, y, col = d), size = 1) +
    scale_color_identity() +
    theme_bw()
 print(p)
}
dev.off()


###################################################################################################
print("Plot gene ranks between AD vs. control fold change before vs. after accounting for RIN.")

pdf("Figure_5d_scatterplotDementiaUncorrectedVsCorrectedRanks.pdf",height=7,width=7)
for (r in regions){
 xd = rank(fcDemn[,r])
 xa = rank(fcADn[,r])
 yd = rank(fcDemu[,r])
 ya = rank(fcADu[,r])
 verboseScatterplot(xd,yd,xlab="Rank Dementia vs. control (RIN corrected)",ylab="Rank Dementia vs. control (uncorrected)",main=r,pch=19,cex=0.3)
 df <- data.frame(x = xd, y = yd,
    d = densCols(xd, yd, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
 p <- ggplot(df) +
    geom_point(aes(xd, yd, col = d), size = 1) +
    scale_color_identity() +
    theme_bw()
 print(p)
 verboseScatterplot(xa,ya,xlab="Rank AD vs. control (RIN corrected)",ylab="Rank AD vs. control (uncorrected)",main=r,pch=19,cex=0.3)
 df <- data.frame(x = xa, y = ya,
    d = densCols(xa, ya, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
 p <- ggplot(df) +
    geom_point(aes(xa, ya, col = d), size = 1) +
    scale_color_identity() +
    theme_bw()
 print(p)
}
dev.off()


###################################################################################################
print("Compare genes differentially expressed in AD from previous studies with our RIN correlations.")
print("-----For studies in hippocampus, compare with HIP; for all other studies, compare with TCx.")

ADGenes  = read.csv(paste(extraFolder,"AlzheimersGeneLists.csv",sep=""))
ADInfo   = read.csv(paste(extraFolder,"AlzhiemersDataSets.csv",sep=""),row.names=1)

plotInfo = list()
for (r in regions){
 plotInfo[[r]] = list()
 rinRank2 = rank(-rinCors[,r]);   plotInfo[[r]][["RIN"]] = rinRank2
 rankTmp  = rank(fcDemn[,r]);     plotInfo[[r]][["Dementia_Cor"]] = rankTmp[names(rinRank2)]
 rankTmp  = rank(fcDemu[,r]);     plotInfo[[r]][["Dementia_Unc"]] = rankTmp[names(rinRank2)]
 rankTmp  = rank(fcADn[,r]);      plotInfo[[r]][["AD_Cor"]] = rankTmp[names(rinRank2)]
 rankTmp  = rank(fcADu[,r]);      plotInfo[[r]][["AD_Unc"]] = rankTmp[names(rinRank2)]
}

print("-----Get the relevant statistics.")
## Function for collecting the statistics for comparing one gene list against another ranked list
getStatsForApply2 <- function(x,ranks){   # For ranks, lower numbers are BETTER
    gnIn = intersect(names(ranks),unique(x))
    rkIn = as.numeric(ranks[gnIn])
    per5 = round(length(ranks)/20)
    len  = length(ranks)
    num  = length(gnIn)
    quan = quantile(rkIn,probs=c(0.10,0.25,0.50,0.75,0.90))/len
    top5 = sum(rkIn<=per5)
    if(length(rkIn)>=5){
      wcox = wilcox.test(rkIn,setdiff(1:length(ranks),rkIn),alternative="two.sided",exact=FALSE,correct=TRUE)$p.value
    } else { wcox=1 }
    out  = c(num,quan,top5,wcox)
    names(out) = c("#Genes","10%","25%","Median","75%","90%","#top5%","pvalWilcox")
    return(out) 
}
## END OF FUNCTION

nm       = names(plotInfo[[1]])
outStats = list()
for (l in nm){
 outStats[[l]] = NULL
 for (n in rownames(ADInfo)){
  r            = as.character(ADInfo[n,"ComparisonRegion"])
  testGenes    = as.character(ADGenes[ADGenes[,2]==n,1])
  testCatStats = getStatsForApply2(testGenes,plotInfo[[r]][[l]])
  outStats[[l]] = rbind(outStats[[l]],testCatStats)
 }
 rownames(outStats[[l]]) = rownames(ADInfo)
 outStats[[l]][,2:6] = 1-outStats[[l]][,2:6]
}

print("-----Generate the plots.")

pdf("bbFigure_5e_ComparisonWithPreviousStudiesOfAD.pdf",height=7,width=8)
cols = c("red","green","blue","green","blue")
names(cols) = names(plotInfo[[1]])
x = ADInfo$xLocation-1
cexN = ifelse(ADInfo$Direction=="Down",19,5)
ltyN = ifelse(ADInfo$Direction=="Down","solid","dashed")
txtN = paste(ADInfo[,"Study"],"-",ADInfo[,"BrainRegion"])
txtN[ADInfo$Direction=="Down"] = ""
dff  = 0.2

nm = c("RIN","Dementia_Unc","Dementia_Cor")
plot(0,0,col="white",xlab="Normalized gene rank",ylab="",main="RIN and Dementia",xlim=c(0,max(x)),ylim=c(-1,1.1))
abline(h=c(1,0)); abline(h=0.5,lty="dotted")
for (j in 1:length(nm)){  
 l = nm[j]
 xo = x+(j-2)*dff
 points(xo,outStats[[l]][,"Median"],col=cols[l],pch=cexN)
 for(i in 1:dim(ADInfo)[1]){
  segments(xo[i],y0=outStats[[l]][i,"25%"],y1=outStats[[l]][i,"75%"],lty=ltyN[i],col=cols[l])
  segments(y0=outStats[[l]][i,"25%"],x0=xo[i]-0.05,x1=xo[i]+0.05,col=cols[l])
  segments(y0=outStats[[l]][i,"75%"],x0=xo[i]-0.05,x1=xo[i]+0.05,col=cols[l])
 }
}
text(x,rep(-0.5,length(x)),txtN,srt=90,cex=0.75)

nm = c("RIN","AD_Unc","AD_Cor")
plot(0,0,col="white",xlab="Normalized gene rank",ylab="",main="RIN and AD",xlim=c(0,max(x)),ylim=c(-1,1.1))
abline(h=c(1,0)); abline(h=0.5,lty="dotted")
for (j in 1:length(nm)){  
 l = nm[j]
 xo = x+(j-2)*dff
 points(xo,outStats[[l]][,"Median"],col=cols[l],pch=cexN)
 for(i in 1:dim(ADInfo)[1]){
  segments(xo[i],y0=outStats[[l]][i,"25%"],y1=outStats[[l]][i,"75%"],lty=ltyN[i],col=cols[l])
  segments(y0=outStats[[l]][i,"25%"],x0=xo[i]-0.05,x1=xo[i]+0.05,col=cols[l])
  segments(y0=outStats[[l]][i,"75%"],x0=xo[i]-0.05,x1=xo[i]+0.05,col=cols[l])
 }
}
text(x,rep(-0.5,length(x)),txtN,srt=90,cex=0.75)
dev.off()


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
print("Run several additional analyses comparing control and dementia using different statistical tests and data sets.")

##### ---------------------------------------------------------------------------------------------
print("===== SVA in normalized data, Bonferroni corrected")
coln    = c("act_demented","ADorControl")
svaOut  = list() 
corMet  = "bonferroni"
comparisonInfo$RIN = sampleInfo$RIN
for (cn in coln) {
 svaTmp = datExprRg[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprRg)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(as.matrix(datExprRg[,kp]),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOut[[cn]] = svaTmp
 isPN = svaOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   1 
# HIP TCx PCx FWM 
#   0   0   0   0 
signif(svaOut[[1]][apply(svaOut[[1]],1,min)<pThresh,],2)
#         HIP TCx PCx   FWM
# CBFA2T3   1   1   1 0.047


##### ---------------------------------------------------------------------------------------------
print("===== SVA in web site normalized data, Bonferroni corrected")
coln    = c("act_demented","ADorControl")
svaOut  = list() 
corMet  = "bonferroni"
comparisonInfo$RIN = sampleInfo$RIN
for (cn in coln) {
 svaTmp = datExprNp[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprNp)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(as.matrix(datExprNp[,kp]),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOut[[cn]] = svaTmp
 isPN = svaOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   0 
# HIP TCx PCx FWM 
#   0   0   0   0 


##### ---------------------------------------------------------------------------------------------
print("===== SVA in non RIN-corrected data, but for a subset of AD and control donors that are RIN-matched")

print("Read in data from the RIN, sex, and dementia matched cohort, then repeat DEX analysis on data uncorrected for RIN.")
matched  = scan(paste(extraFolder,"RIN_sex_dementia_matched_cohort.txt",sep=""),what="character",sep="\n")
isMatched= is.element(sampleInfo$donor_name,matched)&(sampleRIN>4.5)
datExprM = datExprTb2[,isMatched] 
comparisonInfoM = comparisonInfo[isMatched,]
regionM  = region[isMatched]

print("Identify genes significantly associated with all of the biological variables in each region using SVA as above.")
coln    = c("act_demented","ADorControl")
svaOut  = list() 
corMet  = "bonferroni"
for (cn in coln) {
 svaTmp = datExprM[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (regionM==r)&(!is.na(comparisonInfoM[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfoM[,cn])),comparisonInfoM[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprM)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(as.matrix(datExprM[,kp]),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOut[[cn]] = svaTmp
 isPN = svaOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   2   2   0   2 
# HIP TCx PCx FWM 
#   1   0   0   1 
signif(svaOut[[1]][apply(svaOut[[1]],1,min)<pThresh,],2)
#                HIP    TCx PCx   FWM
# ASPH         1.000 0.0041   1 1.000
# CHST6        0.014 1.0000   1 1.000
# DDX56        1.000 0.0440   1 1.000
# LOC105378049 0.033 1.0000   1 1.000
# TRBC2        1.000 1.0000   1 0.019
# VWF          1.000 1.0000   1 0.019


##### ---------------------------------------------------------------------------------------------
print("=====Student T.test in normalized data, Bonferroni corrected")

ttOut  = list() 
t.test.p = function(x,i1,i2){
 tt = t.test(x[i1],x[i2])
 return(tt$p.value)
} 
corMet  = "bonferroni"
for (cn in coln) {
 svaTmp = datExprRg[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  pValuesBatch    = apply(as.matrix(datExprRg[,kp]),1,t.test.p,siTmp[kp,2]==1,siTmp[kp,2]==2)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 ttOut[[cn]] = svaTmp
 isPN = ttOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   0 
# HIP TCx PCx FWM 
#   0   0   0   0 


##### ---------------------------------------------------------------------------------------------
print("=====Student T.test in normalized data, Bonferroni corrected")

ttOut  = list() 
t.test.p = function(x,i1,i2){
 tt = t.test(x[i1],x[i2])
 return(tt$p.value)
} 
corMet  = "bonferroni"
for (cn in coln) {
 svaTmp = datExprRg[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  pValuesBatch    = apply(as.matrix(datExprRg[,kp]),1,t.test.p,siTmp[kp,2]==1,siTmp[kp,2]==2)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 ttOut[[cn]] = svaTmp
 isPN = ttOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   0 
# HIP TCx PCx FWM 
#   0   0   0   0 


##### ---------------------------------------------------------------------------------------------
print("===== Limma in normalized data, Bonferroni corrected")
coln    = c("act_demented","ADorControl")
svaOut  = list() 
corMet  = "bonferroni"
comparisonInfo$RIN = sampleInfo$RIN
for (cn in coln) {
 svaTmp = datExprRg[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp    = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprRg)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~0+factor(siTmp[kp,2]))
  colnames(modBatch) = c("a","b")
  cont.matrix    <- makeContrasts(ab="a-b",levels=modBatch)
  fit   <- lmFit(as.matrix(datExprRg[,kp]),modBatch) 
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  pValuesBatch    = fit2$p.value[,"ab"]
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOut[[cn]] = svaTmp
 isPN = svaOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   1 
# HIP TCx PCx FWM 
#   0   0   0   0 
signif(svaOut[[1]][apply(svaOut[[1]],1,min)<pThresh,],2)
#         HIP TCx PCx   FWM
# CBFA2T3   1   1   1 0.039


##### ---------------------------------------------------------------------------------------------
print("===== ANOVA in normalized data, Bonferroni corrected")
coln    = c("act_demented","ADorControl")
svaOut  = list() 
corMet  = "bonferroni"
comparisonInfo$RIN = sampleInfo$RIN
for (cn in coln) {
 svaTmp = datExprRg[,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  pValuesBatch    = apply(as.matrix(datExprRg[,kp]),1,getAnovaPvalforApply,siTmp[kp,2])
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 ttOut[[cn]] = svaTmp
 isPN = ttOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   1 
# HIP TCx PCx FWM 
#   0   0   0   0 
signif(ttOut[[1]][apply(ttOut[[1]],1,min)<pThresh,],2)
#         HIP TCx PCx   FWM
# CBFA2T3   1   1   1 0.047


##### ---------------------------------------------------------------------------------------------
print("===== SVA on top 25 PCs in normalized data, Bonferroni corrected")
coln    = c("act_demented","ADorControl")
svaOut  = list() 
corMet  = "bonferroni"
numPCs  = 25
comparisonInfo$RIN = sampleInfo$RIN
for (cn in coln) {
 svaTmp = datExprRg[1:numPCs,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  datUse = datExprRg[,kp]
  pcaUse = prcomp(t(datUse))
  varExp = (floor(1000*(pcaUse$sdev)^2 / sum(pcaUse$sdev^2))/10)[1:numPCs]
  usePCs = pcaUse$x[,1:numPCs]
  rownames(svaTmp) <- colnames(usePCs)
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprRg)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(t(as.matrix(usePCs)),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOut[[cn]] = svaTmp
 isPN = svaOut[[cn]]<pThresh;  count = colSums(isPN)
 print(count)
}
# HIP TCx PCx FWM 
#   0   0   0   0 
# HIP TCx PCx FWM 
#   0   0   0   0 

##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------

print("===== Repeat SVA on top 25 PCs for both normalized and data, and for all of the various quantitative metrics, Bonferroni corrected")
coln    = colnames(comparisonInfo)
svaOutN = list() 
corMet  = "bonferroni"
numPCs  = 25
outN    = NULL
for (cn in coln) {
 svaTmp = datExprRg[1:numPCs,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  datUse = datExprRg[,kp]
  pcaUse = prcomp(t(datUse))
  varExp = (floor(1000*(pcaUse$sdev)^2 / sum(pcaUse$sdev^2))/10)[1:numPCs]
  usePCs = pcaUse$x[,1:numPCs]
  rownames(svaTmp) <- colnames(usePCs)
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprRg)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(t(as.matrix(usePCs)),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOutN[[cn]] = svaTmp
 isPN = svaOutN[[cn]]<pThresh;  count = colSums(isPN)
 outN = rbind(outN,count)
}

svaOutU = list() 
outU    = NULL
for (cn in coln) {
 svaTmp = datExprTb2[1:numPCs,1:length(regions)]*0
 colnames(svaTmp) = regions
 for (r in regions){
  kp = (region==r)&(!is.na(comparisonInfo[,cn]))
  datUse = datExprTb2[,kp]
  pcaUse = prcomp(t(datUse))
  varExp = (floor(1000*(pcaUse$sdev)^2 / sum(pcaUse$sdev^2))/10)[1:numPCs]
  usePCs = pcaUse$x[,1:numPCs]
  rownames(svaTmp) <- colnames(usePCs)
  siTmp = cbind(rep(1,length(comparisonInfo[,cn])),comparisonInfo[,cn])
  colnames(siTmp) = c("tmp","VAR")
  rownames(siTmp) = colnames(datExprTb2)
  siTmp = as.data.frame(siTmp)
  modBatch        = model.matrix(~VAR,data=siTmp[kp,])
  mod0Batch       = model.matrix(~1,data=siTmp[kp,])
  pValuesBatch    = f.pvalue(t(as.matrix(usePCs)),modBatch,mod0Batch)
  qValuesBatch    = p.adjust(pValuesBatch,method=corMet) 
  qValuesBatch[is.na(qValuesBatch)] = 1
  qValuesBatch[qValuesBatch>1]      = 1
  svaTmp[,r]     = qValuesBatch
 }
 svaOutU[[cn]] = svaTmp
 isPN = svaOutU[[cn]]<pThresh;  count = colSums(isPN)
 outU = rbind(outU,count)
}

out = cbind(outU,outN)
rownames(out) = coln
out[rowSums(out)>1,]
#                          HIP TCx PCx FWM HIP TCx PCx FWM
#Batch                       2   2   1   1   3   3   1   2
#sex                         4   3   3   3   3   1   3   4
#act_demented                1   1   1   1   0   0   0   0
#dsm_iv_clinical_diagnosis   0   0   0   2   0   1   0   0
#braak                       1   1   0   0   0   0   0   0
#ihc_at8                     1   0   0   0   1   0   0   0
#ihc_tau2_ffpe               1   0   0   0   1   0   0   0
#ihc_a_beta_ffpe             0   0   0   1   0   0   0   1
#ihc_iba1_ffpe               1   1   1   0   1   0   0   0
#ptau_ng_per_mg              0   1   0   1   0   0   0   1
#a_syn_pg_per_mg             0   0   0   1   0   0   0   1
#mcp_1_pg_per_mg             0   0   0   0   1   1   0   1
#il_10_pg_per_mg             1   0   0   0   2   0   0   0
#il_6_pg_per_mg              1   2   2   2   1   1   2   1
#mip_1a_pg_per_mg            0   1   1   0   0   0   0   0
#ADorControl                 0   1   1   1   0   0   0   0
#RIN                         1   1   1   1   1   1   3   4






##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------
##### ---------------------------------------------------------------------------------------------

print("Comparing dementia status with RNA quality in several previous data sets, for table 2, matching temporal cortex when possible.")

print("=====First compare against the MSBB data set.")

mc = read.csv(paste0(extraFolder,"MSBB_clinical.csv"))
mr = read.csv(paste0(extraFolder,"MSBB_RNAseq_covariates.csv"))
mm = cbind(mc[match(mr$individualIdentifier,mc$individualIdentifier),],mr)
mAge = as.character(mm$AOD)
mAge = gsub("\\+","",mAge)
mAge = as.numeric(mAge)
mAge[is.na(mAge)] = 0
kp36 = (mAge>=0)&grepl("accepted_hits",mm$fileName)&(mm$BrodmannArea=="BM36")
kp22 = (mAge>=0)&grepl("accepted_hits",mm$fileName)&(mm$BrodmannArea=="BM22")

pdf("zNotAFigure_MSBB_vs_RIN_PMI.pdf",height=8,width=10)
par(mfrow=c(2,2))
verboseScatterplot(jitter(mm$CDR[kp36]),mm$RIN[kp36],main="MSBB, Brodmann 36",xlab="CDR",ylab="RIN",pch=19,cex=0.5)
verboseScatterplot(jitter(mm$CDR[kp22]),mm$RIN[kp22],main="MSBB, Brodmann 22",xlab="CDR",ylab="RIN",pch=19,cex=0.5)
verboseBarplot(mm$RIN[kp36],mm$CDR[kp36],main="MSBB, Brodmann 36",xlab="CDR",ylab="RIN")
verboseBarplot(mm$RIN[kp22],mm$CDR[kp22],main="MSBB, Brodmann 22",xlab="CDR",ylab="RIN")
verboseScatterplot(jitter(mm$CDR[kp36]),mm$PMI[kp36]/60,main="MSBB, Brodmann 36",xlab="CDR",ylab="PMI",pch=19,cex=0.5)
verboseScatterplot(jitter(mm$CDR[kp22]),mm$PMI[kp22]/60,main="MSBB, Brodmann 22",xlab="CDR",ylab="PMI",pch=19,cex=0.5)
verboseBarplot(mm$PMI[kp36]/60,mm$CDR[kp36],main="MSBB, Brodmann 36",xlab="CDR",ylab="PMI")
verboseBarplot(mm$PMI[kp22]/60,mm$CDR[kp22],main="MSBB, Brodmann 22",xlab="CDR",ylab="PMI")
dev.off()
# There is a borderline significant decrease in RIN with increasing CDR (e.g., more severe AD) in this study

ctr = mm$RIN[kp36][mm$CDR[kp36]<2]
length(ctr) # 98
ctr = ctr[!is.na(ctr)]
adr = mm$RIN[kp36][mm$CDR[kp36]>=2]
length(adr) # 169
adr = adr[!is.na(adr)]
c(mean(ctr),sd(ctr),mean(adr),sd(adr))
# [1] 6.513 1.298 5.863 1.592
t.test(ctr,adr)$p.value
# [1] 0.001094


##### ---------------------------------------------------------------------------------------------
print("=====Next compare against the Mayo RNA-Seq data set.")

md    = read.csv(paste0(extraFolder,"MayoRNAseq_RNAseq_TCX_covariates.csv"),sep="\t")
kpMBB = md$Source=="MayoBrainBank_Dickson"
ctr   = md$RIN[(md$Diagnosis=="Control")&kpMBB]
length(ctr) # 33
ctr   = ctr[!is.na(ctr)]
adr   = md$RIN[(md$Diagnosis=="AD")&kpMBB]
length(adr) # 84
adr   = adr[!is.na(adr)]
c(mean(ctr),sd(ctr),mean(adr),sd(adr))
# [1] 8.5865854 0.5485699 7.6419355 1.2184072
t.test(ctr,adr)$p.value  
# [1] 0.0001982523

pdf("zNotAFigure_MayoRNAseq_vs_RIN_PMI.pdf",height=8,width=10)
par(mfrow=c(2,2))
dg = (md$Diagnosis=="AD")+2*(md$Diagnosis=="PSP")+3*(md$Diagnosis=="Pathologic Aging")
verboseScatterplot(jitter(dg),md$RIN,main="Mayo RNA-seq (TCx)",xlab="",ylab="RIN",pch=19,cex=0.5)
verboseScatterplot(jitter(dg),md$PMI,main="Mayo RNA-seq (TCx)",xlab="",ylab="PMI",pch=19,cex=0.5)
verboseBarplot(md$RIN,factor(md$Diagnosis,levels=levels(md$Diagnosis)[c(2,1,4,3)]),main="Mayo RNA-seq (TCx)",xlab="Diagnosis",ylab="RIN")
verboseBarplot(md$PMI,factor(md$Diagnosis,levels=levels(md$Diagnosis)[c(2,1,4,3)]),main="Mayo RNA-seq (TCx)",xlab="Diagnosis",ylab="PMI")
dev.off()
# This data set shows a highly significant INCREASE in RIN for AD donors.


##### ---------------------------------------------------------------------------------------------
print("=====Now repeat the analysis in our data set using a consistent statistical test, and only in temporal cortex.")

ctr = sampleRIN[(DementedYesNo=="N")&(region=="TCx")]
length(ctr)
ctr = ctr[!is.na(ctr)]
adr = sampleRIN[(comparisonInfo$dsm_iv_clinical_diagnosis=="Alzheimer's Disease Type")&(region=="TCx")]
length(adr)
adr = adr[!is.na(adr)]
c(mean(ctr),sd(ctr),mean(adr),sd(adr))
# [1] 6.8660000 0.9219788 6.1827586 1.2714078
t.test(ctr,adr)$p.value  
# [1] 0.01483126


##### ---------------------------------------------------------------------------------------------
print("=====Finally, repeat this analysis looking at the ROS and MAP data sets.")

mr = read.csv(paste0(extraFolder,"ROS_MAP_RNAquality.csv"))

# ROS
ctr = mr$rin_num[(mr$study=="ROS")&(mr$cogdxstr=="NCI")]
length(ctr)  # 107
ctr = ctr[!is.na(ctr)]
adr = mr$rin_num[(mr$study=="ROS")&(mr$cogdxstr=="Dementia-AD")]
length(adr)  # 136
adr = adr[!is.na(adr)]
c(mean(ctr),sd(ctr),mean(adr),sd(adr))
# [1] 7.2073895 1.0148662 7.0550668 0.9516606
t.test(ctr,adr)$p.value  
# [1] 0.2339039

# MAP
ctr = mr$rin_num[(mr$study=="MAP")&(mr$cogdxstr=="NCI")]
length(ctr)  # 94
ctr = ctr[!is.na(ctr)]
adr = mr$rin_num[(mr$study=="MAP")&(mr$cogdxstr=="Dementia-AD")]
length(adr)  # 120
adr = adr[!is.na(adr)]
c(mean(ctr),sd(ctr),mean(adr),sd(adr))
# [1] 7.2517859 1.0715720 6.7868043 0.9580844
t.test(ctr,adr)$p.value  
# [1] 0.001180904

# No signifincant difference in PMI in ROS or MAP (not shown).  Trend toward lower PMI in AD cases.  No difference in pH as well, although ~60% of the data points are missing.

