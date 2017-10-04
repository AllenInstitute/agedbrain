#################################################################################################
##  These scripts take data from the the Aging, Dementia, and TBI website along with supplementary 
##  materials and uses it to perform all analyses for the TBI manuscript 
#################################################################################################

# R version 3.2.5 (2016-04-14) -- "Very, Very Secure Dishes" was used for the analysis, but the code should
# be compatible with most of the recent versions of R.

#################################################################################################
print("***** ------------------------------------------------------------------------")
print("***** Code #3: Produce figures and plots related to pathology for Figures 1-3.")
print("***** ------------------------------------------------------------------------")


#################################################################################################
print("Create box plots summarizing the donor info for presentation in Figure 1a.")
print("-----Note that text for this figure was added after the fact in Photoshop.")

isFirst <- function(x){
 out = rep(TRUE,length(x));
 for (i in 2:length(x)) if(is.element(x[i],x[1:(i-1)])) out[i]=FALSE
 return(out)
}

di2  = donorInfo[isFirst(donorInfo$name),]

pdf("Figure_1a_donorBarPlots.pdf")

# General demographic variables
x = table(di2$sex) 
 barplot(cbind(x,x),beside=FALSE, main="Gender",col=c("pink","blue"))
x = di2$ageNumeric; #x = table(x[x>0])
 par(mfrow=c(3,1)); hist(x,breaks=26,xlim=c(76.5,102.5),col="black",main="Age"); par(mfrow=c(1,1)); 
x = di2$education_years; #x = table(x)
  par(mfrow=c(3,1)); hist(x,breaks=15,xlim=c(6,21),col="black",main="Education (years)"); par(mfrow=c(1,1)); 

# Dementia / pathology variables
x = table(di2$act_demented);  
 barplot(cbind(x,x),beside=FALSE, main="ACT dementia diag",col=c("white","red"))
x = table(di2$dsm_iv_clinical_diagnosis);  
 barplot(cbind(x,x),beside=FALSE, main="DSM diag",col=c("white","red","red","red","grey","grey"))
x = table(di2$nincds_arda_diagnosis);  
 barplot(cbind(x,x),beside=FALSE, main="NINCDS diag",col=c(blueWhiteRed(6)[4:6],"grey"))
x = table(di2$braak);  
 barplot(cbind(x,x),beside=FALSE, main="Braak stage",col=c(blueWhiteRed(12)[6:12]))
x = table(di2$nia_reagan);  
 barplot(cbind(x,x),beside=FALSE, main="NIA Reagan",col=blueWhiteRed(8)[5:8])
x = table(di2$cerad);  
 barplot(cbind(x,x),beside=FALSE, main="CERAD",col=blueWhiteRed(8)[5:8])
x = table(di2$apo_e4_allele);  
 barplot(cbind(x,x),beside=FALSE, main="APO_E4",col=c("white","red","grey"))

# TBI variables
x = table(di2$num_tbi_w_loc);  
 barplot(cbind(x,x),beside=FALSE, main="Number of TBIs",col=blueWhiteRed(8)[5:8])
x = di2$age_at_first_tbi; x = x[x>0] #table(x[x>0])
 par(mfrow=c(3,1)); hist(x,breaks=18,xlim=c(0,90),col="black",main="Age at first TBI"); par(mfrow=c(1,1)); 
x = di2$longest_loc_duration;  x = table(x);  x = x[names(x)!="Unknown or N/A"]
 barplot(cbind(x,x),beside=FALSE, main="Longest LOC duration",col=blueWhiteRed(16)[10:16])

dev.off()


#################################################################################################
print("Summarize values above by dementia vs. control and TBI vs. non-TBI.")

numCn  <- c("ageNumeric","education_years","num_tbi_w_loc","age_at_first_tbi","braak","nia_reagan","cerad")
catCn  <- c("sex","apo_e4_allele","longest_loc_duration","act_demented","dsm_iv_clinical_diagnosis","nincds_arda_diagnosis")
isA    <- di2$act_demented=="Dementia"
isT    <- di2$num_tbi_w_loc>0
outTab <- catTab <- NULL

for (n in numCn){
 dat = as.numeric(di2[,n])
 outTab <- rbind(outTab,
   c(mean(dat[!isA]),sd(dat[!isA]),mean(dat[isA]),sd(dat[isA]),t.test(dat[!isA],dat[isA])$p.value,
   mean(dat[!isT]),sd(dat[!isT]),mean(dat[isT]),sd(dat[isT]),t.test(dat[!isT],dat[isT])$p.value))
}
colnames(outTab) <- c("Mean_NonDem","SD_NonDem","Mean_Dem","SD_Dem","PvalueDem","Mean_NonTBI","SD_NonTBI","Mean_TBI","SD_TBI","PvalueTBI")
rownames(outTab) <- numCn

for (n in catCn){
 dat = di2[,n]
 tmp <- cbind(table(dat),table(dat[!isA]),table(dat[isA]),table(dat[!isT]),table(dat[isT]))
 rownames(tmp) = paste(n,"-",names(table(dat)))
 catTab = rbind(catTab,tmp)
}
colnames(catTab) = c("Count_all","Count_NonDem","Count_Dem","Count_NonTBI","Count_TBI")

## Manual calculation of p-values to add to figure

phyper2 <- function (total, group1, group2, overlap, verySig=TRUE ,lt=TRUE){
  # This function is the same is phyper, just allows for more sensible input values
  q = overlap
  m = group1
  n = total-group1
  k = group2
  prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
  if (verySig) return(-prob)
  return(1-prob)
}

phyper2(106,50,43,23,FALSE)  # More females in dementia group relative to nondementia:  [1] 0.1011616
phyper2(106,56,43,22,FALSE)  # More females in TBI group relative to non TBI:           [1] 0.5345084
phyper2(106,50,20,13,FALSE)  # More >0 APOE in dementia group relative to nondementia:  [1] 0.02116414
phyper2(106,56,20,10,FALSE)  # More >0 APOE in TBI group relative to non TBI:           [1] 0.5142519
phyper2(106,56,50,28,FALSE)  # More dementia donors in TBI group relative to non TBI:   [1] 0.2083164

write.csv(outTab,"DemographicGroupStatistics_numeric.csv")
write.csv(catTab,"DemographicGroupStatistics_categorical.csv")


#################################################################################################
print("Create plots showing distribution densities for IHC and Luminex data for Figure 2a and 3a.")

print("-----Subset the appropriate metric data.")
luminex = colnames(metricInfo)[14:33]
luminex = setdiff(luminex,c("ab42_over_ab40_ratio","isoprostane_pg_per_mg","ptau_over_tau_ratio"))
lumText = gsub("_pg_per_mg","",luminex)
lumText = gsub("_ng_per_mg","",lumText)
names(lumText) = luminex
ihc     = colnames(metricInfo)[5:13]
ihcText = gsub("ihc_","i.",ihc)
names(ihcText) = ihc
ilText = c(ihcText,lumText)
il = names(ilText)

metricInfo2 = metricInfo
metricInfo2$ptau_ng_per_mg = metricInfo2$ptau_ng_per_mg*1000
metricInfo2$tau_ng_per_mg  = metricInfo2$tau_ng_per_mg*1000

kpDementia=c("ihc_at8","ihc_at8_ffpe","ihc_tau2_ffpe","ihc_a_beta","ihc_a_beta_ffpe",
  "ptau_ng_per_mg","tau_ng_per_mg","ab42_pg_per_mg")
  # Pathology excluded from plot: "ab40_pg_per_mg","ihc_a_syn","a_syn_pg_per_mg","ihc_ptdp_43_ffpe"

kpInflam=c("ihc_gfap_ffpe","ihc_iba1_ffpe",
  "tnf_a_pg_per_mg","il_6_pg_per_mg","mip_1a_pg_per_mg",
  "il_7_pg_per_mg","mcp_1_pg_per_mg","rantes_pg_per_mg")
  # Pathology excluded from plot: "il_1b_pg_per_mg","il_10_pg_per_mg","il_4_pg_per_mg","ifn_g_pg_per_mg"
  
# Information for the plots
getBins <- function(minVal,testVector,interval=1) return(sum((testVector>=minVal)&(testVector<(minVal+interval))))
nbin    = 25   # The number of bins for the histogram plot (more = more accurate, but noisier)
outlier = 3    # Values more than 3 SD above mean are counted as outliers and not plotted on the line
  
print("-----Generate plots for pathology metrics.")

pdf("Figure_2a_densityPlotsByRegion_dementia_values.pdf",height=9,width=5)
par(mar=c(0,6,0,0))
par(mfrow=c(12,1)) 
yMax = c(0.09,0.09,0.02,0.07,0.07,5000,5000,650) 
names(yMax) = kpDementia
for (l in kpDementia){
 tmp  = metricInfo2[,l]
 tmp2 = (tmp-mean(tmp,na.rm=TRUE))/sd(tmp,na.rm=TRUE)
 cnt  = tmp > yMax[l] # tmp2 > outlier
 cnt[is.na(cnt)] = FALSE
 tmp2 = tmp
  qvl  = yMax[l]  # max(tmp2[!cnt],na.rm=TRUE)
 plot(0,0,axes=FALSE,lwd=8,xlab="",ylab=ilText[l],main="",cex.lab=2,ylim=c(0,1.2),xlim=c(0,qvl*1.2),col="white")
 text(qvl*1.125,0.3,signif(qvl,2),cex=1.5)
 for (r in regions){
  tmpr = tmp2[(region==r)&(!is.na(tmp2))]
  xqvl = qvl*(0:(nbin-1))/nbin
  yqvl = as.numeric(lapply(xqvl,getBins,tmpr,qvl/nbin))
  yqvl = (yqvl+yqvl[c(1,1:(nbin-1))]+yqvl[c(2:nbin,nbin)])+0.5*(yqvl[c(1,1,1:(nbin-2))]+yqvl[c(3:nbin,nbin,nbin)])
  yqvl = yqvl/max(yqvl)
  xqvl = xqvl+0.5*qvl/nbin  
  lines(xqvl,yqvl,lwd=5,col=regCols[r])
 }
 abline(v=c(0,qvl),lwd=2,col=c("black","grey"))
 abline(h=c(0,1.2),lwd=1,col="grey")
 for (r in regions){
  tmpr = tmp2[(region==r)&(!is.na(tmp2))] 
  points(median(tmpr),1,cex=2,col="black",bg=regCols[r],pch=25)
  text(qvl*(1+0.05*which(regions==r)),0.9,sum(cnt[region==r]),cex=2,col=regCols[r])
 }
}
dev.off()


pdf("Figure_3a_densityPlotsByRegion_inflammation_values.pdf",height=9,width=5)
par(mar=c(0,6,0,0))
par(mfrow=c(12,1))
yMax = c(rep(0.15,2),rep(100,6))
names(yMax) = kpInflam
for (l in kpInflam){
 tmp  = metricInfo[,l]
 tmp2 = (tmp-mean(tmp,na.rm=TRUE))/sd(tmp,na.rm=TRUE)
 cnt  = tmp > yMax[l] # tmp2 > outlier
 cnt[is.na(cnt)] = FALSE
 tmp2 = tmp
 qvl  = yMax[l]  # max(tmp2[!cnt],na.rm=TRUE)
 plot(0,0,axes=FALSE,lwd=8,xlab="",ylab=ilText[l],main="",cex.lab=2,ylim=c(0,1.2),xlim=c(0,qvl*1.2),col="white")
 text(qvl*1.125,0.3,signif(qvl,2),cex=1.5)
 for (r in regions){
  tmpr = tmp2[(region==r)&(!is.na(tmp2))]
  xqvl = qvl*(0:(nbin-1))/nbin
  yqvl = as.numeric(lapply(xqvl,getBins,tmpr,qvl/nbin))
  yqvl = (yqvl+yqvl[c(1,1:(nbin-1))]+yqvl[c(2:nbin,nbin)])+0.5*(yqvl[c(1,1,1:(nbin-2))]+yqvl[c(3:nbin,nbin,nbin)])
  yqvl = yqvl/max(yqvl)
  xqvl = xqvl+0.5*qvl/nbin  
  lines(xqvl,yqvl,lwd=5,col=regCols[r])
 }
 abline(v=c(0,qvl),lwd=2,col=c("black","grey"))
 abline(h=c(0,1.2),lwd=1,col="grey")
 for (r in regions){
  tmpr = tmp2[(region==r)&(!is.na(tmp2))] 
  points(median(tmpr),1,cex=2,col="black",bg=regCols[r],pch=25)
  text(qvl*(1+0.05*which(regions==r)),0.9,sum(cnt[region==r]),cex=2,col=regCols[r])
 }
}
dev.off()


#################################################################################################
print("Explore the relationships between age and dementia pathologies for Figure 2b.")

lcasi   = read.csv(paste(extraFolder,"LastCASI110.csv",sep=""),row.names=1)   # Information on cognitive status
lcasi   = lcasi[as.character(donorInfo$name),]
rownames(lcasi) = rownames(donorInfo)
abetaInfo$casi_irt = lcasi$casi_irt
abetaInfo$CASI_irt = c("a--","b-","c0","d+")[pmax(1,floor(lcasi$casi_irt+0.5)+3)]  # <-1.5, -1.5:-0.5, -0.5:0.5, >0.5
abetaInfo$Age = abetaInfo$ageNumeric
abetaInfo$act_demented2 = as.character(abetaInfo$act_demented)
abetaInfo$act_demented2[abetaInfo$act_demented2=="Dementia"] = "Yes Dementia"
abetaInfo$act_demented2 = factor(abetaInfo$act_demented2)
abetaInfo$dsm_iv_hasAD = NA
abetaInfo$dsm_iv_hasAD[donorInfo$dsm_iv_clinical_diagnosis=="No Dementia"] = "No AD"
abetaInfo$dsm_iv_hasAD[donorInfo$dsm_iv_clinical_diagnosis=="Alzheimer's Disease Type"] = "Yes AD"
abetaInfo$dsm_iv_hasAD = factor(abetaInfo$dsm_iv_hasAD)
if(length(intersect(colnames(lcasi),"Age"))>0) abetaInfo$Age = lcasi$Age
sn      = "Age"
cnsTmp  = c("cerad","braakCombine")
outInfo = abetaInfo[isFirst(donorInfo$name),c(sn,cnsTmp)]

getAnovaPvalforApply <- function(x,varLabels,varWeights=NULL){
  anovadat  = as.data.frame(cbind(varLabels,x))
  aov.out   = summary(aov(as.numeric(anovadat[,2])~anovadat[,1],data=anovadat,weights=varWeights))
  return(aov.out[[1]]$'Pr(>F)'[1])
}

pdf("Figure_2b_plotsForAgeVsPathologyMetrics.pdf",height=4,width=11)
grobb = list()
for (cn in cnsTmp) {
 kp  = (!is.na(outInfo[,cn]))&(!is.na(outInfo[,sn]))
 pval= getAnovaPvalforApply(outInfo[kp,sn],outInfo[kp,cn])
 ddf = data.frame(NUMS = outInfo[kp,sn], GRP = as.character(outInfo[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("p =",signif(pval,3))) + labs(x=cn,y=sn) +
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()
print("Note that the horizonal positions of dots within each bar plot is randomly chosen and may differ from run to run.")

pdf("SupFigure_XX_plotsForAgeVsIHCMetrics.pdf")
yl = "Age"
xl = "ihc_at8"
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="HIP")
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="Tau (AT8) in HIP",pch=19,xlab=xl,ylab=yl,col=regCols["HIP"])
xl = "ihc_a_beta"
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="PCx")
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="A-beta in PCx",pch=19,xlab=xl,ylab=yl,col=regCols["PCx"])
dev.off()


#################################################################################################
print("Repeat age-related comparisons for dementia and control groups separately.")

pdf("zNotAFigure_plotsForAgeVsPathologyMetrics_controlDementia.pdf",height=4,width=11)
sn      = "Age"
cnsTmp  = c("cerad","braakCombine")
outInfo = abetaInfo[isFirst(donorInfo$name),c(sn,cnsTmp)]
grobb   = list()
for (d in c("N","Y")) for (cn in cnsTmp) {
 kp  = (!is.na(outInfo[,cn]))&(!is.na(outInfo[,sn]))&(DementedYesNo[isFirst(donorInfo$name)]==d)
 pval= getAnovaPvalforApply(outInfo[kp,sn],outInfo[kp,cn])
 ddf = data.frame(NUMS = outInfo[kp,sn], GRP = as.character(outInfo[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("Dementia =",d,"/ p =",signif(pval,3))) + labs(x=cn,y=sn) +
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()

pdf("zNotAFigure_plotsForAgeVsIHCMetrics_controlDementia.pdf")
yl = "Age"
for (d in c("N","Y")) {
 xl = "ihc_at8"
 kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="HIP")&(DementedYesNo==d)
 verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="Tau (AT8) in HIP",pch=19,xlab=paste("ihc_at8, Dementia =",d),ylab=yl,col=regCols["HIP"])
 xl = "ihc_a_beta"
 kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="PCx")&(DementedYesNo==d)
 verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="A-beta in PCx",pch=19,xlab=paste("ihc_a_beta, Dementia =",d),ylab=yl,col=regCols["PCx"])
}
dev.off()



#################################################################################################
print("Explore the relationships between dementia/AD status and dementia pathologies for Figure 2d.")
cnsTmp  = c("act_demented2","dsm_iv_hasAD")
sns     = c("cerad","braak","ihc_at8","ihc_a_beta","ihc_a_beta_ffpe")
outInfo = abetaInfo[,c(sns,cnsTmp)]

pdf("Figure_2c_plotsForDementiaOrADVsPathologyMetrics.pdf",height=4,width=11)
grobb = list()
for (cn in cnsTmp) for (sn in sns) {
 kp  = (!is.na(outInfo[,cn]))&(!is.na(outInfo[,sn]))
 if(is.element(sn,c("cerad","braak"))) kp = kp&isFirst(donorInfo$name)
 if(is.element(sn,c("ihc_at8"))) kp = kp&(region=="HIP")
 if(is.element(sn,c("ihc_a_beta","ihc_a_beta_ffpe"))) kp = kp&(region=="PCx")
 pval= getAnovaPvalforApply(outInfo[kp,sn],outInfo[kp,cn])
 ddf = data.frame(NUMS = outInfo[kp,sn], GRP = as.character(outInfo[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("p =",signif(pval,3))) + labs(x=cn,y=sn) +
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()
print("Note that the horizonal positions of dots within each bar plot is randomly chosen and may differ from run to run.")

#################################################################################################
print("Explore the relationships between CASI and dementia pathologies for Figure 2d.")
sn      = "casi_irt"
cnsTmp  = c("cerad","braakCombine")
outInfo = abetaInfo[isFirst(donorInfo$name),c(sn,cnsTmp)]

pdf("Figure_2d_plotsForCASIVsPathologyMetrics.pdf",height=4,width=11)
grobb = list()
for (cn in cnsTmp) {
 kp  = (!is.na(outInfo[,cn]))&(!is.na(outInfo[,sn]))
 pval= getAnovaPvalforApply(outInfo[kp,sn],outInfo[kp,cn])
 ddf = data.frame(NUMS = outInfo[kp,sn], GRP = as.character(outInfo[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("p =",signif(pval,3))) + labs(x=cn,y=sn) +
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()
print("Note that the horizonal positions of dots within each bar plot is randomly chosen and may differ from run to run.")

pdf("Figure_2d_plotsForCASIVsIHCMetrics.pdf")
yl = "casi_irt"
xl = "ihc_at8"
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="HIP")
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="Tau (AT8) in HIP",pch=19,xlab=xl,ylab=yl,col=regCols["HIP"])
xl = "ihc_a_beta"
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="PCx")
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="A-beta in PCx",pch=19,xlab=xl,ylab=yl,col=regCols["PCx"])
dev.off()


#################################################################################################
print("Repeat the comparisons between pathology and dementia in <90 vs. 90+ age groups to see if correlations are weaker in oldest old.")

isOO    = abetaInfo$Age>=90
cnsTmp  = c("act_demented2","dsm_iv_hasAD")
sns     = c("cerad","braak","ihc_at8","ihc_a_beta","ihc_a_beta_ffpe")
outInfo = abetaInfo[,c(sns,cnsTmp)]

pdf("SupFigure_XX_plotsForDementiaOrADVsPathologyMetrics_ByAge.pdf",height=4,width=14)
grobb = list()
for (cn in cnsTmp) for (sn in sns) {
 kp  = (!is.na(outInfo[,cn]))&(!is.na(outInfo[,sn]))
 if(is.element(sn,c("cerad","braak"))) kp = kp&isFirst(donorInfo$name)
 if(is.element(sn,c("ihc_at8"))) kp = kp&(region=="HIP")
 if(is.element(sn,c("ihc_a_beta","ihc_a_beta_ffpe"))) kp = kp&(region=="PCx")
 yl  = range(outInfo[kp,sn]);
 kp  = kp&(!isOO)
 pval= getAnovaPvalforApply(outInfo[kp,sn],outInfo[kp,cn])
 ddf = data.frame(NUMS = outInfo[kp,sn], GRP = as.character(outInfo[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("<90: p =",signif(pval,3))) + labs(x=cn,y=sn) + 
 ylim(yl[1], yl[2]) + geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
for (cn in cnsTmp) for (sn in sns) {
 kp  = (!is.na(outInfo[,cn]))&(!is.na(outInfo[,sn]))
 if(is.element(sn,c("cerad","braak"))) kp = kp&isFirst(donorInfo$name)
 if(is.element(sn,c("ihc_at8"))) kp = kp&(region=="HIP")
 if(is.element(sn,c("ihc_a_beta","ihc_a_beta_ffpe"))) kp = kp&(region=="PCx")
 yl  = range(outInfo[kp,sn]) 
 kp  = kp&(isOO)
 pval= getAnovaPvalforApply(outInfo[kp,sn],outInfo[kp,cn])
 ddf = data.frame(NUMS = outInfo[kp,sn], GRP = as.character(outInfo[kp,cn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste("90+: p =",signif(pval,3))) + labs(x=cn,y=sn) +
 ylim(yl[1], yl[2]) + geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()

pdf("zNotAFigure_plotsForCASIVsIHCMetrics_ByAge.pdf")
yl = "casi_irt"
xl = "ihc_at8"
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="HIP")&(!isOO)
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="Tau (AT8) in HIP, <90yo",pch=19,xlab=xl,ylab=yl,col=regCols["HIP"])
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="HIP")&(isOO)
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="Tau (AT8) in HIP, 90+ yo",pch=19,xlab=xl,ylab=yl,col=regCols["HIP"])
xl = "ihc_a_beta"
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="PCx")&(!isOO)
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="A-beta in PCx, <90yo",pch=19,xlab=xl,ylab=yl,col=regCols["PCx"])
kp = (!is.na(abetaInfo[,yl]))&(!is.na(abetaInfo[,xl]))&(region=="PCx")&(isOO)
verboseScatterplot(abetaInfo[kp,xl],abetaInfo[kp,yl],main="A-beta in PCx, 90+ yo",pch=19,xlab=xl,ylab=yl,col=regCols["PCx"])
dev.off()
# No difference in CASI vs. pathology wrt age.


#################################################################################################
print("Show that there is consistant pathology metrics between IHC and Braak/CERAD for Supplemental Figure.")

pdf("SupFigure_XX_FFPEvsFreshFrozen_pathologyQuants.pdf")
verboseScatterplot(abetaInfo[,"ihc_at8_ffpe"],abetaInfo[,"ihc_at8"],main="Tau (AT8)",pch=19,xlab="FFPE",ylab="Fresh Frozen",col=regionColors)
verboseScatterplot(abetaInfo[,"ihc_a_beta_ffpe"],abetaInfo[,"ihc_a_beta"],main="A-beta",pch=19,xlab="FFPE",ylab="Fresh Frozen",col=regionColors)
for (i in 1:4) text(0.055,i*0.005,regions[i],col=regCols[regions[i]])
dev.off()

sn      = "braakCombine"
cnsTmp  = c("ihc_at8","ihc_at8_ffpe","ihc_tau2_ffpe","ptau_ng_per_mg","tau_ng_per_mg")
outInfo = abetaInfo[,c(sn,cnsTmp)]

pdf("SupFigure_XX_plotsForTauVsBraak.pdf",height=4,width=11)
grobb = list()
for (r in regions) for (cn in cnsTmp) {
 kp  = (!is.na(outInfo[,sn]))&(!is.na(outInfo[,cn]))&(region==r)
 pval= getAnovaPvalforApply(outInfo[kp,cn],outInfo[kp,sn])
 ddf = data.frame(NUMS = outInfo[kp,cn], GRP = as.character(outInfo[kp,sn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste(r,"- p =",signif(pval,3))) + labs(x=sn,y=cn) +
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()

sn      = "cerad"
cnsTmp  = c("ihc_a_beta","ihc_a_beta_ffpe","ab42_pg_per_mg","ab40_pg_per_mg","braak")
outInfo = abetaInfo[,c(sn,cnsTmp)]

pdf("SupFigure_XX_plotsForAbetaVsCERAD.pdf",height=4,width=11)
grobb = list()
for (r in regions) for (cn in cnsTmp) {
 kp  = (!is.na(outInfo[,sn]))&(!is.na(outInfo[,cn]))&(region==r)
 pval= getAnovaPvalforApply(outInfo[kp,cn],outInfo[kp,sn])
 ddf = data.frame(NUMS = outInfo[kp,cn], GRP = as.character(outInfo[kp,sn]))
 plt = ggplot(ddf, aes(x=GRP, y=NUMS)) + ggtitle(paste(r,"- p =",signif(pval,3))) + labs(x=sn,y=cn) +
 geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
 geom_jitter(position=position_jitter(width=.5, height=0))
 grobb[[length(grobb)+1]] = ggplotGrob(plt)
}
print(marrangeGrob(grobb,ncol=5,nrow=1))
dev.off()


#################################################################################################
print("Compare inflammation vs. dementia pathology metrics.")

kpDementia=c("ihc_at8","ihc_at8_ffpe","ihc_tau2_ffpe","ihc_a_beta","ihc_a_beta_ffpe",
  "ptau_ng_per_mg","tau_ng_per_mg","ab42_pg_per_mg")
  # Pathology excluded from plot: "ab40_pg_per_mg","ihc_a_syn","a_syn_pg_per_mg","ihc_ptdp_43_ffpe"
kpInflam=c("ihc_gfap_ffpe","ihc_iba1_ffpe",  "tnf_a_pg_per_mg","il_6_pg_per_mg","mip_1a_pg_per_mg",
  "il_7_pg_per_mg","mcp_1_pg_per_mg","rantes_pg_per_mg")
# Pathology excluded from plot: "il_1b_pg_per_mg","il_10_pg_per_mg","il_4_pg_per_mg","ifn_g_pg_per_mg"

pdf("Figure_4c_InflammationVsDementiaValues.pdf",width=10,height=10)
par(mfrow=c(2,2))
for (inf in kpInflam) for (dem in kpDementia) for (r in regions){
 kp = (region==r)&(!is.na(metricInfo2[,dem]))&(!is.na(metricInfo2[,inf]))
 verboseScatterplot(metricInfo2[kp,dem],metricInfo2[kp,inf],xlab=dem,ylab=inf,main=r,pch=19,col=regCols[r])
} 
dev.off()
## GFAP and Tau in Tcx only ==> No other relationships between pathology and inflammation
## The two outlier samples on the inflmaation panel C do not show differential regionality of disease metrics
## --> That said, the high microglia one also is higher in tau

