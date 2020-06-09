cancer="GBM"
library("GSVA")
library(GSEABase)
geneSets <- getGmt("example_data\\cibersort_LM22_LM7_ImSigNano.gmt", collectionType=BroadCollection(category="c7"), geneIdType=SymbolIdentifier())
pathmatex <- read.table(paste("example_data\\",cancer,"_data_expression_median_Zscores.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
pathmatex <- pathmatex[,-2]
rownames(pathmatex) <- pathmatex[,1]
pathmatex <- pathmatex[,-1]

genesetsample <- gsva(as.matrix(pathmatex), geneSets, method="ssgsea", abs.ranking=FALSE, ssgsea.norm=TRUE, verbose=TRUE)
write.table(genesetsample,paste("example_data\\",cancer,"_TCGA_pooledimmuneset_ssGSEAnew.txt"), sep="\t", quote=FALSE)

datamatclin <- read.table(paste("example_data\\",cancer,"_TCGA_pooledimmuneset_ssGSEAnew.txt", sep=""), sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
datamat <- na.omit(t(datamatclin))
for(j in 1:nrow(datamat)){
alls <- unlist(strsplit(rownames(datamat)[j],"[.]"))
rownames(datamat)[j] <- paste(alls[1],alls[2],alls[3],sep="-")
}
survdata <- data.frame(GeneSet="",SurvPval=0)
clin <- read.table(paste("example_data\\",cancer,"_clinical_data.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
for(k in 1:ncol(datamat)){
if(nrow(datamat)>2){
gene2 <- colnames(datamat)[k]
gene2mutexp <- subset(datamat, select=colnames(datamat)[k])
geneup2 <- row.names(subset(gene2mutexp, gene2mutexp[,1]>=median(gene2mutexp[,1])))
genedown2 <- row.names(subset(gene2mutexp, gene2mutexp[,1]<median(gene2mutexp[,1])))
datac <- subset(clin, select=c("PATIENT_ID","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","AGE"))
rownames(datac) <- datac[,1]
datac <- datac[,-1]
if("DFS_MONTHS" %in% colnames(datac) && length(geneup2)>=2 && length(genedown2)>=2){
clinical <- subset(datac, select=c(DFS_STATUS,DFS_MONTHS,AGE))
mutcondition <- character(nrow(clinical))
for(i in 1:nrow(clinical)){
if(row.names(clinical)[i] %in% geneup2)
mutcondition[i] <- "down"
else if(row.names(clinical)[i] %in% genedown2)
mutcondition[i] <- "up"
else
mutcondition[i] <- "NA"
}
clinical <- data.frame(clinical, C=mutcondition, stringsAsFactors=FALSE)
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- '0')
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == 'NA'] <- '0')
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, DFS_STATUS!="NA")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 0)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 1)
clin1 <- subset(clinical, C=="down"|C=="up", select=c("DFS_MONTHS","DFS_STATUS","AGE","C"))
if(nrow(clin1)>2){
pb <- list()
blank <- grid::rectGrob(gp=grid::gpar(col="white"))
tryCatch(
{mutsurv2 <- survival::survfit(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ C, data=clin1)
pb[[2]] <- survminer::ggsurvplot(mutsurv2, data=clin1, palette=c("blue","red"), pval=TRUE, risk.table=TRUE, legend="none", xlab=NULL, ylab=NULL, font.tickslab=c(18,"plain","black"))
diffboth <- try(survival::survdiff(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ C, data=clin1))
pvalboth <- pchisq(diffboth$chisq, length(diffboth$n)-1, lower.tail=FALSE)
sdata <- data.frame(GeneSet=gene2,SurvPval=pvalboth)
survdata <- rbind(survdata,sdata)
png(paste("example_data\\DFS_","_",gene2,"_",cancer,"_Enrichment_Expression.png",sep=""), width=2000, height=1500, res=300)
par(mar=c(4,3,4,6)+0.1)
print(pb[[2]])
dev.off()
},
error=function(cond){
message(paste("Cannot plot Survival curves for",gene2,"in",cancer,"\n",cond,"\n"))
},
warning=function(cond){
message(paste("Warning: for",gene2,"in",cancer,"\n",cond,"\n"))
}
)
}
else{
message(paste("Cannot plot Survival curves for",gene2,"in",cancer,"\n","Not enough clinical data in mutation subgroups"))
}
}
else{
message(paste("Cannot plot Survival curves for","in",cancer,"\n","Not enough clinical data in mutation subgroups"))
}
}
}
survdata <- subset(survdata, SurvPval<0.05)
write.table(survdata,paste("example_data\\",cancer,"_Immune_Survival_significant.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

