library(glmnet)
library(pROC)
cancer="GBM"
#read immune genesets associated with disease-free survival
survgbm <- read.table(paste("example_data\\",cancer,"_Immune_Survival_significant.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

#pull out genes in these genesets
library(GSEABase)
gene <- read.table("example_data\\human_gene_symbol_geneid.txt", sep=",", header=TRUE, stringsAsFactors=FALSE)
gene <- gene[,-1]
colnames(gene) <- c("Gene","geneId")
geneSets <- getGmt("example_data\\cibersort_LM22_LM7_ImSigNano.gmt", collectionType=BroadCollection(category="c7"), geneIdType=SymbolIdentifier())
survdatad <- unique(as.character(survgbm[,1]))
genesup <- ""
for(i in 1:length(survdatad)){
set <- survdatad[i]
for(j in 1:length(geneSets)){
if(geneSets[[j]]@setName==set){
genesup <- c(genesup,geneSets[[j]]@geneIds)
}
}
}
par3d <- dplyr::inner_join(data.frame(Gene=genesup[-1]),gene)
survgenegbm <- unique(par3d$Gene)

#select immune genes important for patient's disease-free survival in a specific cancer type using elastic net model
mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
datac <- cgdsr::getClinicalData(mycgds,mycaselist)
datac <- data.frame(Patient_ID=rownames(datac),datac, stringsAsFactors=FALSE)
clinical <- subset(datac, select=c("Patient_ID","DFS_STATUS","DFS_MONTHS"))
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- subset(clinical, DFS_MONTHS!="")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- 'NA')
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 0)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 1)
dataexpr <- data.frame(Patient_ID=colnames(datamat),t(datamat))
dataexpr <- subset(dataexpr, select=colnames(dataexpr) %in% c("Patient_ID",survgenegbm))
clindataexp <- dplyr::inner_join(clinical[,c(1,3)],dataexpr)
rownames(clindataexp) <- clindataexp[,1]
clindataexp <- clindataexp[,-1]
clindataexp[is.na(clindataexp)] <- 0
clindataexp[,1] <- as.numeric(clindataexp[,1])

clindataexp$DFS_MONTHS[clindataexp$DFS_MONTHS < 6] <- 0
clindataexp$DFS_MONTHS[clindataexp$DFS_MONTHS >= 12] <- 1
clindataexp <- subset(clindataexp, DFS_MONTHS==0 | DFS_MONTHS==1)
n <- nrow(clindataexp)
sample <- sample(seq(n), size = n * 0.6, replace = FALSE)
trainexp <- clindataexp[sample,]
testexp <- clindataexp[-sample,]
surv1exp <- subset(trainexp, DFS_MONTHS==1)
surv0exp <- subset(trainexp, DFS_MONTHS==0)
n1 <- nrow(surv1exp)
n0 <- nrow(surv0exp)
sample1 <- sample(seq(n0), size = n1, replace = FALSE)
surv0expdata <- surv0exp[sample1,]
surv0testdata <- surv0exp[-sample1,]
trainexp <- rbind(surv1exp,surv0expdata)
testexp <- rbind(testexp,surv0testdata)
a <- seq(0.1, 0.9, 0.05)
#Optimize lambda and alpha hyperparameters by 10-fold cross-validation
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- glmnet::cv.glmnet(as.matrix(clindataexp[,2:ncol(clindataexp)]), as.numeric(clindataexp[,1]), family = "gaussian", nfold = 10, type.measure = "deviance", parallel = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(as.matrix(clindataexp[,2:ncol(clindataexp)]), as.numeric(clindataexp[,1]), family = "gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
pathselect <- coef(md3)
pathpargbm <- as.matrix(pathselect)
write.table(pathpargbm,paste("example_data\\ElNet_model_",cancer,"_immunegenes_selected.txt", sep=""), sep="\t", quote=FALSE)
selectlincgbm <- unique(rownames(subset(as.data.frame(pathpargbm), s0!=0))[-1])
predicttest <- as.numeric(predict(md3, as.matrix(testexp[,2:ncol(testexp)]), type = "class")[,1])
classtest <- as.numeric(testexp[,1])
elnetroc <- roc(as.numeric(classtest), as.numeric(predicttest))
elnetauc <- elnetroc$auc

#calculate cancer-specific immune score using genes selected by the elastic net model
cancer <- "GBM"
ridgecoeff <- read.table(paste("example_data\\ElNet_model_",cancer,"_immunegenes_selected.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
ridgecoeff <- subset(ridgecoeff, s0 != 0)
datamat <- read.table(paste("example_data\\",cancer,"_data_expression_median_Zscores.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
datamat <- datamat[,-2]
colnames(datamat)[1] <- "Gene"
datamat <- subset(datamat, Gene %in% ridgecoeff[,1])
datamatridge <- dplyr::inner_join(ridgecoeff,datamat)
datakm <- read.table(paste("example_data\\DFS_KMpvalues_",cancer,".txt",sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
datakm <- subset(datakm, DFS_pval<0.05)
datamatridge <- subset(datamatridge, Gene %in% datakm[,1])

immuneScore <- numeric(ncol(datamatridge)-2)
for(i in 3:ncol(datamatridge)){
#immuneScore[i-2] <- 0
for(j in 1:nrow(datamatridge)){
immuneScore[i-2] <- immuneScore[i-2] + (as.numeric(datamatridge[j,i]) * as.numeric(datamatridge[j,2]))
}
}
imscoregbm <- data.frame(PATIENT_ID=colnames(datamatridge)[3:ncol(datamatridge)], immunescore=(scale(immuneScore)[,1]), stringsAsFactors=FALSE)

#check the survival significance of the calculated immune score for patients stratified by high (>median) or low (<median) immune score
imscoregbm[,1] <- as.character(imscoregbm[,1])
for(i in 1:nrow(imscoregbm)){
alls <- unlist(strsplit(imscoregbm[i,1],"[.]"))
imscoregbm[i,1] <- paste(alls[1],alls[2],alls[3],sep="-")
}
clin <- read.table(paste("example_data\\",cancer,"_clinical_data.txt", sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
clin <- subset(clin, select=c("PATIENT_ID","OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS"))
clin <- within(clin, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- 'NA')
clin$DFS_MONTHS <- as.numeric(clin$DFS_MONTHS)
clinclust <- dplyr::inner_join(imscoregbm,clin)
clinical <- clinclust
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, DFS_STATUS!="NA")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 0)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 1)
imcut <- character(nrow(clinical))
for(k in 1:nrow(clinical)){
if(clinical$immunescore[k] > median(clinical$immunescore))
imcut[k] <- "high"
if(clinical$immunescore[k] < median(clinical$immunescore))
imcut[k] <- "low"
}
clinical <- data.frame(clinical, ImmuneScore_GBM=imcut)
clinical <- subset(clinical, ImmuneScore_GBM=="high" | ImmuneScore_GBM=="low")
diffpi <- survival::survdiff(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ ImmuneScore_GBM, data=clinical)
pvalpi <- pchisq(diffpi$chisq, length(diffpi$n)-1, lower.tail=FALSE)
if(pvalpi < 0.0001)
pvalpi <- "p <0.0001"
mutsurv <- survival::survfit(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ ImmuneScore_GBM, data=clinical)
p <- survminer::ggsurvplot(mutsurv, pval=pvalpi, risk.table=TRUE, pval.size=12, font.tickslab=c(18,"plain","black"))
cancer <- "GBM"
png(paste("example_data\\DFSgenes_Array_median_ImmuneScore_ElNet_",cancer,"_KM.png",sep=""), width=4600, height=3000, res=600)
par(mar=c(4,3,4,6)+0.1)
print(p)
dev.off()
#OS
clin <- within(clin, OS_MONTHS[OS_MONTHS == '[Not Available]'] <- 'NA')
clin$OS_MONTHS <- as.numeric(clin$OS_MONTHS)
clinclust <- dplyr::inner_join(imscoregbm,clin)
clinical <- clinclust
clinical <- subset(clinical, OS_STATUS!="")
clinical <- within(clinical, OS_STATUS[OS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, OS_STATUS!="NA")
clinical <- within(clinical, OS_STATUS[OS_STATUS == 'LIVING'] <- 0)
clinical <- within(clinical, OS_STATUS[OS_STATUS == 'DECEASED'] <- 1)
imcut <- character(nrow(clinical))
for(k in 1:nrow(clinical)){
if(clinical$immunescore[k] > median(clinical$immunescore))
imcut[k] <- "high"
if(clinical$immunescore[k] < median(clinical$immunescore))
imcut[k] <- "low"
}
clinical <- data.frame(clinical, ImmuneScore_GBM=imcut)
clinical <- subset(clinical, ImmuneScore_GBM=="high" | ImmuneScore_GBM=="low")
diffpi <- survival::survdiff(survival::Surv(as.numeric(OS_MONTHS), as.numeric(OS_STATUS)) ~ ImmuneScore_GBM, data=clinical)
pvalpi <- pchisq(diffpi$chisq, length(diffpi$n)-1, lower.tail=FALSE)
if(pvalpi < 0.0001)
pvalpi <- "p <0.0001"
mutsurv <- survival::survfit(survival::Surv(as.numeric(OS_MONTHS), as.numeric(OS_STATUS)) ~ ImmuneScore_GBM, data=clinical)
p <- survminer::ggsurvplot(mutsurv, pval=pvalpi, risk.table=TRUE, pval.size=12, font.tickslab=c(18,"plain","black"))
cancer <- "GBM"
png(paste("example_data\\OSgenes_Array_median_ImmuneScore_ElNet_",cancer,"_KM.png",sep=""), width=4600, height=3000, res=600)
par(mar=c(4,3,4,6)+0.1)
print(p)
dev.off()
