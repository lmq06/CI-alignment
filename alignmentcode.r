library(limma)
library(survival)
library(survminer)
library(tidyverse)
# read in RNAseq data

rna.data <- read.table("TCGA.PAAD.HiSeqV2", sep = "\t", header = TRUE, row.names = 1)
dim(rna.data)
rna.data[1:3, 1:3]
rna.data <- as.matrix(rna.data)

# read in survival data

surv.data <- read.table("PAAD_survival.txt", sep = "\t", header = T, row.names = 1)
rownames(surv.data) <- gsub(rownames(surv.data), pattern = "-", replace = ".")
surv.data[1:3, ]


# read in clinical data

clin.data <- read.table("PAAD_clinicalMatrix", sep = "\t", header = T, row.names = 1)
rownames(clin.data) <- gsub(rownames(clin.data), pattern = "-", replace = ".")
clin.data[1:3, ]

# remove normal sample
mete <- data.frame(colnames(rna.data))
for (i in 1:length(mete[, 1])) {
  num <- as.numeric(as.character(substring(mete[i, 1], 14, 15))) # The 14th and 15th digits represent whether it is a cancer sample, 01~09 are cancer samples
  if (num < 10) {
    mete[i, 2] <- "T"
  } else {
    mete[i, 2] <- "N"
  }
}

names(mete) <- c("id", "group")
mete$group <- as.factor(mete$group)
mete <- subset(mete, mete$group == "T")

mete <- mete[which(mete$id %in% rownames(surv.data)), ]

rna.data <- rna.data[, which(colnames(rna.data) %in% mete$id)]
surv.data <- surv.data[match(mete$id, rownames(surv.data)), ]
clin.data <- clin.data[match(mete$id, rownames(clin.data)), ]

# fit a multivariate model for each genes's expression, and find those most significantly associated to patient outcome

paad.os <- Surv(surv.data$OS.time, surv.data$OS)

# First find out which clinical data are associated with survival

age <- as.numeric(clin.data$age_at_initial_pathologic_diagnosis)

histology <- as.factor(clin.data$histological_type)

smoke <- clin.data$number_pack_years_smoked
table(clin.data$tobacco_smoking_history)
smoke[which(clin.data$tobacco_smoking_history == 1)] <- 0

smoke.cat <- rep(NA, nrow(clin.data))
smoke.cat[which(clin.data$tobacco_smoking_history == 1)] <- 0
smoke.cat[which(clin.data$tobacco_smoking_history > 1)] <- 1


drink <- rep(NA, nrow(clin.data))
drink[which(clin.data$alcohol_history_documented == "NO")] <- 0
drink[which(clin.data$alcohol_history_documented == "YES")] <- 1


summary(coxph(paad.os ~ age))$coef
summary(coxph(paad.os ~ histology))$coef
summary(coxph(paad.os ~ smoke))$coef
summary(coxph(paad.os ~ smoke.cat))$coef
summary(coxph(paad.os ~ drink))$coef

results.multivariate <- array(NA, c(nrow(rna.data), 4))

colnames(results.multivariate) <- c("HR", "LCI", "UCI", "PVAL")
rownames(results.multivariate) <- rownames(rna.data)

results.multivariate <- as.data.frame(results.multivariate)



for (i in 1:nrow(rna.data))

{
  coxphmodel <- coxph(paad.os ~ rna.data[i, ] + age)

  results.multivariate$HR[i] <- summary(coxphmodel)$coef[1, 2]

  results.multivariate$LCI[i] <- summary(coxphmodel)$conf.int[1, 3]

  results.multivariate$UCI[i] <- summary(coxphmodel)$conf.int[1, 4]

  results.multivariate$PVAL[i] <- summary(coxphmodel)$coef[1, 5]
}



results.multivariate <- as.data.frame(results.multivariate)

results.multivariate$FDR <- p.adjust(results.multivariate$PVAL, method = "fdr")



results.multivariate <- results.multivariate[order(results.multivariate$FDR, decreasing = F), ]

results.multivariate[1:10, ]


USP20.high <- as.numeric(rna.data["USP20", ] > median(rna.data["USP20", ]))
data1 <- data.frame(os.time=surv.data$OS.time, os=surv.data$OS.time, USP20=USP20.high)

fit <- survfit(paad.os ~USP20, data=data1)
png("PAAD_OS_byUSP20.png",width=6,height=6,units='in',res=300)
ggsurvplot(fit,
  pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  ggtheme = theme_bw(), # Change ggplot2 theme
  surv.median.line = "hv",
  palette = c("#E7B800", "#2E9FDF")
)
print(myplot)
dev.off()

#<U+7814><U+7A76> USP20

# Methylation

# getting the 450k annotation file
meth.annot <- read.table("illuminaMethyl450_hg19_GPL16304_TCGAlegacy", sep = "\t", head = FALSE)

# load in the methylation data
meth.data <- read.table("HumanMethylation450", sep = "\t", header = T)
rownames(meth.data) <- meth.data$sample
meth.data <- meth.data[, -1]

rna.data2 <- rna.data[, which(is.element(colnames(rna.data), colnames(meth.data)))]
meth.data2 <- meth.data[, which(is.element(colnames(meth.data), colnames(rna.data2)))]
surv.data2 <- surv.data[which(is.element(rownames(surv.data), colnames(rna.data2))), ]

rna.data2 <- as.matrix(rna.data2[, order(colnames(rna.data2))])
meth.data2 <- as.matrix(meth.data2[, order(colnames(meth.data2))])
surv.data2 <- as.data.frame(surv.data2[colnames(surv.data2), ])

meth.USP20 <- rownames(meth.annot[which(meth.annot$V2 == "USP20"), ])
meth.annot[meth.USP20, ]

meth.data.USP20 <- meth.data2[as.numeric(meth.USP20), ]
rownames(meth.data.USP20) <- meth.annot[meth.USP20, ]$V1
na.count <- apply(meth.data.USP20, 1, function(x) sum(as.numeric(is.na(x))))

exclude <- as.numeric(na.count > 0.5 * ncol(meth.data.USP20))

meth.data.USP20 <- meth.data.USP20[which(exclude == 0), ]



results <- array(NA, c(nrow(meth.data.USP20), 4))
rownames(results) <- rownames(meth.data.USP20)
colnames(results) <- c("Cor.USP20", "Cor.test.USP20", "Mean.high.USP20", "Mean.low.USP20")

USP20.high2 <- as.numeric(as.numeric(rna.data2["USP20", ]) > median(as.numeric(rna.data2["USP20", ])))

# Confirm whether the data obey a normal distribution
shapiro.test(as.numeric(meth.data.USP20[6, ]))

for (i in 1:nrow(meth.data.USP20))
{
  results[i, 1] <- cor.test(as.numeric(rna.data2["USP20", ]), as.numeric(meth.data.USP20[i, ]), use = "c")$est
  results[i, 2] <- cor.test(as.numeric(rna.data2["USP20", ]), as.numeric(meth.data.USP20[i, ]), use = "c")$p.value
}

results[, 3] <- apply(meth.data.USP20[, which(USP20.high2 == 1)], 1, mean, na.rm = T)
results[, 4] <- apply(meth.data.USP20[, which(USP20.high2 == 0)], 1, mean, na.rm = T)

results$FDR <- p.adjust(results$Cor.test.USP20, method = "fdr")
results <- results[order(results$FDR, decreasing = F), ]

# plot of methylation vs expression
png("PAAD_USP20_ExpVsMeth.png")
plot(as.numeric(meth.data2["cg02278574", ]), as.numeric(rna.data2["USP20", ]), xlab = "cg02278574", ylab = "USP20 RNAseq")
abline(lm(rna.data2["USP20", ] ~ meth.data2["cg02278574", ]))
dev.off()



#CNV
cnv.data <- read.delim("Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
cnv.USP20 <- cnv.data[which(cnv.data$Gene.Symbol=="USP20"),]
as.numeric(cnv.USP20)
table(cnv.USP20)

#mirna
mirna.data <- read_tsv("TCGA-PAAD.mirna.tsv")
mirna.data <- as.data.frame(mirna.data)
rownames(mirna.data) <- mirna.data$miRNA_ID
mirna.data <- mirna.data[,-1]

colnames(mirna.data) <- gsub(colnames(mirna.data), pattern="-", replace=".")
colnames(mirna.data) <- gsub('.{1}$', '', colnames(mirna.data))

rna.data2 <- rna.data[, which(is.element(colnames(rna.data), colnames(mirna.data)))]
mirna.data2 <- mirna.data[, which(is.element(colnames(mirna.data), colnames(rna.data2)))]


rna.data2 <- as.matrix(rna.data2[, order(colnames(rna.data2))])
mirna.data2 <- as.matrix(mirna.data2[, order(colnames(mirna.data2))])


results <- array(NA, c(nrow(mirna.data2), 4))
rownames(results) <- rownames(mirna.data2)
colnames(results) <- c("Cor.USP20", "Cor.test.USP20", "Mean.high.USP20", "Mean.low.USP20")

USP20.high2 <- as.numeric(as.numeric(rna.data2["USP20", ]) > median(as.numeric(rna.data2["USP20", ])))

# Confirm whether the data obey a normal distribution
shapiro.test(as.numeric(meth.data.USP20[6, ]))

for (i in 1:nrow(mirna.data2))
{
  results[i, 1] <- cor.test(as.numeric(rna.data2["USP20", ]), as.numeric(mirna.data2[i, ]), use = "c")$est
  results[i, 2] <- cor.test(as.numeric(rna.data2["USP20", ]), as.numeric(mirna.data2[i, ]), use = "c")$p.value
}

results[, 3] <- apply(mirna.data2[, which(USP20.high2 == 1)], 1, mean, na.rm = T)
results[, 4] <- apply(mirna.data2[, which(USP20.high2 == 0)], 1, mean, na.rm = T)

results<-as.data.frame(results)

results$FDR<-p.adjust(results$Cor.test.USP20,method="fdr")
results<-results[order(results$FDR, decreasing=F),]
results <- na.omit(results)
results[1:100,]


#Prognostic biomarker related differentially expressed genes
rna.data <- read.table("TCGA.PAAD.HiSeqV2", sep = "\t", header = TRUE, row.names = 1)
dim(rna.data)
rna.data[1:3, 1:3]
rna.data <- as.matrix(rna.data)

surv.data <- read.table("PAAD_survival.txt", sep = "\t", header = T, row.names = 1)
rownames(surv.data) <- gsub(rownames(surv.data), pattern = "-", replace = ".")
surv.data[1:3, ]

mete <- data.frame(colnames(rna.data))
for (i in 1:length(mete[, 1])) {
  num <- as.numeric(as.character(substring(mete[i, 1], 14, 15))) # The 14th and 15th digits represent whether it is a cancer sample, 01~09 are cancer samples
  if (num < 10) {
    mete[i, 2] <- "T"
  } else {
    mete[i, 2] <- "N"
  }
}

names(mete) <- c("id", "group")
mete$group <- as.factor(mete$group)
mete <- subset(mete, mete$group == "T")
rna.data <- rna.data[, which(colnames(rna.data) %in% mete$id)]


USP20.high3 <- as.numeric(as.numeric(rna.data["USP20", ]) > median(as.numeric(rna.data["USP20", ])))
group <- cbind(USP20.high3,1-USP20.high3)
rownames(group) <- colnames(rna.data)
colnames(group) <- c('USP20.high','USP20.low')

fit <- lmFit(rna.data, design = group)
fit <- eBayes(fit)
#top 20 genes where the model fit had an adjusted p-value less than 0.05 and a fold-change of at least 2
topTable(fit,coef=2,number=20,p.value=0.05,lfc=log(2,base=2))
