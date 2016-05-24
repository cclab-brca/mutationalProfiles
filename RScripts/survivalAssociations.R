library(survival)
samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)
samples$daysFollowup <- as.Date(samples$dateLastFollowup, format="%Y-%m-%d")-as.Date(samples$dateOfDiagnosis, format="%Y-%m-%d")
posSamples <- samples[which(samples$erStatus=='POS'),]
negSamples <- samples[which(samples$erStatus=='NEG'),]

source('getFunctionalMutations.R')
genes <- unique(driverMuts$gene)

##Set up structures to hold results
erPosHR <- rep(NA, length(genes)); names(erPosHR) <- genes
erPosPval <- rep(NA, length(genes)); names(erPosPval) <- genes
erPosLower <- rep(NA, length(genes)); names(erPosLower) <- genes
erPosUpper <- rep(NA, length(genes)); names(erPosUpper) <- genes
erNegHR <- rep(NA, length(genes)); names(erNegHR) <- genes
erNegPval <- rep(NA, length(genes)); names(erNegPval) <- genes
erNegLower <- rep(NA, length(genes)); names(erNegLower) <- genes
erNegUpper <- rep(NA, length(genes)); names(erNegUpper) <- genes


for (gene in genes){
thisGene <- driverMuts[which(driverMuts$gene==gene),]
posSamples$geneMut <- posSamples$sample%in%thisGene$sample
negSamples$geneMut <- negSamples$sample%in%thisGene$sample

posSurv <- Surv((as.numeric(posSamples$daysFollowup)/30), event=(posSamples$patientStatus=='DDS'))
negSurv <- Surv((as.numeric(negSamples$daysFollowup)/30), event=(negSamples$patientStatus=='DDS'))

posModel <- coxph(posSurv ~ geneMut + grade + ageGreater55 + sizeGreater50 + lnPositive, data=posSamples)
negModel <- coxph(negSurv ~ geneMut + grade + ageGreater55 + sizeGreater50 + lnPositive, data=negSamples)

posResult <- summary(posModel)
negResult <- summary(negModel)

erPosHR[gene] <- posResult$conf.int['geneMutTRUE',1]
erPosLower[gene] <- posResult$conf.int['geneMutTRUE',3]
erPosUpper[gene] <- posResult$conf.int['geneMutTRUE',4]
erPosPval[gene] <- coefficients(posResult)['geneMutTRUE',5]

erNegHR[gene] <- negResult$conf.int['geneMutTRUE',1]
erNegLower[gene] <- negResult$conf.int['geneMutTRUE',3]
erNegUpper[gene] <- negResult$conf.int['geneMutTRUE',4]
erNegPval[gene] <- coefficients(negResult)['geneMutTRUE',5]
}

#Samples with minimum info
numberPosSamples <- posSamples$n
numberNegSamples <- negSamples$n
