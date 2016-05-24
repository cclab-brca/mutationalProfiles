library(survival)
library(multcomp)

samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)
samples$daysFollowup <- as.Date(samples$dateLastFollowup, format="%Y-%m-%d")-as.Date(samples$dateOfDiagnosis, format="%Y-%m-%d")
posSamples <- samples[which(samples$erStatus=='POS'),]

source('getFunctionalMutations.R')
driverMuts$intClust <- samples$integrativeCluster[match(driverMuts$sample, samples$sample)]
pik3ca <- driverMuts[which(driverMuts$gene=='PIK3CA'),]
pik3ca <- pik3ca[pik3ca$sample%in%posSamples$sample,]

intClusts <- c(1,2,3,'4+', 5,6,7,8,9,10)
interactionPvals <- numeric(10); names(interactionPvals) <- intClusts
results <- list()

for (ic in intClusts){
	posSamples$geneMut <- posSamples$sample%in%pik3ca$sample
    posSamples$thisIc <- posSamples$integrativeCluster==ic
    survivalObject <- Surv((as.numeric(posSamples$daysFollowup)/30), event=(posSamples$patientStatus=='DDS'))
    model1 <- coxph(survivalObject ~ geneMut + thisIc + geneMut*thisIc, data=posSamples)
	interactionPvals[ic] <- coefficients(summary(model1))['geneMutTRUE:thisIcTRUE',5]

	#To obtain simultaneous confidence intervals
	X <- matrix(c(1,0,0,
    0,1,0,
    1,1,1), nrow=3, byrow=T, dimnames=list(1:3, c('geneMut', 'thisIc','geneMut:thisIc')))
    results[[ic]] <- confint(summary(glht(model1, linfct=X)))
		
}

