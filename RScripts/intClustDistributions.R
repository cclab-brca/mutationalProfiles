source('getFunctionalMutations.R')  #We are looking only at 'functional' ,mutations
samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)
intClustLevels <- c('1','2', '3', '4+', '4-', '5', '6', '7', '8', '9', '10')
erPosSamples <- samples$sample[which(samples$erStatus=='POS')]

#Get integrative clusters.
#Uncomment the following line to get to Fig 6A (ER+ only), otherwise all samples will be considered
mutMatrix <- mutMatrix[, colnames(mutMatrix)%in%erPosSamples]
sampleIcs <- samples$integrativeCluster[match(colnames(mutMatrix), samples$sample)]
mutMatrix <- mutMatrix[,!is.na(sampleIcs)] #Remove samples without IntClust assigned
sampleIcs <- sampleIcs[!is.na(sampleIcs)]

##Run over/under representation analysis
pvals <- matrix(nrow=length(intClustLevels),ncol= nrow(mutMatrix), dimnames=list(intClustLevels, rownames(mutMatrix)))
oddsR <- matrix(nrow=length(intClustLevels),ncol=nrow(mutMatrix), dimnames=list(intClustLevels, rownames(mutMatrix)))
mutRate <- matrix(nrow=length(intClustLevels), ncol=nrow(mutMatrix), dimnames=list(intClustLevels, rownames(mutMatrix)))

genes <- rownames(mutMatrix)
for (gene in genes) {
    for (ic in intClustLevels) {
	tmp2 <- factor(mutMatrix[gene,], levels=c(1, 0), labels=c("YES", "NO"))
	tmp1 <- factor(I(1 * (sampleIcs==ic)), levels=c(1, 0), labels=c("IC", "Rest"))
	if (table(tmp2)['YES'] > 0) {
		testMatrix <- table(tmp1, tmp2)
		res <- fisher.test(testMatrix)
		pvals[ic,gene] <- res$p.value
		oddsR[ic,gene] <- res$estimate
		mutRate[ic,gene] <- testMatrix[1,1]/table(sampleIcs)[ic]
  		}	 
	}
}

pcorr <- matrix(p.adjust(pvals, method='fdr'), nrow=length(intClustLevels), dimnames=list(intClustLevels, rownames(mutMatrix)))

