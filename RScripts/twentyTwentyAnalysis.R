mutations <- read.table('../Data/somaticMutations_incNC.txt', sep='\t', as.is=T, header=T)
samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)

erPatients <- table(samples$erStatus)
erPosMutations <- mutations[mutations$sample%in%samples$sample[which(samples$erStatus=='POS')],]
erNegMutations <- mutations[mutations$sample%in%samples$sample[which(samples$erStatus=='NEG')],]


##################################################################COMPUTE SCORES
computeScores <- function(mutations, oncThreshold, noSamples){
#mutations is a dataframe with columns sample, gene, location, mutationType
#oncThreshold is minimum number of mutations in a codon to be considered recurrently mutated
#noSamples is number of samples in cohort to compute mutation prevalence

genesToProcess <- unique(mutations$gene)
result <- matrix(nrow=length(genesToProcess), ncol=5, dimnames=list(sort(genesToProcess), c('ONC','noONC', 'TSG', 'noTSG', 'percentCoding')))
for (thisGene in genesToProcess){
	theseMuts <- mutations[which(mutations$gene==thisGene),]

	#oncogene score
	countCodons <- table(theseMuts$codon[which(theseMuts$mutationType%in%c('missense SNV', 'inframe indel'))])
	oncContributors <- countCodons[countCodons>=oncThreshold]
	result[thisGene, 'ONC'] <- round(sum(oncContributors)/nrow(theseMuts),2)
	result[thisGene, 'noONC'] <- round(sum(oncContributors),2)

		
	#tumour suppressor gene score
	inactivating <- sum(theseMuts$mutationType%in%c('nonsense SNV', 'frameshift indel', 'stoploss SNV')|theseMuts$location%in%c('splicing', 'exonic;splicing'), na.rm=T)
	result[thisGene, 'TSG'] <- round(inactivating/nrow(theseMuts),2)
	result[thisGene, 'noTSG'] <- inactivating

	
	#percent coding mutations in cohort
	codingMutations <- theseMuts[which(theseMuts$mutationType!='silent SNV'|theseMuts$location%in%c('splicing', 'exonic;splicing')),]
	result[thisGene, 'percentCoding'] <- round(length(unique(codingMutations$sample))/noSamples,2)
	}
	
	return(result)
}
#################################################################################

##################################################################FIND DRIVERS
findDrivers <- function(scores, threshold, oncExclude, noOnc, noTSG){
#result is from computeScores()
#threshold is threshold for calling drivers (0.2 in Vogelstein et al.)
#oncExclude is the threshold required to call a gene a TSG even if the ONC score is above the threshold
#noOnc is minimum number of recurrent mutations required to call oncogene
#noTSG is minimum number of inactivating mutations required to call tumour suppressor gene


oncogenes <- cbind(scores[which(scores[,'ONC']>=threshold & scores[,'TSG']<oncExclude & scores[,'noONC']>=noOnc),], data.frame(type='oncogene'))
tumourSupps <- cbind(scores[which((scores[,'TSG']>=threshold|(scores[,'TSG']>=oncExclude&scores[,'ONC']>=threshold)) & scores[,'noTSG']>=noTSG),], data.frame(type='TSG'))

result <- rbind(oncogenes, tumourSupps)
return(result)
}
#Reference: Vogelstein et al. (2013) Science
#################################################################################

erPosScores <- computeScores(erPosMutations, 3, erPatients['POS'])
erNegScores <- computeScores(erNegMutations, 3, erPatients['NEG'])

posDrivers <- findDrivers(erPosScores, 0.2, 0.05,5,5)
negDrivers <- findDrivers(erNegScores, 0.2, 0.05,5,5)

