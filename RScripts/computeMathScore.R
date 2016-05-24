library(GenomicRanges)
library(Hmisc)
samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)
allMutations <- read.table('../Data/somaticMutations_incNC.txt', sep='\t', as.is=T,  header=T)

ascat <- read.table('../Data/ascatSegments_withoutCNVs.txt', sep='\t', as.is=T, header=T)
ascatGRanges <- GRanges(seqnames=ascat$chr, ranges=IRanges(start=ascat$start, end=ascat$end))
ascatGRanges$sample <- ascat$sample
ascatGRanges$nMajor <- ascat$nMajor
ascatGRanges$nMinor <- ascat$nMinor
ascatGRanges$total <- ascat$nMinor+ascat$nMajor
ascatGRanges$purity <- ascat$purity
ascatGRanges$ploidy <- ascat$ploidy
ascatGRanges <- ascatGRanges[ascatGRanges$purity!=1]
samplesWcopyNumber <- unique(ascatGRanges$sample)

#Correct calls for ploidy & take LOH into account
ascatGRanges$correctedCall <- round((ascatGRanges$total/ascatGRanges$ploidy)*2)
lohIndices <- which(ascatGRanges$correctedCall>=2&(ascatGRanges$nMajor*ascatGRanges$nMinor==0))
#ascatGRanges$correctedCall[lohIndices] <- ascatGRanges$correctedCall[lohIndices]-0.5
ascatGRanges <- ascatGRanges[ascatGRanges$correctedCall!=2]

###Chromosomal instability score
#This is % of genome altered
cin <- numeric(length(samplesWcopyNumber))
mathScore <- numeric(length(samplesWcopyNumber))
numberMutations <- numeric(length(samplesWcopyNumber))
names(cin) <- samplesWcopyNumber
names(mathScore) <- samplesWcopyNumber
names(numberMutations) <- samplesWcopyNumber


for (i in 1:length(samplesWcopyNumber)){
	thisCN <- ascatGRanges[which(ascatGRanges$sample==samplesWcopyNumber[i])]
    cin[samplesWcopyNumber[i]] <- sum(as.numeric(width(thisCN)))
    
    thisMuts <- allMutations[allMutations$sample==samplesWcopyNumber[i],]
    madVaf <- mad(thisMuts$vaf)
    mathScore[samplesWcopyNumber[i]] <- madVaf/median(thisMuts$vaf, na.rm=T)
    numberMutations[samplesWcopyNumber[i]] <- nrow(thisMuts)
    
}

ics <- c(1,2,3,'4+','4-',5,6,7,8,9,10)
cin <- cin/3000000000
cin <- cin[numberMutations>=5]
mathScore <- mathScore[numberMutations>=5]
sampleIcs <- samples$integrativeCluster[match(names(cin), samples$sample)]
cinScoreMedian <- sapply(ics, function(m) median(cin[which(sampleIcs==m)], na.rm=T))
mathScoreMedian <- sapply(ics, function(m) median(mathScore[which(sampleIcs==m)], na.rm=T))
