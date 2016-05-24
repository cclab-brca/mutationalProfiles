library(GenomicRanges)
library(Hmisc)
samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)

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

source('getFunctionalMutations.R')
driverGRanges <- GRanges(seqnames=driverMuts$chr, ranges=IRanges(start=driverMuts$start, end=driverMuts$end))
driverGRanges$sample <- driverMuts$sample
driverGRanges$vaf <- driverMuts$vaf
driverGRanges$depth <- driverMuts$reads
driverGRanges$gene <- driverMuts$gene
driverGRanges$location <- driverMuts$location
driverGRanges$class <- driverMuts$class
driverGRanges$intClust <- samples$integrativeCluster[match(driverGRanges$sample,samples$sample)]
driverGRanges <- driverGRanges[!is.na(driverGRanges$intClust)]
driverGRanges <- driverGRanges[driverGRanges$sample%in%samplesWcopyNumber]

###Add CN info
driverGRanges$CN <- 2
driverGRanges$purity <- NA

for (i in 1:length(samplesWcopyNumber)){
theseVariants <- which(driverGRanges$sample==samplesWcopyNumber[i])
theseCN <- ascatGRanges[ascatGRanges$sample==samplesWcopyNumber[i]]
ovlp <- findOverlaps(driverGRanges[theseVariants], theseCN)
driverGRanges$purity[theseVariants] <- theseCN$purity[1]

if (length(ovlp)>0){
driverGRanges$CN[theseVariants][queryHits(ovlp)]  <- theseCN$total[subjectHits(ovlp)]

}
}

##Set up results
driverGRanges$ccf <- rep(NA, length(driverGRanges))
driverGRanges$upperLimit <- rep(NA, length(driverGRanges))
driverGRanges$clonalType <- 'CLONAL'

##Compute CCF
for (i in 1:length(driverGRanges)){
confInterval <- binconf(x=driverGRanges$vaf[i]*driverGRanges$depth[i], n=driverGRanges$depth[i])

driverGRanges$ccf[i] <- min(1,round((confInterval[1]/driverGRanges$purity[i]) * ( driverGRanges$purity[i]*driverGRanges$CN[i] + (1-driverGRanges$purity[i])*2    ),2))

driverGRanges$lowerLimit[i] <- min(1,round((confInterval[2]/driverGRanges$purity[i]) * ( driverGRanges$purity[i]*driverGRanges$CN[i] + (1-driverGRanges$purity[i])*2    ),2))

driverGRanges$upperLimit[i] <- min(1,round((confInterval[3]/driverGRanges$purity[i]) * ( driverGRanges$purity[i]*driverGRanges$CN[i] + (1-driverGRanges$purity[i])*2    ),2))
}

##By Integrative Cluster
ics <- c(1,2,3,'4+','4-',5,6,7,8,9,10)
result <- sapply(ics, function(ic) sapply(unique(driverGRanges$gene), function(gene) median(driverGRanges$ccf[which(driverGRanges$intClust==ic&driverGRanges$gene==gene)], na.rm=T)))
