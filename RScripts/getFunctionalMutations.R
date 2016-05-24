##Produces a binary matrix of gene vs sample
##0 - indicates no FUNCTIONAL mutations
##1 - FUNCTIONAL mutation present
##'Non-functional' mutations may be present even if entry is 0
##Functional mutations defined as those that contribute to ONC score for oncogenes and TSG score for tumour suppressors
##Both ONC and TSG mutations considered for TP53, as it has an ONC score >0.2 and a TSG score >0.05

mutations <- read.table('../Data/somaticMutations_incNC.txt', sep='\t', as.is=T, header=T)
samples <- read.table('../Data/patientData.txt', sep='\t', as.is=T, header=T)

###Identify functional mutations
source('twentyTwentyAnalysis.R')
mutations <- mutations[which(mutations$gene%in%c(rownames(posDrivers), rownames(negDrivers))),]


####ONCOGENES
#Use 0.2 as we want tumour suppressors with high ONC scores too (eg.TP53)
oncogenes <- c(rownames(posDrivers[posDrivers$ONC>=0.2,]), rownames(negDrivers[negDrivers$ONC>=0.2,]))
oncThreshold <- 3
inframe_missense <- mutations[which(mutations$mutationType%in%c('inframe indel','missense SNV')&!mutations$location%in%c('splicing', 'exonic;splicing')),]

inframe_missense <- inframe_missense[which(inframe_missense$gene%in%oncogenes),]
recurrentCodons <- table(paste(inframe_missense$gene, inframe_missense$codon, sep='_'))
recurrentCodons <- recurrentCodons[recurrentCodons>=oncThreshold]
oncMutations <- inframe_missense[paste(inframe_missense$gene, inframe_missense$codon, sep='_')%in%names(recurrentCodons),]
oncMutations <- oncMutations[which(oncMutations$gene%in%oncogenes),]

####TUMOUR SUPPRESOR GENES
tsgs <- c(rownames(posDrivers[posDrivers$TSG>=0.2|(posDrivers$TSG>=0.05&posDrivers$ONC>=0.2),]), rownames(negDrivers[negDrivers$TSG>=0.2|(negDrivers$TSG>=0.05&negDrivers$ONC>=0.2),]))
tsgMutations <- mutations[which(mutations$location%in%c('splicing', 'exonic;splicing')|mutations$mutationType%in%c('nonsense SNV', 'frameshift indel')),]
tsgMutations <- tsgMutations[which(tsgMutations$gene%in%tsgs),]


#####MERGE AND MAKE MATRIX
driverMuts <- rbind(oncMutations, tsgMutations)
mutMatrix <- table(driverMuts$gene,driverMuts$sample)
mutMatrix <- mutMatrix>0

wtSamples <- matrix(0, ncol=length(setdiff(samples$sample, colnames(mutMatrix))), nrow=nrow(mutMatrix), dimnames=list(rownames(mutMatrix), setdiff(samples$sample, colnames(mutMatrix)))) #no mutations in driver genes 
mutMatrix <- cbind(mutMatrix, wtSamples)


