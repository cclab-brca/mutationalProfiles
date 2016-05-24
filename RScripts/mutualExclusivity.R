source('../Data/getFunctionalMutations.R')

testMatrix <- mutMatrix
testMatrix <- testMatrix[rowSums(testMatrix)>=0.005*ncol(testMatrix),]
features <- rownames(testMatrix)

pvals <- matrix(nrow=length(features), ncol=length(features), dimnames=list(features, features))
oddsRatio <- matrix(nrow=length(features), ncol=length(features), dimnames=list(features, features))
confint <- matrix(nrow=length(features), ncol=length(features), dimnames=list(features, features))

for (i in 1:length(features)){
	for (j in 1:length(features)){
	if (i !=j){

	toTest <- table(testMatrix[i,], testMatrix[j,])
	if (ncol(toTest)>1&nrow(toTest)>1){
	fisherRes <- fisher.test(toTest)
	pvals[i,j] <- fisherRes$p.value
	oddsRatio[i,j] <- fisherRes$estimate
	confint[i,j] <- paste(signif(fisherRes$conf.int,2), collapse='-')
}
	}}
}


###FDR Correction
p.upper <- pvals
p.upper[lower.tri(p.upper)] <- NA
p.upper.corr <- matrix(p.adjust(p.upper, method='fdr'), ncol=ncol(p.upper), nrow=nrow(p.upper), dimnames=dimnames(p.upper))

###LOG ODDS
oddsRatio[is.na(oddsRatio)] <-1
logOR <- log(oddsRatio+0.01)