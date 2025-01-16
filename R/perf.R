## PERF.R: Functions for computing various performance measures.
##
## Author: Vlad Popovici
##         Bioinformatics Core Facility
##         Swiss Institute of Bioinformatics
##

## PERF.MATRIX
## Constructs a confusion matrix given the true labels and
## the predicted labels, FOR BINARY CLASSIFICATION PROBLEMS.
##
## The label of the first class ('negative') is lneg, while
## the label of the other class ('positive') is lpos. By
## default, lneg==0 and lpos==1.
##
## Parameters:
##   true.y  - the true labeling of the data
##   pred.y  - the predicted labels
##   lneg, lpos - the labels of the two classes
##
## Value:
##  the confusion matrix
##
perf.matrix = function(true.y, pred.y, lneg=0, lpos=1)
{
	if (length(true.y) != length(pred.y))
		stop("The two label vectors must have the same length!")
	
	r = matrix(0, nrow=2, ncol=2)
	colnames(r) = c(paste("true", lpos, sep='.'), paste("true", lneg, sep='.'))
	rownames(r) = c(paste("pred", lpos, sep='.'), paste("pred", lneg, sep='.'))
	
	r[1,1] = sum(pred.y[true.y == lpos] == lpos)  # correct positive predictions
	r[1,2] = sum(pred.y[true.y == lneg] == lpos)  # incorrect pos. predictions
	r[2,1] = sum(pred.y[true.y == lpos] == lneg)  # incorrect neg. predictions
	r[2,2] = sum(pred.y[true.y == lneg] == lneg)  # correct negative predictions
	
	return (r)
}
## end PERF.MATRIX

## The following 6 functions are computing some performance
## indicators based on the confusion matrix. The names of
## the functions are self-explanatory.
perf.sensitivity = function(R)
{
	return (R[1,1] / (R[1,1]+R[2,1]))
}

perf.specificity = function(R)
{
	return (R[2,2] / (R[1,2]+R[2,2]))
}

perf.error = function(R)
{
	return ((R[1,2]+R[2,1]) / sum(R))
}

perf.ppv = function(R)
{
	return (R[1,1]/(R[1,1]+R[1,2]))
}

perf.npv = function(R)
{
	return (R[2,2]/(R[2,1]+R[2,2]))
}

perf.mcc = function(R) # MCC = Matthew Correlation Coefficient
{
	return ((R[1,1]*R[2,2]-R[1,2]*R[2,1])/sqrt((R[1,1]+R[2,1])*(R[1,1]+R[1,2])*
							(R[2,2]+R[2,1])*(R[2,2]+R[1,2])))
}

## AUC: Area under the curve
##
## It assumes that y.true contains the true class posterior probabilities
## as either 0 or 1 (corresponding to the "negative" and "positive" classes
## and that y.prob contains the predicted probabilities that an object
## belongs to the "positive" class (so that all objects having a probability
## higher than a threshold will be labeled as "positive").
perf.auc = function(y.prob, y.true, lneg=0, lpos=1)
{
	i = is.na(y.prob)
	y.true = y.true[!i]
	y.prob = y.prob[!i]
	
	n0  = sum(y.true == lneg)
	n1  = sum(y.true == lpos)
	r   = rank(y.prob, ties.method="random")
	s   = sum(r[y.true==lpos])
	auc = (s - 0.5*n1*(n1+1)) / (n0*n1)
	
	return (auc)
}

## RMSE: Root mean square error.
##
perf.rmse = function(y.prob, y.true)
{
	i = is.na(y.prob)
	y.true = y.true[!i]
	y.prob = y.prob[!i]
	
	rmse = sqrt(mean((y.true - y.prob)^2, na.rm=TRUE))
	
	return (rmse)
}

