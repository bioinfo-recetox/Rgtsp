## KCV.R: generates the folds for k-fold cross-validation.
##
## Author: Vlad Popovici
##         Bioinformatics Core Facility
##         Swiss Institute of Bioinformatics
##

## KCV.SPLITS
## Generates a list of index sets, such that each set
## corresponds to a fold (split) in the cross-validation.
##
## Parameters:
##   labels - a vector of labels (usually {0,1} or {-1,1})
##   k      - number of folds in validation
##   stratified - should the sampling be stratified? (i.e.
##            preserving the relative proportion of the
##            two classes)
##
## Values:
##   a list with 2 components: "cv.train" and "cv.valid"
##   $cv.train is a list of length k containing k lists
##     of indexes for training sets
##   $cv.valid is a list of length k containing k lists
##     of indexes for validation sets
##
## Example:
##   cv = kcv.splits(c(1,2,1,2,1,2,1,2,1,2), k=5, stratified=TRUE)
## Then
##   cv$cv.train
## contains the training sets, and
##   cv$cv.valid
## contains the validation sets:
##
## cv$cv.train[[1]] is the first training set and
## cv$cv.valid[[1]] is the first validation set
###############################################################################

kcv.splits = function(labels, k=10, stratified=FALSE) 
{
	n  = length(labels)
	nc = length(unique(labels))
	
	if (nc != 2) {
		stop('Currently, only 2-class settings are supported')
	}
	lneg = unique(labels)[1]  # first class label ("negative")
	lpos = unique(labels)[2]  # second class label ("positive")
	
	if (stratified) {
		# Initially, we consider the indexes of the positive and
		# negative class as being 0..np-1 and 0..nn-1, respectively.
		# After partitioning in folds, we map these indexes
		# onto original ones.
		np = sum(labels == lpos)    # count the "positive" class
		nn = sum(labels == lneg)    # count the "negative" class
		fp = floor(np / k)          # size of a positive fold
		fn = floor(nn / k)          # size of a negative fold
		rp = np - k*fp              # positive leftovers
		rn = nn - k*fn              # negative leftovers
		
		p.idx = (1:n)[labels == lpos]; p.idx = sample(p.idx)[1:(np-rp)]
		n.idx = (1:n)[labels == lneg]; n.idx = sample(n.idx)[1:(nn-rn)]
		train.idx = list()
		valid.idx = list()
		for (j in 1:k) {
			# find the indexes for the validation set
			vd = c(p.idx[((j-1)*fp+1):(j*fp)], n.idx[((j-1)*fn+1):(j*fn)])
			vd = sample(vd)
			
			# the training indexes are all but the validation ones:
			tr = setdiff(c(p.idx,n.idx), vd)
			tr = sample(tr)
			train.idx = c(train.idx, list(tr))
			valid.idx = c(valid.idx, list(vd))
		}
	}
	else {
		# Don't bother with stratification: just shuffle data
		# and produce folds.
		f   = floor(n / k)
		
		train.idx = list()
		valid.idx = list()
		folds = split(sample(seq(n)), rep(1:k, length=n))
		
		for (j in 1:k) {
			vd = folds[[k]]
			tr = setdiff(1:n, vd)
			tr = sample(tr)
			train.idx = c(train.idx, list(tr))
			valid.idx = c(valid.idx, list(vd))      
		}
	}
	
	return (list(cv.train=train.idx, cv.valid=valid.idx))
}

## end kcv.splits

