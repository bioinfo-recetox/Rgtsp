################################################################################
## Rgtsp library                                                              ##
##                                                                            ##
## Author: Vlad Popovici, vlad.popovici@isb.sib.ch                            ##
################################################################################


### Utility functions for converting labels to numerical labels (forth and
### back).
build.labels.map = function(y)
{
    original.class = class(y)
    y = as.character(y)
    
    labels.map = 0 : (length(unique(y))-1)
    names(labels.map) = sort(unique(y))
    
    return (list(labels.map=labels.map, original.class=original.class))
}

labels.to.numeric = function(y, map)
{
    y = as.character(y)
    yn  = rep(map$labels.map[1], length(y))
    names(yn) = names(y)
    
    for (l in map$labels.map) { # l = 0, 1, ...
        yn[y == names(map$labels.map)[l+1]] = l
    }
    
    return (yn)
}

numeric.to.labels = function(yn, map)
{
    if (!is.null(map)) {
        y0 = rep("", length(yn))
        names(y0) = names(yn)
        
        for (l in map$labels.map) {
            y0[yn == l] = names(map$labels.map)[l+1]
        }
        
        if (map$original.class == "factor") {
            y = as.factor(y0)
        }
        else {
            y = as(y0, map$original.class)
        }
        names(y) = names(yn)
    }
    else
        y = yn
    
    return (y)
}

##
## TSP.N: Return the top N scoring pairs.
##
## Parameters:
##  X       - data matrix: features by columns, samples by rows
##  y       - data labels: 0/1
##  n       - number of top pairs to return
##  score   - scoring method to be used: only one value, for the moment, "geman"
##  weights - if provided, allows weighting the samples
##
## Returns:
##
tsp.n = function(X, y, n,
  score="geman",
  weights=NULL,
  remap.labels=TRUE)
{
    if (nrow(X) != length(y))
        stop('Number of rows in X must equal the number of elements in y!')
  
    if (sum(rownames(X) == names(y)) != length(y))
        stop('X row names must match the names of y!')
  
    if (!is.null(weights) && (length(weights) != length(y)))
        stop('Weights vector length must equal the number of samples!')
      
    if (!(score %in% c("geman")))
        stop('Unknown score specified!')
      
    if (any(is.na(X)) || any(is.na(y)))
        stop("NAs detected in you data. Please correct.")
      
    if (!is.null(weights) && any(is.na(weights)))
        stop('NAs detected in weights vector! Please correct.')
      
    if (remap.labels) {
        lm = build.labels.map(y)
        y = labels.to.numeric(y, lm)         # encode the labels to be 0/1

        if (length(lm$labels.map) < 2) {
            stop('Exactly 2 classes must be given.')
        }
        if (length(lm$labels.map) > 2) {
            stop('Use mtsp() for multiclass problems!')
        }
    }
    else lm = NULL
    
    n = min(max(n,1), 0.5*ncol(X)*(ncol(X)-1))
      
    i = rep(0, n)
    j = rep(0, n)
    s = rep(0, n)
    z = list(I=i, J=j, S=s)
      
    if (is.null(weights)) {
        if (score == "geman") {
            z = .C("RGTSP_tsp_N",
                as.double(X),             # X is stored in column-major order (so C gets t(X))
                as.integer(y),
                as.integer(ncol(X)),      # number of rows in C <-> number of cols in R
                as.integer(nrow(X)),      # number of cols in C <-> number of rows in R
                i=as.integer(i),
                j=as.integer(j),
                s=as.double(s),
                n=as.integer(n))
        }
        #else { # others to come
        #}
    }
    else {
        if (score == "geman") {
            z = .C("RGTSP_tsp_NW",
                as.double(X),             # X is stored in column-major order (so C gets t(X))
                as.integer(y),
                as.double(weights),
                as.integer(ncol(X)),      # number of rows in C <-> number of cols in R
                as.integer(nrow(X)),      # number of cols in C <-> number of rows in R
                i=as.integer(i),
                j=as.integer(j),
                s=as.double(s),
                n=as.integer(n))
        }
        #else { # others to come
        #}
  }
  
  TSP = list(I=z$i[1:z$n], J=z$j[1:z$n], S=z$s[1:z$n], labels.map=lm)
  class(TSP) = "tspair"
  return (TSP)
}

##
## TSP.S: Return the top scoring pairs with a score larger than a given minimum.
##
## Parameters:
##  X       - data matrix: features by columns, samples by rows
##  y       - data labels: 0/1
##  ms      - minimum score required
##  n       - maximum number of pairs to return (all having a score >= ms)
##  score   - scoring method to be used: "geman", "lreg"
##  weights - if provided, allows weighting the samples
##
## Returns:
##
tsp.s = function(X, y, ms, n, score="geman", weights=NULL, remap.labels=TRUE)
{
    if (nrow(X) != length(y))
      stop('Number of rows in X must equal the number of elements in y!')
    
    if (sum(rownames(X) == names(y)) != length(y))
      stop('X row names must match the names of y!')
    
    if (!is.null(weights) && (length(weights) != length(y)))
      stop('Weights vector length must equal the number of samples!')
    
    if (!(score %in% c("geman", "lreg")))
      stop('Unknown score specified!')
    
    if (any(is.na(X)) || any(is.na(y)))
      stop("NAs detected in you data. Please correct.")
    
    if (!is.null(weights) && any(is.na(weights)))
      stop('NAs detected in weights vector! Please correct.')
  
    n = min(max(n,1), 0.5*ncol(X)*(ncol(X)-1))
    
    ms = min(max(0, ms), 1)       # the score must be between 0 and 1

    if (remap.labels) {
        lm = build.labels.map(y)
        y = labels.to.numeric(y, lm)         # encode the labels to be 0/1

        if (length(lm$labels.map) < 2) {
            stop('Exactly 2 classes must be given.')
        }
        if (length(lm$labels.map) > 2) {
            stop('Use mtsp() for multiclass problems!')
        }
    }
    else lm = NULL

    i = rep(0, n)
    j = rep(0, n)
    s = rep(0, n)
    z = list(I=i, J=j, S=s)
    
    if (is.null(weights)) {
        if (score == "geman") {
            z = .C("RGTSP_tsp_S",
                as.double(X),             # X is stored in column-major order (so C gets t(X))
                as.integer(y),
                as.integer(ncol(X)),      # number of rows in C <-> number of cols in R
                as.integer(nrow(X)),      # number of cols in C <-> number of rows in R
                as.double(ms),
                i=as.integer(i),
                j=as.integer(j),
                s=as.double(s),
                n=as.integer(n))
        }
        #else { # others to come
        #}
    }
    else {
        if (score == "geman") {
            z = .C("RGTSP_tsp_SW",
                as.double(X),             # X is stored in column-major order (so C gets t(X))
                as.integer(y),
                as.double(weights),
                as.integer(ncol(X)),      # number of rows in C <-> number of cols in R
                as.integer(nrow(X)),      # number of cols in C <-> number of rows in R
                as.double(ms),
                i=as.integer(i),
                j=as.integer(j),
                s=as.double(s),
                n=as.integer(n))
        }
        #else { # others to come
        #}
    }
  
    TSP = list(I=z$i[1:z$n], J=z$j[1:z$n], S=z$s[1:z$n], labels.map=lm)
    class(TSP) = "tspair"
    
    return (TSP)
}


##
## TSP.HUB: try to construct a number of hubs from a list of TSPs
##
## Parameters:
##
##   tsp - a tspair object, with pairs of indexes, from which the hubs
##         are to be constructed
##   min.hub.size - minimum number of pairs in a hub
##
tsp.hub = function(tsp, min.hub.size=5)
{
    if (attr(tsp, 'class') != 'tspair')
        stop('First argument must be a tspair object!')
        
    # allocate memory:
    m = length(tsp$I)
    hc = rep(0, m)
    hn = rep(0, m)
    hs = rep(0, m)
    
    z = .C("RGTSP_tsp_hub",
           as.integer(tsp$I),
           as.integer(tsp$J),
           as.double(tsp$S),
           as.integer(length(tsp$I)),
           hc = as.integer(hc),
           hn = as.integer(hn),
           hs = as.integer(hs),
           m = as.integer(m),
           as.integer(min.hub.size)
       )
    
    hubs = NULL

    cnt = setdiff(unique(z$hc), c(0))  # unique hub centers, without 0
    
    for (k in cnt) {
        s = z$hs[z$hc==k][1] # get the sign
        h = list(center=k, sign=s, neighbors=z$hn[z$hc == k])
        hubs = c(hubs, h)
    }
    hubs = c(hubs, labels.map=tsp$labels.map)
    class(hubs) = 'tsphub'
    
    return (hubs)
}


##
## IDX2NAMES: instead of variable indexes, use variable names.
##
## Parameters:
##   tsp.obj - a tspair object
##   var.names - a list of variable names (usually colnames(X), where X
##            was used to build the TSP object)
##
idx2names = function(tsp.obj, var.names)
{
    tsp.obj$I = var.names[tsp.obj$I]
    tsp.obj$J = var.names[tsp.obj$J]
    
    return (tsp.obj)
}


##
## PREDICT.TSPAIR: predict the labels of a data set
##
## Parameters:
##   tsp.obj - a TSP object, either tspair or tsphub
##   X  - data matrix with samples by rows, features by columns
##   combiner - "none", "majority", "wmajority"
##   weghted - [boolean] use weighted prediction?
##   min.score - use only the pairs having at least <min.score> score
##
predict.tspair = function(tsp.obj, X,
                       combiner="none",
                       weighted=FALSE,
                       min.score=min(tsp.obj$S),
                       ...)
{
    weights = NULL
    
    if (!(combiner %in% c('none', 'majority', 'wmajority', 'average')))
        stop(' combiner must be one of "none", "majority", "wmajority" or "average"')
        
    p.maj = function(y) {
        z = as.numeric(rowSums(y) > ncol(y)/2)
        return (z)
    }
    p.wmaj = function(y, w) {
        z = as.numeric(rowSums(y*matrix(w, nrow=nrow(y),ncol=ncol(y),byrow=TRUE)) > sum(w)/2)
        return (z)
    }
    p.avg = function(tsp.list, X) {
        z = rep(0, nrow(X))
        z[rowSums(X[,tsp.list$J]) > rowSums(X[,tsp.list$I])] = 1
        return (z)
    }
    
    if (attr(tsp.obj, 'class') == 'tspair') {
        np = length(tsp.obj$I)
        if (np == 0) stop('No TSPs in tspair object!')
        if (tsp.obj$I[1] == 0 || tsp.obj$J[1] == 0) stop('No TSPs in tspair object!')
   
        if (class(tsp.obj$I) == 'character') {   # use column names, not indexes
            # check for those genes that are present:
            k = intersect(which(tsp.obj$I %in% colnames(X)), 
                which(tsp.obj$J %in% colnames(X)))
            idx = intersect(k, which(tsp.obj$S >= min.score))
        }
        else {
            idx = which(tsp.obj$S >= min.score)
        }
        tsp.obj$I = tsp.obj$I[idx]
        tsp.obj$J = tsp.obj$J[idx]
        tsp.obj$S = tsp.obj$S[idx]
        np = length(tsp.obj$I)
        
        y = matrix(0, ncol=np, nrow=nrow(X))
        rownames(y) = rownames(X)
        colnames(y) = paste('P', colnames(X)[tsp.obj$I], colnames(X)[tsp.obj$J], sep='.')
        weights = tsp.obj$S
        
        if (weighted) {
            for (k in 1:length(tsp.obj$I)) {
                y[,k] = as.numeric(X[,tsp.obj$I[k]] < tsp.obj$S[k]*X[,tsp.obj$J[k]])
            }
        }
        else {
            for (k in 1:length(tsp.obj$I)) {
                y[,k] = as.numeric(X[,tsp.obj$I[k]] < X[,tsp.obj$J[k]])
            }
        }
    }
    else stop('Unknown object!')
    
    y = switch (combiner,
                "none" = y,
                "majority" = p.maj(y),
                "wmajority" = p.wmaj(y, weights),
                "average" = p.avg(tsp.obj, X)
               )
    
    if (combiner == "none") {
        y0 = NULL
        for (k in 1:ncol(y)) {
            y0 = cbind(y0, numeric.to.labels(y[,k], tsp.obj$labels.map))
        }
        colnames(y0) = colnames(y)
        rownames(y0) = rownames(y)
    }
    else {
        if (!is.null(tsp.obj$labels.map)) {
            y0 = numeric.to.labels(y, tsp.obj$labels.map)
        }
        else {
            y0 = y
        }
        names(y0) = rownames(X)
    }
    
    return (y0)
}

##
## PREDICT.TSHUB: predict the labels of a data set
##
## Parameters:
##   tsp.obj - a TSP object, either tspair or tsphub
##   X  - data matrix with samples by rows, features by columns
##   combiner - "none", "majority", "wmajority"
##   weghted - [boolean] use weighted prediction?
##
predict.tsphub = function(tsp.obj, X,
                       combiner="none",
                       weighted=FALSE,...)
{
    weights = NULL
    
    p.maj = function(y) {
        z = as.numeric(rowSums(y) > ncol(y)/2)
        return (z)
    }
    p.wmaj = function(y,w) {
        z = as.numeric(rowSums(y*matrix(w, nrow=nrow(y),ncol=ncol(y),byrow=TRUE)) > sum(w)/2)
        return (z)
    }
    
    if (attr(tsp.obj, 'class') == 'tsphub') {
        # each of the hubs makes its own predictionas that may be combined at
        # the end
        nh = length(tsp.obj) / 3   # each TSP hub takes 3 components
        y = matrix(0, ncol=nh, nrow=nrow(X))
        rownames(y) =  rownames(X)
        colnames(y) = paste('hub', 1:nh, sep='.')
        weights = rep(1, nh)       # not using weights for multiple hubs
        for (h in 0:(nh-1)) {
            i = tsp.obj[[3*h+1]][1]
            s = tsp.obj[[3*h+2]][1]
            n = length(tsp.obj[[3*(h+1)]])
            ytmp = matrix(0, nrow=nrow(X), ncol=n)
            for (k in 1:n) {
                j = tsp.obj[[3*(h+1)]][k]               
                ytmp[,k] = as.numeric(s*X[,i] > s*X[,j])
            }
            y[,h+1] = p.maj(ytmp)
        }
    }
    else stop('Unknown object!')
    
    y = switch (combiner,
                "none" = y,
                "majority" = p.maj(y),
                "wmajority" = p.wmaj(y, weights)
               )
    
    if (!is.null(tsp.obj$labels.map)) {
        y0 = numeric.to.labels(y, tsp.obj$labels.map)
    }
    else {
        y0 = y
    }
    
    names(y0) = rownames(X)
    
    return (y0)
}


print.tspair = function(tsp.obj,...)
{
    if (attr(tsp.obj, 'class') == 'tspair') {
        cat('\t\tIndex var 1\t\tIndex var 2\n')
        for (k in 1:length(tsp.obj$I)) {
            cat(paste('Pair:', k, tsp.obj$I[k], '<', tsp.obj$J[k], 'score = ', round(tsp.obj$S[k],4), '\n', sep='\t'))
        }
    }
    else stop('Unknown object!')
}

print.tsphub = function(tsp.obj,...)
{
    if (attr(tsp.obj, 'class') == 'tsphub') {
        nh = length(tsp.obj) / 3
        for (k in 0:(nh-1)) {
            i = tsp.obj[[3*k+1]][1]
            s = tsp.obj[[3*k+2]][1]
            n = length(tsp.obj[[3*(k+1)]])
            cat(paste('Hub ', k+1, ': ', n, ' pairs\n', sep=''))
            if (s == 1) ss = '>'
            else ss = '<'
            cat(paste('\tCenter:', i, ss, '\n'))
            cat(paste(tsp.obj[[3*(k+1)]], collapse=' '))
            cat('\n')
        }
    }
    else stop('Unknown object!')
}


####
#### Multi-class extension.
####

## Class labels are supposed to be numeric and of the form 0,1,2,3,... so the
## number of classes is given by max_class_label+1.
mtsp = function(X, y, min.score=0.75, max.pairs=10, remap.labels=TRUE, ...)
{
    if (remap.labels) {
        lm = build.labels.map(y)
        yn = labels.to.numeric(y, lm)         # encode the labels to be 0/1

        if (length(lm$labels.map) <= 2) {
            stop('mtsp() is to be used on multiclass problems!')
        }
    }
    else lm = NULL
    
    nc = max(yn) + 1                      # number of classes
    Z = NULL
    all.tsp = list()
    
    # pairwise TSPs:
    for (y1 in 1:(nc-1)) {
        for (y2 in 0:(y1-1)) {
            #cat(paste('Class', y1, 'vs Class', y2, '\n'))
            
            # build TSPs to discriminate y1 from y2:
            i = which(yn == y1 | yn == y2)
            y.tmp = yn[i]
            y.new = rep(-1, length(y.tmp))
            y.new[y.tmp == y1] = 0
            y.new[y.tmp == y2] = 1
            names(y.new) = names(y)[i]
            X.new = X[i,]
            m = tsp.s(X.new, y.new, min.score, max.pairs, remap.labels=FALSE)
            all.tsp$I = c(all.tsp$I, m$I)
            all.tsp$J = c(all.tsp$J, m$J)
            all.tsp$S = c(all.tsp$S, m$S)
        }
    }
    
    all.tsp$labels.map = NULL
    attr(all.tsp, 'class') = 'tspair'        # make sure it's a tspair object
    Z = predict(all.tsp, X)
    #print(colnames(Z))
    i = which(duplicated(colnames(Z)))       # get rid of duplicated TSPs...
    if (length(i) > 0) {
        all.tsp$I = all.tsp$I[-i]
        all.tsp$J = all.tsp$J[-i]
        all.tsp$S = all.tsp$S[-i]
        Z = Z[,-i]
    }
    
    # build the classification tree:
    tr = ctree(yn ~ ., data=data.frame(yn=as.factor(yn), Z), ...)
    
    m = list(TSP=all.tsp, TR=tr, labels.map=lm)
    attr(m, 'class') = 'mtsp'
    
	return (m)
}


predict.mtsp = function(obj, X, ...)
{
    if (attr(obj, 'class') == 'mtsp') {
        Z = predict(obj$TSP, X, combiner="none")
        y = predict(obj$TR, as.data.frame(Z))
    }
    else stop('Unknown object!')
    
    y = as.numeric(as.character(y))
    
    if (is.null(obj$labels.map)) {
        y0 = y
    }
    else {
        y0 = numeric.to.labels(y, obj$labels.map)
    }
    
    names(y0) = rownames(X)
    
    return (y0)
}


## CV.TSP - use k-fold cross-validation to assess the model.
## The actual computation of performance metrics is done elsewhere.
cv.tsp = function(X, y, kfold=5, which.tsp="tsp.n", max.pairs=10, min.score=0.7,
									combiner="majority", weighted=FALSE)
{  
    d = nrow(X)          # total number of features
    lneg = as.character(levels(as.factor(y))[1])   # negative label
    lpos = as.character(levels(as.factor(y))[2])   # positive label
    which.tsp = tolower(which.tsp)
    
    if (!(which.tsp %in% c("tsp.n", "tsp.s"))) {
        stop('Unknown TSP model!')
    }
    r = matrix(0, 2, 2)                 # confusion matrix
    p.tr = matrix(0, kfold, 4); colnames(p.tr) = c("Error.rate", "Sensitivity", "Specificity", "AUC")
    p.vd = matrix(0, kfold, 4); colnames(p.vd) = c("Error.rate", "Sensitivity", "Specificity", "AUC")
    
    # Generate the folds:
		cv = kcv.splits(y, k=kfold, stratified=TRUE)
		
		for (k in 1:kfold) {
        itr = cv$cv.train[[k]]          # index of elements for training
        ivd = cv$cv.valid[[k]]          # index of elements for validation
        Xtr = X[itr,]; ytr = y[itr]
        Xvd = X[ivd,]; yvd = y[ivd]
			
        cat(paste("    Fold:", k, "\n"))
        
        m = switch (which.tsp,          
            "tsp.n" = tsp.n(Xtr, ytr, max.pairs),
            "tsp.s" = tsp.s(Xtr, ytr, min.score, max.pairs)
        )
        
				yp.tr = predict(m, Xtr, combiner=combiner, weighted=weighted)
				yp.vd = predict(m, Xvd, combiner=combiner, weighted=weighted)
        
        r = perf.matrix(ytr, yp.tr, lneg=lneg, lpos=lpos)
        p.tr[k,] = c(perf.error(r), perf.sensitivity(r), perf.specificity(r), 0)
        p.tr[k,4] = 0.5*(p.tr[k,2] + p.tr[k,3])  # AUC is only approximated
        
        r = perf.matrix(yvd, yp.vd, lneg=lneg, lpos=lpos)
        p.vd[k,] = c(perf.error(r), perf.sensitivity(r), perf.specificity(r), 0)
        p.vd[k,4] = 0.5*(p.vd[k,2] + p.vd[k,3])  # AUC is only approximated
		}

    return (list(tr.m=apply(p.tr, 2, mean, na.rm=TRUE), tr.s=apply(p.tr, 2, sd, na.rm=TRUE),
                 vd.m=apply(p.vd, 2, mean, na.rm=TRUE), vd.s=apply(p.vd, 2, sd, na.rm=TRUE)))
}

