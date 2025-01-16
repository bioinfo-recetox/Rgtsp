% File man/cv.tsp.Rd
\name{cv.tsp}
\alias{cv.tsp}
\title{K-fold Cross Validation for TSP models}
\description{
    Performs a k-fold (stratified) cross-validation for TSP models.
}
\usage{
cv.tsp(X, y, kfold=5, which.tsp="tsp.n", max.pairs=10, min.score=0.7,
				combiner="majority", weighted=FALSE)
}
\arguments{
    \item{X}{a data matrix, with features given by columns and samples by rows.}
    \item{y}{a vector of labels: 0,1}
    \item{kfold}{(optional) number of folds to split the data set into.}
    \item{which.tsp}{(optional) which version of the TSP to use: "tsp.n" or
    "tsp.s".}
    \item{max.pairs}{(optional) maximum number of pairs - parameter for
    tsp.n() and tsp.s().}
    \item{min.score}{(optional) minimum score for TSPs - parameter for
    tsp.s().}
    \item{combiner}{(optional) how to combiner the TSPs for predicting the
    labels ("majority", "wmajority") - parameter for predict.tsp().}
    \item{weigthed}{(optional) should the pairs' individual predictions be
    weighted? - parameter for predict.tsp()}
}
\value{
    A list with four components (vectors with values for error rate,
		sensitivity, specificity and AUC):
		\item{$tr.m}{the mean values for training performance}
		\item{$tr.s}{standard deviation for training performance}
		\item{$vd.m}{the mean values for validation performance}
		\item{$vd.s}{standard deviation for validation performance}
}
\seealso{
    \code{\link{tsp.n}}, \code{\link{tsp.s}}, \code{\link{predict.tspair}}
}
\author{Vlad Popovici, \email{popovici@recetox.muni.cz}}
