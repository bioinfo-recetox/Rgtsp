% File man/predict.tspair.Rd
\name{predict.tspair}
\alias{predict.tspiar}
\title{Predict the labels with Top Scoring Pairs}
\description{
    Once a model is built, there are several ways of predicting the labels. This
    function predicts the labels using either simple TSPs or TSP hubs.
}
\usage{
predict.tspair(tsp.obj, X, combiner="none", weighted=FALSE,
    min.score=min(tsp.obj$S))
}
\arguments{
    \item{tsp.obj}{a TSP object, as returned by tsp.n()/tsp.s() or tsp.hub()
    functions}
    \item{X}{data matrix whose labels are to be predicted (samples by rows,
    features by columns).}
    \item{combiner}{(optional) how the individual preditions of the TSPs or
    hubs are to be combined.
    If combiner == "none", the individual predictions are not combined, and
    are returned as such.
    If combiner == "majority", the individual predictions are combined through
    majority vote.
    If combiner == "wmajority", the individual predictions are combined through
    weighted majority vote.}
    \item{weigthed}{if true, the individual decisions could be biased by the
    use of weights. This is currently implemented only for TSPs, and the
    weights are taken from the $S component of the tsp.obj object.}
    \item{min.score}{for predicting the labels, use only those pairs with a
    score at least <min.score>.}
}
\value{
    Either a vector (one value for each row in X) or a matrix (one column for
    each TSP or hub, one row per sample) with 0/1 labels. In the case the returned
    value is a matrix, the individual predictions were not combined.
}
\seealso{
    \code{\link{tsp.s}}, \code{\link{tsp.n}}, \code{\link{tsp.hub}}
}
%\examples{
%}
\author{Vlad Popovici, \email{vlad.popovici@isb-sib.ch}}
