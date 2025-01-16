% File man/tsp.n.Rd
\name{tsp.n}
\alias{tsp.n}
\title{Top Scoring N Pairs}
\description{
    Return the n highest scoring pairs of features, of the form (i,j) such that
    X[,i] < X[,j] for class "+1", and X[,i] > X[,j] for class "0".
}
\usage{
tsp.n(X, y, n, score="geman", weights=NULL, remap.labels=TRUE)
}
\arguments{
    \item{X}{a data matrix, with features given by columns and samples by rows.}
    \item{y}{a vector of labels (preferably with 0/1 coding).}
    \item{n}{the maximum number of pairs to return.}
    \item{score}{the scoring method. If score=="geman", then the classical
    version of the TSP is used. If score=="lreg" then the logistic regression
    is used for scoring the pairs.}
    \item{weights}{(optional) a vector containing sample weights.}
    \item{remap.labels}{(optional) if labels are not 0/1 coded, then a remapping
    is needed.}
}
\value{
    A list with three components:
    \item{I}{the vector of first feature indexes}
    \item{J}{the vector of second feature indexes}
    \item{S}{the vector of scores}
    \item{labels.map}{a hash table with the remapped labels}
    
    A feature pair has the form (I[k], J[k]) with the corresponding score S[k].
}
\seealso{
    \code{\link{tsp.s}}
}
\examples{
library(Rgtsp)
data(mdabr)
z = tsp.n(X, y.erpos, 10)   # X is a data matrix and y a labels vector
data.frame(I=colnames(X)[z$I], J=colnames(X)[z$J])
}
\references{
    Geman, D. and d'Avignon, C. and Naiman, D.Q. and
    Winslow, R.L. (2004) Classifying Gene Expression
    Profiles from Pairwise mRNA Comparisons,
    \emph{Stat Appl Genet Mol Biol}, \bold{3}, Article 19
}
\author{Vlad Popovici, \email{vlad.popovici@isb-sib.ch}}
