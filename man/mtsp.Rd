% File man/mtsp.Rd
\name{mtsp}
\alias{mtsp}
\title{Multiclass Top Scoring Pairs}
\description{
    A multi-class classification problem can be decomposed into a series of
    binary classification problems, for each of which a TSP classifier can be
    built. The predictions of these TSPs (0s or 1s) can be then combined in a
    tree-like structure to predict the original labels.
    
    We use the \code{ctree()} function from the party package for building the
    classification tree.    
}
\usage{
mtsp(X, y, min.score=0.75, max.pairs=10, remap.labels=TRUE, ...)
}
\arguments{
    \item{X}{a data matrix, with features given by columns and samples by rows.}
    \item{y}{a vector of labels, with class labels in sequence: 0,1,2,...}
    \item{min.score}{(optional) the minimum score a TSP has to achieve on a binary
    problem to be further considered.}
    \item{max.pairs}{(optional) maximum number of pairs to be kept from any
    binary problem.}
    \item{remap.labels}{(optional) if labels are not in sequence 0,1,2...,
    then a remapping is needed.}
    
    All extra arguments are passed directly to \code{ctree()} function from the
    \code{party} package.
}
\value{
    An object of class mtsp, with three components:
    \item{$TSP}{the list of all TSPs to be used by the classification tree}
    \item{$TR}{the actual classification tree}
    \item{$labels.map}{the hash table used for remapping the labels}
}
\seealso{
    \code{\link{tsp.s}}, \code{\link{tsp.n}}, \code{\link{predict.mtsp}}
}
\author{Vlad Popovici, \email{vlad.popovici@isb-sib.ch}}
