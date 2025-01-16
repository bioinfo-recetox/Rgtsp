% File man/predict.mtsp.Rd
\name{predict.mtsp}
\alias{predict.mtsp}
\title{Predict function for Multiclass Top Scoring Pairs}
\description{
    Given an mtsp object, predicts the labels for a new data set.
}
\usage{
predict.mtsp(object, X)
}
\arguments{
    \item{object}{an mtsp object.}
    \item{X}{a data matrix, with features given by columns and samples by rows.}
}
\value{
    A vector of labels, with class labels in sequence: 0,1,2,...
}
\seealso{
    \code{\link{mtsp}}
}
\author{Vlad Popovici, \email{vlad.popovici@isb-sib.ch}}
