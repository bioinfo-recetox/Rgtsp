\name{Rgtsp-package}
\alias{Rgtsp-package}
\alias{Rgtsp}
\docType{package}
\title{
Generalized Top Scoring Pairs in R
}
\description{
This package implements the Top Scoring Pairs algorithm and a number of
    extensions, to allow different scoring criteria and sample weighting. All
    the functions are designed to take full advantage of the multicore/
    multi-processor capabilities of the host machine, by using the OpenMP
    library. This means that the package needs to be compiled with a reasonably
    recent C++ compiler.
    
    The basic idea is to score all the feature pairs (i,j) and to return the
    top ones (those with a highest score). There are many ways in which the
    score can be computed:
        -by counting how may times X(:,i) > X(:,j) in the first class and
        X(:,i) < X(:,j) in the second class (as in Geman et al.)
        -by fitting a bi-variate logistic regression (with features i and j as
        predictors) and estimate its suitabilty for class discrimination
        
    All these scores accept also sample weights, so the pairs can be used in a
    meta-classifier setting.
    
    While the package is designed to be used for binary classification problems,
    by implementing different scoring functions, one can use the same machinery
    for other pairwise interactions of the features. This, however, requires
    extending the C++ code.
}
\details{
\tabular{ll}{
Package: \tab Rgtsp\cr
Version: \tab 0.0-2\cr
Date: \tab 2010-04-27\cr
License: \tab GPL (>= 2)\cr
Packaged: \tab 2010-04-28, vlad\cr
Built: \tab R 2.10.1; x86_64-unknown-linux-gnu; 2010-04-29 11:59:10 UTC; unix\cr
}

Index:
\preformatted{
                        Generalized Top Scoring Pairs in R.
tsp.n                   Top Scoring N Pairs
tsp.s                   Top Scoring Pairs with a minimum guaranteed
                        score.
}
}
\author{
Vlad Popovici <vlad.popovici@isb-sib.ch>

Maintainer: Vlad Popovici <vlad.popovici@isb-sib.ch>
}
\references{
    Geman, D. and d'Avignon, C. and Naiman, D.Q. and
    Winslow, R.L. (2004) Classifying Gene Expression
    Profiles from Pairwise mRNA Comparisons,
    \emph{Stat Appl Genet Mol Biol}, \bold{3}, Article 19
}
\keyword{ package }
\seealso{
\code{\link[ts.pair:tspair-package]{tspair}} Package
}
\examples{
}
