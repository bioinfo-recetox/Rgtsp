% File man/Rgtsp.Rd
\name{Rgtsp}
\title{Generalized Top Scoring Pairs in R}
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
\seealso{
    \code{\link{tsp.s}}, \code{\link{tsp.n}}
}
\references{
    D. Geman, C. d'Avignon, D.Q. Naiman, R.L. Winslow, "Classifying
    Gene Expression Profiles from Pairwise mRNA Comparisons",
    Stat Appl Genet Mol Biol, 2004:3:Article 19
}
\author{Vlad Popovici, \email{vlad.popovici@isb-sib.ch}}