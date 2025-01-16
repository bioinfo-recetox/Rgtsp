% File man/tsp.hub.Rd
\name{tsp.hub}
\alias{tsp.hub}
\title{Hubs of Top Scoring Pairs}
\description{
    Given a set of top scoring pairs, finds 'hubs' of features. A hub is defined
    as a 1-level multi-tree, in which one feature (the center of the hub) has
    either a higher or a lower value than all its neighbors for one class, and
    the reverse for the second class.
    
    A hub of TSPs has the form \{(i, j_1, s_1),..., (i, j_n, s_n)\}, each triplet
    predicting the label of a sample by the rule
    
    if s_k*X[,i] > s_k*X[,j_k] then predict "+1"; othewise predict "-1"
    
    s_k are +/-1s and they are used to indicate the 'direction' of the decision
    in the above rule. These individual decisions are combined by majority vote.
}
\usage{
tsp.hub(tsp, min.hub.size=5)
}
\arguments{
    \item{tsp}{a TSP object, as returned by tsp.n() or tsp.s() functions}
    \item{min.hub.size}{(optional) minimum number of pairs in a hub.}
}
\value{
    A list of the form \{(center_1, sign_1, neighbors_1), ...,
    (center_p, sign_p, neighbors_p)\}
}
\seealso{
    \code{\link{tsp.s}}, \code{\link{tsp.n}}
}
\author{Vlad Popovici, \email{vlad.popovici@isb-sib.ch}}
