#ifndef __TSP_HUB_H__
#define __TSP_HUB_H__

#include <vector>
#include <ostream>
#include <list>

#include "ngraph.hpp"

using namespace std;

/**
 * TSP_HUB: a collection of TSPs of the form
 * (i, j1), (i, j2), ..., (i, jn)
 * or (equivalently)
 * (j1, i), (j2, i), ..., (jn, i)
 * such that each pair predicts class "+1" when
 * x[i] > x[j1], x[i] > x[j2], ... (for the first form),
 * or
 * -x[i] > -x[j1], -x[i] > -x[j2],...
 */

class TSP_hub
{
    public:
        unsigned i;             // center of the hub
        vector<unsigned> j;     // all pairing features
        int sign;               // direction of the decision rule
        vector<double> w;       // each pair (i,jk) may have a weight (e.g. its score)
        
    public:
        TSP_hub();
        TSP_hub(const TSP_hub&);
        ~TSP_hub();
        
        void clear();
        
        void center(unsigned _i); // set the "center" of the hub
        unsigned center() const;  // returns the "center" of the hub
        
        void orientation(int s);  // out-going/in-coming edges
        int orientation();
        
        void add_neighbor(unsigned,double=1.0); // adds a new neighbor in the hub
        
	// Predicts either by majority vote or by weighted mean of all the
	// pairs in the hub.
        void predict(const double* const X, unsigned* y, unsigned m, unsigned n,
                     unsigned agg_method=TSP_hub::MAJORITY) const;
        
        TSP_hub& operator =(const TSP_hub&);
        
        friend ostream& operator <<(ostream&, const TSP_hub&);
        
        enum agg_methods { MAJORITY=1, WEIGHTED_MEAN };
        
    protected:
        void predict_majority(const double* X, unsigned* y, unsigned m, unsigned n) const;
        void predict_weighted_mean(const double* X, unsigned* y, unsigned m, unsigned n) const;
};

ostream& operator <<(ostream&, const TSP_hub&);

typedef list<TSP_hub> TSP_hub_list;
typedef map< pair<unsigned,unsigned>, double> edge_weights;

/**
 * Given a network of TSPs, build a list of "hubs". A hub is defined as a
 * subgraph (a multi-tree) in which one node i has either all edges going
 * out or comming in. In terms of TSP, it means that the feature i is either
 * higher or lower than all its neighbors.
 *
 * Parameters:
 *   g - a graph (NGraph) (IN)
 *   hubs - a list of TSP_hub objects (OUT)
 *   min_hub_size - the minimum number of nodes in a hub (IN)
 *
 * Return:
 *   number of hubs found
 */
int build_hubs(NGraph::Graph& g, TSP_hub_list& hubs, unsigned min_hub_size=5);


/**
 * Given a network of TSPs, build a list of "hubs". A hub is defined as a
 * subgraph (a multi-tree) in which one node i has either all edges going
 * out or comming in. In terms of TSP, it means that the feature i is either
 * higher or lower than all its neighbors.
 *
 * Parameters:
 *   g - a graph (NGraph) (IN)
 *   w - an edge_weights object, containing the weights for every pair of features
 *   hubs - a list of TSP_hub objects (OUT)
 *   min_hub_size - the minimum number of nodes in a hub (IN)
 *
 * Return:
 *   number of hubs found
 */
int build_hubs(NGraph::Graph& g, edge_weights& w,
               TSP_hub_list& hubs, unsigned min_hub_size=5);


/**
 * Given a list of hubs, predict the labels of a data set, by aggregating the
 * individual predictions of each hub (majority vote).
 *
 * Parameters:
 *   hubs - a list of TSP_hub objects (IN)
 *   X - data matrix [m x n] with samples by columns (IN)
 *   y - predicted labels vector [n], Caller must manage the memory. (OUT)
 *   m - number of features (IN)
 *   n - number of samples (IN)
 *   yhub - predicted labels matrix: each individual hub generates a vector
 *      of labels, the rows of this matrix. Caller must manage the memory. (OPTIONAL)
 *
 * Return:
 *   0 - if OK
 *  -1 - for error
 */
int meta_gtsp_hub_predict(const TSP_hub_list& hubs,
						  const double* const X, unsigned* y,
						  unsigned m, unsigned n,
                          unsigned* yhub=0,
                          unsigned agg_method=TSP_hub::MAJORITY);


#endif
