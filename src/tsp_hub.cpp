#include "tsp_hub.h"
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <ostream>
#include <cstring>

using namespace std;

void TSP_hub::center(unsigned _i)
{
    i = _i;
}

unsigned TSP_hub::center() const
{
    return i;
}
        
void TSP_hub::orientation(int s)
{
    sign = s;
}

int TSP_hub::orientation()
{
    return sign;
}
        
TSP_hub::TSP_hub()  : i(0), sign(1)
{
    j.reserve(10);
    w.reserve(10);
}

TSP_hub::~TSP_hub() {}

TSP_hub::TSP_hub(const TSP_hub& other)
{
    *this = other;
}

void TSP_hub::clear()
{
    i = 0; sign = 1;
    j.clear();
    w.clear();
}

TSP_hub& TSP_hub::operator =(const TSP_hub& other)
{
    if (this != &other) {
        i = other.i;
        sign = other.sign;
        j.clear();
        if (!other.j.empty()) {
            for (vector<unsigned>::const_iterator p = other.j.begin(); p != other.j.end(); p++) {
                j.push_back(*p);
            }
        }
        w.clear();
        if (!other.w.empty()) {
            for (vector<double>::const_iterator p = other.w.begin(); p != other.w.end(); p++) {
                w.push_back(*p);
            }
        }
    }
    
    return *this;
}


void TSP_hub::add_neighbor(unsigned _j, double _w)
{
    j.push_back(_j);
    w.push_back(_w);
}

void TSP_hub::predict(const double* const X, unsigned* y, unsigned m, unsigned n,
                     unsigned agg_method) const
{
    switch(agg_method) {
        case (MAJORITY):
            predict_majority(X, y, m, n);
            break;
        case (WEIGHTED_MEAN):
            predict_weighted_mean(X, y, m, n);
            break;
        default: // should use some exceptions
            cerr << "Unknown aggregation method!\n";
            exit(1);
    }
    
}

void TSP_hub::predict_majority(const double* X, unsigned* y, unsigned m, unsigned n) const
{
    if (j.empty()) return;  // nothing to do

	assert(i < m);
    vector<unsigned> c1(n); // count the votes for class "+1"
    const double* xi = static_cast<const double*>(X + i*n);
    
    for (unsigned k=0; k < j.size(); ++k) {
        assert(j[k] < m);
        const double* xj = static_cast<const double*>(X + j[k]*n);
    
    	for (unsigned l=0; l < n; ++l) {
        	c1[l] += (sign*xi[l] > sign*xj[l] ? 1 : 0);
        }
    }
    
    for (unsigned l=0; l < n; ++l) {
        y[l] = 0;
        if (c1[l] >= j.size()/2+1) y[l] = 1;
    }
}

void TSP_hub::predict_weighted_mean(const double* X, unsigned* y, unsigned m, unsigned n) const
{
    if (j.empty()) return;  // nothing to do
    if (w.empty()) return;

	assert(i < m);
    vector<double> wm(n); // weighted mean
    const double* xi = static_cast<const double*>(X + i*n);
    double sw = 0.0;      // sum of all weights
    
    for (unsigned k=0; k < j.size(); ++k) {
        assert(j[k] < m);
        const double* xj = static_cast<const double*>(X + j[k]*n);
    
    	for (unsigned l=0; l < n; ++l) {
        	wm[l] += (sign*xi[l] > sign*xj[l] ? w[k] : 0);
        }
        sw += w[k];
    }
    
    for (unsigned l=0; l < n; ++l) {
        y[l] = 0;
        if (wm[l] >= sw / 2.0) y[l] = 1;
    }
}

ostream& operator <<(ostream& out, const TSP_hub& h)
{
    out << h.i;
    out << (h.sign > 0 ? " > " : " < ");
    for (vector<unsigned>::const_iterator k=h.j.begin(); k != h.j.end(); k++)
        out << (*k) << ' ';
    
    return out;
}


int meta_gtsp_hub_predict(const TSP_hub_list& hubs,
						  const double* const X, unsigned* y,
						  unsigned m, unsigned n,
                          unsigned* yhub,
                          unsigned agg_method)
{
	if (hubs.empty()) return -1;
	
	unsigned* ytmp = new unsigned[n*hubs.size()];
	unsigned i, j, k;
	
	k = 0;
    for (TSP_hub_list::const_iterator h=hubs.begin(); h != hubs.end(); h++) {
        (*h).predict(X, &ytmp[k], m, n, agg_method);
        k += n;
    }
		
	for (i = 0; i < n; ++i) {
		unsigned s = 0;
		for (k = 0; k < hubs.size(); ++k) {
			s += ytmp[i+k*n];
		}
		y[i] = 0;
		if (s >= hubs.size()/2+1) y[i] = 1;
	}
	if (yhub) {         // copy the individual predictions of each hub:
        memcpy(yhub, ytmp, n*hubs.size()*sizeof(unsigned));
    }
    
	delete[] ytmp;
	
	return 0;
}

int build_hubs(NGraph::Graph& g, TSP_hub_list& hubs, unsigned min_hub_size)
{
	NGraph::Graph gtmp;
	
	unsigned int max_in_degree = 0, max_out_degree = 0;
	unsigned int in_max = 0, out_max = 0;
	TSP_hub hub;
	
	// find all 'outer-hubs' (i > j1,j2,...)	
	gtmp = g;
	
	while (gtmp.num_edges() > 0) {
		max_out_degree = 0;
		for (NGraph::Graph::const_iterator p=gtmp.begin(); p != gtmp.end(); p++) {
			unsigned int k = gtmp.node(p);
			if (gtmp.out_degree(k) > max_out_degree) {
				max_out_degree = gtmp.out_degree(k);
				out_max = k;
			}
		}
		if (max_out_degree < min_hub_size) break;
		
		NGraph::Graph::vertex_set vs = gtmp.out_neighbors(out_max);
		
		hub.clear();
		hub.center(out_max);
		for (NGraph::Graph::vertex_set::const_iterator p=vs.begin(); p != vs.end(); p++) {
			hub.add_neighbor(*p);
			gtmp.remove_edge(out_max, *p);
		}
		hubs.push_back(hub);
	}

	// find all 'inner-hubs' (i < j1,j2,...)	
	gtmp = g;
	
	while (gtmp.num_edges() > 0) {
		max_in_degree = 0;
		for (NGraph::Graph::const_iterator p=gtmp.begin(); p != gtmp.end(); p++) {
			unsigned int k = gtmp.node(p);
			if (gtmp.in_degree(k) > max_in_degree) {
				max_in_degree = gtmp.in_degree(k);
				in_max = k;
			}
		}
		if (max_in_degree < min_hub_size) break;
		
		NGraph::Graph::vertex_set vs = gtmp.in_neighbors(in_max);
		TSP_hub hub;
		hub.center(in_max);
		hub.orientation(-1);
		for (NGraph::Graph::vertex_set::const_iterator p=vs.begin(); p != vs.end(); p++) {
			hub.add_neighbor(*p);
			gtmp.remove_edge(*p, in_max);
		}
		hubs.push_back(hub);
	}

	return hubs.size();
}


int build_hubs(NGraph::Graph& g, edge_weights& w,
               TSP_hub_list& hubs, unsigned min_hub_size)
{
	NGraph::Graph gtmp;
	
	unsigned int max_in_degree = 0, max_out_degree = 0;
	unsigned int in_max = 0, out_max = 0;
	TSP_hub hub;
	
	// find all 'outer-hubs' (i > j1,j2,...)	
	gtmp = g;
	
	while (gtmp.num_edges() > 0) {
		max_out_degree = 0;
		for (NGraph::Graph::const_iterator p=gtmp.begin(); p != gtmp.end(); p++) {
			unsigned int k = gtmp.node(p);
			if (gtmp.out_degree(k) > max_out_degree) {
				max_out_degree = gtmp.out_degree(k);
				out_max = k;
			}
		}
		if (max_out_degree < min_hub_size) break;
		
		NGraph::Graph::vertex_set vs = gtmp.out_neighbors(out_max);
		
		hub.clear();
		hub.center(out_max);
		for (NGraph::Graph::vertex_set::const_iterator p=vs.begin(); p != vs.end(); p++) {
			hub.add_neighbor(*p, w[pair<unsigned,unsigned>(out_max, *p)]);
			gtmp.remove_edge(out_max, *p);
		}
		hubs.push_back(hub);
	}

	// find all 'inner-hubs' (i < j1,j2,...)	
	gtmp = g;
	
	while (gtmp.num_edges() > 0) {
		max_in_degree = 0;
		for (NGraph::Graph::const_iterator p=gtmp.begin(); p != gtmp.end(); p++) {
			unsigned int k = gtmp.node(p);
			if (gtmp.in_degree(k) > max_in_degree) {
				max_in_degree = gtmp.in_degree(k);
				in_max = k;
			}
		}
		if (max_in_degree < min_hub_size) break;
		
		NGraph::Graph::vertex_set vs = gtmp.in_neighbors(in_max);
		TSP_hub hub;
		hub.center(in_max);
		hub.orientation(-1);
		for (NGraph::Graph::vertex_set::const_iterator p=vs.begin(); p != vs.end(); p++) {
            hub.add_neighbor(*p, w[pair<unsigned,unsigned>(*p, in_max)]);
			gtmp.remove_edge(*p, in_max);
		}
		hubs.push_back(hub);
	}

	return hubs.size();
}

