#include "rgtsp.h"
#include "tsp.h"
#include "btypes.h"
#include "tsp_hub.h"
#include "ngraph.hpp"


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <omp.h>
#include <iostream>

using namespace std;

/*
 * Classical TSP
 */
extern "C"
void RGTSP_tsp_N(double* X, unsigned* y, unsigned* nrow_X, unsigned* ncol_X,
                 unsigned* i, unsigned* j, double* s, unsigned* max_np)
{
	TSP* worker = new TSP(X, y, *nrow_X, *ncol_X);
    FeatureCollection_dec fc;
    unsigned k;
    
    worker->get_topN(fc, *max_np);
	
    for (k = 0; k < *max_np; ++k) { i[k] = j[k] = 0; s[k] = -1.0; }
    
    k = 0;
    for(FeatureCollection_dec::const_iterator p=fc.begin(); p != fc.end(); ++p) {
		i[k] = (*p).i+1;
        j[k] = (*p).j+1;
        s[k] = (*p).score;
        ++k;
    }
    *max_np = k;     // the total number of pairs actually returned
	
	delete worker;
}


extern "C"
void RGTSP_tsp_NW(double* X, unsigned* y, double* w,
				  unsigned* nrow_X, unsigned* ncol_X,
				  unsigned* i, unsigned* j, double* s, unsigned* max_np)
{
	TSP* worker = new TSP(X, y, *nrow_X, *ncol_X);
    FeatureCollection_dec fc;
    unsigned k;
    
    worker->get_topN(fc, w, *max_np);
	
    for (k = 0; k < *max_np; ++k) { i[k] = j[k] = 0; s[k] = -1.0; }
    
    k = 0;
    for(FeatureCollection_dec::const_iterator p=fc.begin(); p != fc.end(); ++p) {
		i[k] = (*p).i+1;
        j[k] = (*p).j+1;
        s[k] = (*p).score;
        ++k;
    }
    *max_np = k;     // the total number of pairs actually returned
	
	delete worker;
}


extern "C"
void RGTSP_tsp_S(double* X, unsigned* y, unsigned* nrow_X, unsigned* ncol_X,
                 double* min_score, unsigned* i, unsigned* j, double* s,
                 unsigned* np)
{
	TSP* worker = new TSP(X, y, *nrow_X, *ncol_X);
    FeatureCollection_dec fc;
    unsigned k;
    
    worker->get_topS(fc, *min_score);
    
    for (k = 0; k < *np; ++k) { i[k] = j[k] = 0; s[k] = -1.0; }
    
    k = 0;
    for(FeatureCollection_dec::const_iterator p=fc.begin(); p != fc.end(); ++p) {
        if (k >= *np) break;
		i[k] = (*p).i+1;
        j[k] = (*p).j+1;
        s[k] = (*p).score;
        ++k;
    }
    *np = k;     // the total number of pairs actually returned
	
	delete worker;
}


extern "C"
void RGTSP_tsp_SW(double* X, unsigned* y, double* w,
				  unsigned* nrow_X, unsigned* ncol_X,
				  double* min_score, unsigned* i, unsigned* j, double* s,
				  unsigned* np)
{
	TSP* worker = new TSP(X, y, *nrow_X, *ncol_X);
    FeatureCollection_dec fc;
    unsigned k;
    
    worker->get_topS(fc, w, *min_score);
    
    for (k = 0; k < *np; ++k) { i[k] = j[k] = 0; s[k] = -1.0; }
    
    k = 0;
    for(FeatureCollection_dec::const_iterator p=fc.begin(); p != fc.end(); ++p) {
        if (k >= *np) break;
		i[k] = (*p).i+1;
        j[k] = (*p).j+1;
        s[k] = (*p).score;
        ++k;
    }
    *np = k;     // the total number of pairs actually returned
	
	delete worker;
}

///////
// X: matrix nrow_X x ncol_X [in]
// y: matrix np x ncol_X [out]
// nrow_X: no. of features
// ncol_X: no. of samples
// i, j: vectors defining the TSPs
// np: no. of pairs
//
// All memory is supposed to be allocated outside this function, by the caller,
// as contiguous blocks.
extern "C"
void RGTSP_tsp_predict(double* X, unsigned* y,
                       unsigned* nrow_X, unsigned* ncol_X,
                       unsigned* i, unsigned* j, unsigned* np)
{
	for (unsigned k=0; k < *np; ++k) {
		TSP_Predictor p(i[k], j[k]);
		p.predict(X, &y[k*(*ncol_X)], *nrow_X, *ncol_X);
	}
}

///////
// X: matrix nrow_X x ncol_X [in]
// y: matrix np x ncol_X [out]
// w: a vector of np elements, weighting the features j's
// nrow_X: no. of features
// ncol_X: no. of samples
// i, j: vectors defining the TSPs
// np: no. of pairs
//
// All memory is supposed to be allocated outside this function, by the caller,
// as contiguous blocks.
extern "C"
void RGTSP_tsp_predict_W(double* X, unsigned* y, double* w,
                       unsigned* nrow_X, unsigned* ncol_X,
                       unsigned* i, unsigned* j, unsigned* np)
{
	for (unsigned k=0; k < *np; ++k) {
		TSP_Predictor p(i[k], j[k]);
		p.predict(X, &y[k*(*ncol_X)], (*nrow_X), (*ncol_X), w[k]);
	}	
}


extern "C"
void RGTSP_tsp_hub(unsigned* i, unsigned* j, double* w, unsigned* np,
				   unsigned* hc, unsigned* hn, int* hs, unsigned* m,
				   unsigned* minhs)
{
	NGraph::Graph g;
	edge_weights weights;
	TSP_hub_list hl;
	unsigned nh;
	
	// construct the graph
	for (unsigned k=0; k < *np; ++k) {
		g.insert_edge(j[k],i[k]); // i < j -> class "+1"
		weights[pair<unsigned,unsigned>(j[k],i[k])] = w[k];
	}
	
	nh = build_hubs(g, weights, hl, *minhs);
	
	if (*m == 0) return;  // found nothing 
	// copy the hubs back:
	unsigned k = 0;
	for (TSP_hub_list::const_iterator h=hl.begin(); h != hl.end(); h++) {
		if (k >= *m) break;        // at most m pairs in the output
		for (unsigned l = 0; l < (*h).j.size(); ++l) {
			hc[k] = (*h).i;       // hub center
			hn[k] = (*h).j[l];
			hs[k] = (*h).sign;
			++k;
			if (k >= *m) break;    // at most m pairs in the output
		}			
	}
	*m = k;
}




/*
 * DLL-related functions.
 */
extern "C"
void R_init_rgtsp(DllInfo* info)
{
//	fprintf(stderr, "Loading Rgtsp.\n");

#ifdef _OPENMP
	fprintf(stderr, "Using OpenMP!\n");
//	fprintf(stderr, "\t---> Maximum number of threads available: %d\n", omp_get_max_threads());
//	fprintf(stderr, "\t---> Number of processing units: %d\n", omp_get_num_procs());
	omp_set_dynamic(1);
//	fprintf(stderr, "\t---> Using %s thread allocation strategy.\n",
//			(omp_get_dynamic() ? "dynamic" : "static"));
#else
	fprintf(stderr, "No OpenMP.\n");
#endif
	fflush(NULL);
}

extern "C"
void R_unload_mylib(DllInfo *info)
{
}

