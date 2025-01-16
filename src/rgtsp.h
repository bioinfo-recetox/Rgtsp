/*
 * Title      : rgtsp.h
 * Author     : Vlad Popovici
 * Description: C-wrappers for the C++ functions, intended to be called from R.
 */
#ifndef __RGTSP_H__
#define __RGTSP_H__

/*
 * Classical Top Scoring Pairs algorithm: returns either the <max_np> top
 * scoring pairs (_N), or all the pairs above a minimal score (_S).
 */
// Return the max_np top scoring pairs, using the classical scoring function
// (as in Geman et al.)
extern "C"
void RGTSP_tsp_N(double* X, unsigned* y, unsigned* nrow_X, unsigned* ncol_X,
                 unsigned* i, unsigned* j, double* s, unsigned* max_np);
extern "C"
void RGTSP_tsp_NW(double* X, unsigned* y, double* w,
				  unsigned* nrow_X, unsigned* ncol_X,
				  unsigned* i, unsigned* j, double* s, unsigned* max_np);

// Return the feature pairs scoring at least min_score, but no more than
// np pairs.
extern "C"
void RGTSP_tsp_S(double* X, unsigned* y, unsigned* nrow_X, unsigned* ncol_X,
                 double* min_score, unsigned* i, unsigned* j, double* s,
                 unsigned* np);
extern "C"
void RGTSP_tsp_SW(double* X, unsigned* y, double* w,
				  unsigned* nrow_X, unsigned* ncol_X,
				  double* min_score, unsigned* i, unsigned* j, double* s,
				  unsigned* np);

// Given a list of TSPs {(i[k],j[k]), k=1...np}, predict the labels for a
// data set X (samples by columns) and return them as a matrix y of np
// rows and ncol_X columns.
extern "C"
void RGTSP_tsp_predict(double* X, unsigned* y,
                       unsigned* nrow_X, unsigned* ncol_X,
                       unsigned* i, unsigned* j, unsigned* np);

// Given a list of TSPs {(i[k],j[k],w[k]), k=1...np}, predict the labels for a
// data set X (samples by columns) and return them as a matrix y of np
// rows and ncol_X columns. The weighted decision rule is:
//
// predict label "+1" when i[k] < w[k]*j[k]
//
extern "C"
void RGTSP_tsp_predict_W(double* X, unsigned* y, double* w,
                       unsigned* nrow_X, unsigned* ncol_X,
                       unsigned* i, unsigned* j, unsigned* np);

// Given a list of TSPs construct a number of hubs.
// TSPs: {(i[k], j[k]), k=1,...,np}
//
// The hubs are returned in hc, hn, hs vectors, with at most m (IN/OUT)
// elements:
// hc - the centers of the hubs
// hn - the neighbors in a hub
// hs - the signs
// The 3 vectors (m=11) may look something like:
//  hc = 1,1,1,1,1, 2, 2, 2, 2,10,10
//  hn = 3,4,5,6,7, 3, 5, 8, 9, 7,11
//  hs = 1,1,1,1,1,-1,-1,-1,-1, 1, 1
// for 3 hubs of the form:
//   1 > {3,4,5,6,7}
//   2 < {3,5,8,9}
//  10 > {7,11}
extern "C"
void RGTSP_tsp_hub(unsigned* i, unsigned* j, double* w, unsigned* np,
				   unsigned* hc, unsigned* hn, int* hs, unsigned* m,
                   unsigned* minhs);


#endif
