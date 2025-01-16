#include "tsp.h"
#include "btypes.h"

#include <list>
#include <cmath>
#include <algorithm>
#include <limits.h>
#include <omp.h>

#include <assert.h>

using namespace std;


TSP::TSP(const double* const _X, const unsigned* const _y,
	 const unsigned _m, const unsigned _n) :
        X(0), y(0), m(_m), n(_n)
{
    X = new double[m*n];
    y = new unsigned[n];
    w0 = new double[n+1];   // last element will hold the sum of the weights
    w1 = new double[n+1];   // last element will hold the sum of the weights
    
    for (unsigned i = 0; i < m*n; ++i) X[i] = _X[i];
    
    n0 = n1 = 0;
    for (unsigned l = 0; l < n; ++l) {
        y[l] = _y[l];
        n1  += y[l];                    // adds 1 for y == +1
        n0  += 1 - y[l];                // adds 1 for y == 0
    }
}

TSP::~TSP()
{
    if (X) delete X;
    if (y) delete y;
    if (w0) delete w0;
    if (w1) delete w1;
}

int TSP::
get_topN(FeatureCollection_dec& final_collection, const unsigned max_pairs) const
{
    /* Algorithm:
     * For each pair of features (i,j), i=0,...,m-1 and j=0,...,m-1, compute
     * the score and return the top scoring S pairs.
     *
     * Instead of using 2 nested loops
     *   for (i=0; i<m; i++) {
     *   	for (j=0; j<i; j++) {
     *   		...compute the score for (i,j)...
     *   	}
     *   }
     *
     * we will use a single loop, which is easy to parallelize:
     *   for (k=0; k<m(m-1)/2; k++) {
     *   	i = some formula of k
     *      j = some formula of k
     *      ...compute the score for (i,j)...
     *   }
     */
    const LU_INT M = m*(m-1L)/2L;    

#ifdef _OPENMP
    const unsigned NC = omp_get_max_threads();
#else
    const unsigned NC = 1;
#endif
    // use several collections, to reduce concurrency:
    FeatureCollection_inc fcoll[NC]; 

#ifdef _OPENMP	
    omp_lock_t locks[NC];             // access to each list is controlled by a lock
	
    for (unsigned i = 0; i < NC; ++i) {
	omp_init_lock(&locks[i]);
    }
#endif

    
#pragma omp parallel shared(fcoll,locks)
{
    #pragma omp for 
    for (LU_INT k = 0; k < M; ++k) {
        unsigned i = static_cast<unsigned>(floor((sqrt(1.0+8.0*k)+1.0)/2.0));
	unsigned j = static_cast<unsigned>(k - i*(i-1L)/2L);

	double sc = score(i,j);             // compute the score for (i,j)
        
	unsigned lk = k % NC;			
#ifdef _OPENMP
	omp_set_lock(&locks[lk]);           // ...CRITICAL STARTS
#endif
	double smallest_score = (*fcoll[lk].begin()).score;
	if ((fcoll[lk].size() < max_pairs) || (sc > smallest_score)) {
            fcoll[lk].insert(FeaturePair(sc, i, j));
            if (fcoll[lk].size() > max_pairs) {
                fcoll[lk].erase(FeaturePair(smallest_score, 0, 0));
            }
	}
#ifdef _OPENMP
	omp_unset_lock(&locks[lk]);         // ...CRITICAL ENDS
#endif

    }
} //#pragma omp

#ifdef _OPENMP
    for (unsigned i = 0; i < NC; ++i) 
	omp_destroy_lock(&locks[i]);
#endif

    final_collection.clear();
    for (unsigned i = 0; i < NC; ++i) {
	final_collection.insert(fcoll[i].begin(), fcoll[i].end());
    }

    while (final_collection.size() > max_pairs) {  // at most max_pairs are kept
	// remove all smallest scores (at the end of the list):
        FeatureCollection_dec::iterator last = final_collection.end();
        last--;
	final_collection.erase(FeaturePair((*last).score, 0, 0));
    }
	
    return final_collection.size();
}

int TSP::
get_topN(FeatureCollection_dec& final_collection,
         const double* const weights,
         const unsigned max_pairs) const
{
    const LU_INT M = m*(m-1L)/2L;    

#ifdef _OPENMP
    const unsigned NC = omp_get_max_threads();
#else
    const unsigned NC = 1;
#endif
    // use several collections, to reduce concurrency:
    FeatureCollection_inc fcoll[NC]; 

#ifdef _OPENMP	
    omp_lock_t locks[NC];             // access to each list is controlled by a lock
	
    for (unsigned i = 0; i < NC; ++i) {
	omp_init_lock(&locks[i]);
    }
#endif

    // prepare the weights for each class:
    w0[n] = w1[n] = 0.0;
    for (unsigned i = 0; i < n; ++i) {
        w0[i] = (1 - y[i]) * weights[i];
        w1[i] =      y[i]  * weights[i];
        w0[n]+= w0[i];
        w1[n]+= w1[i];
    }
    
#pragma omp parallel shared(fcoll,locks)
{
    #pragma omp for 
    for (LU_INT k = 0; k < M; ++k) {
        unsigned i = static_cast<unsigned>(floor((sqrt(1.0+8.0*k)+1.0)/2.0));
	unsigned j = static_cast<unsigned>(k - i*(i-1L)/2L);

	double sc = wscore(i, j);   // compute the score for (i,j)
        
	unsigned lk = k % NC;			
#ifdef _OPENMP
	omp_set_lock(&locks[lk]);           // ...CRITICAL STARTS
#endif
	double smallest_score = (*fcoll[lk].begin()).score;
	if ((fcoll[lk].size() < max_pairs) || (sc > smallest_score)) {
            fcoll[lk].insert(FeaturePair(sc, i, j));
            if (fcoll[lk].size() > max_pairs) {
                fcoll[lk].erase(FeaturePair(smallest_score, 0, 0));
            }
	}
#ifdef _OPENMP
	omp_unset_lock(&locks[lk]);         // ...CRITICAL ENDS
#endif

    }
} //#pragma omp

#ifdef _OPENMP
    for (unsigned i = 0; i < NC; ++i) 
	omp_destroy_lock(&locks[i]);
#endif

    final_collection.clear();
    for (unsigned i = 0; i < NC; ++i) {
	final_collection.insert(fcoll[i].begin(), fcoll[i].end());
    }

    while (final_collection.size() > max_pairs) {  // at most max_pairs are kept
	// remove all smallest scores (at the end of the list):
        FeatureCollection_dec::iterator last = final_collection.end();
        last--;
	final_collection.erase(FeaturePair((*last).score, 0, 0));
    }
	
    return final_collection.size();
}

int TSP::
get_topS(FeatureCollection_dec& final_collection, const double min_score) const
{
    /* Algorithm:
     * For each pair of features (i,j), i=0,...,m-1 and j=0,...,m-1, compute
     * the score and return the top scoring S pairs.
     *
     * Instead of using 2 nested loops
     *   for (i=0; i<m; i++) {
     *   	for (j=0; j<i; j++) {
     *   		...compute the score for (i,j)...
     *   	}
     *   }
     *
     * we will use a single loop, which is easy to parallelize:
     *   for (k=0; k<m(m-1)/2; k++) {
     *   	i = some formula of k
     *      j = some formula of k
     *      ...compute the score for (i,j)...
     *   }
     */
    const LU_INT M = m*(m-1L)/2L;    

#ifdef _OPENMP
    const unsigned NC = omp_get_max_threads();
#else
    const unsigned NC = 1;
#endif
    // use several collections, to reduce concurrency:
    FeatureCollection_inc fcoll[NC]; 

#ifdef _OPENMP	
    omp_lock_t locks[NC];             // access to each list is controlled by a lock
	
    for (unsigned i = 0; i < NC; ++i) {
	omp_init_lock(&locks[i]);
    }
#endif

    
#pragma omp parallel shared(fcoll,locks)
{
    #pragma omp for 
    for (LU_INT k = 0; k < M; ++k) {
        unsigned i = static_cast<unsigned>(floor((sqrt(1.0+8.0*k)+1.0)/2.0));
	unsigned j = static_cast<unsigned>(k - i*(i-1L)/2L);

	double sc = score(i,j);             // compute the score for (i,j)
        
	unsigned lk = k % NC;			
#ifdef _OPENMP
	omp_set_lock(&locks[lk]);           // ...CRITICAL STARTS
#endif
	if (sc >= min_score) {
            fcoll[lk].insert(FeaturePair(sc, i, j));
	}
#ifdef _OPENMP
	omp_unset_lock(&locks[lk]);         // ...CRITICAL ENDS
#endif

    }
} //#pragma omp

#ifdef _OPENMP
    for (unsigned i = 0; i < NC; ++i) 
	omp_destroy_lock(&locks[i]);
#endif

    final_collection.clear();
    for (unsigned i = 0; i < NC; ++i) {
	final_collection.insert(fcoll[i].begin(), fcoll[i].end());
    }
	
    return final_collection.size();
}

int TSP::
get_topS(FeatureCollection_dec& final_collection,
         const double* const weights,
         const double min_score) const
{
    const LU_INT M = m*(m-1L)/2L;    

#ifdef _OPENMP
    const unsigned NC = omp_get_max_threads();
#else
    const unsigned NC = 1;
#endif
    // use several collections, to reduce concurrency:
    FeatureCollection_inc fcoll[NC]; 

#ifdef _OPENMP	
    omp_lock_t locks[NC];             // access to each list is controlled by a lock
	
    for (unsigned i = 0; i < NC; ++i) {
	omp_init_lock(&locks[i]);
    }
#endif

    // prepare the weights for each class:
	w0[n] = w1[n] = 0.0;
    for (unsigned i = 0; i < n; ++i) {
		w0[i] = (1 - y[i]) * weights[i];
		w1[i] =      y[i]  * weights[i];
		w0[n]+= w0[i];
		w1[n]+= w1[i];
    }

#pragma omp parallel shared(fcoll,locks)
{
    #pragma omp for 
    for (LU_INT k = 0; k < M; ++k) {
        unsigned i = static_cast<unsigned>(floor((sqrt(1.0+8.0*k)+1.0)/2.0));
	unsigned j = static_cast<unsigned>(k - i*(i-1L)/2L);

	double sc = wscore(i, j);           // compute the score for (i,j)
        
	unsigned lk = k % NC;			
#ifdef _OPENMP
	omp_set_lock(&locks[lk]);           // ...CRITICAL STARTS
#endif
	if (sc >= min_score) {
            fcoll[lk].insert(FeaturePair(sc, i, j));
	}
#ifdef _OPENMP
	omp_unset_lock(&locks[lk]);         // ...CRITICAL ENDS
#endif

    }
} //#pragma omp

#ifdef _OPENMP
    for (unsigned i = 0; i < NC; ++i) 
	omp_destroy_lock(&locks[i]);
#endif

    final_collection.clear();
    for (unsigned i = 0; i < NC; ++i) {
	final_collection.insert(fcoll[i].begin(), fcoll[i].end());
    }
	
    return final_collection.size();
}


double TSP::score(unsigned& i, unsigned& j) const
{
    const double *x1, *x2;
    unsigned p1, p0;
    double s, dp1, dp0;

    x1 = feature(i);                       // the i-th feature (row)
    x2 = feature(j);                       // the j-th feature (row)
    p1 = p0 = 0;
    for (unsigned l = 0; l < n; ++l) {
        p1 +=      y[l]  * (x1[l] < x2[l]);
        p0 += (1 - y[l]) * (x1[l] < x2[l]);
    }
    
    dp1 = static_cast<double>(p1)/n1;
    dp0 = static_cast<double>(p0)/n0;
    // Make sure that score(feature(i)) > score(feature(j)) -> class "+1":
    if (dp1 >= dp0) {
        s = dp1 - dp0;
    }
    else {
        s = dp0 - dp1;
        swap(i, j);
    }
    
    return s;
}

double TSP::wscore(unsigned& i, unsigned& j) const
{
    const double *x1, *x2;
    double s, dp1, dp0;

    x1 = feature(i);                      // the i-th feature (row)
    x2 = feature(j);                      // the j-th feature (row)
    dp1 = dp0 = 0.0;
    
    for (unsigned l = 0; l < n; ++l) {
        dp1 += w1[l] * (x1[l] < x2[l]);
        dp0 += w0[l] * (x1[l] < x2[l]);
    }
    
    dp1 /= w1[n];
    dp0 /= w0[n];
    
    // Make sure that score(feature(i)) < score(feature(j)) -> class "+1":
    if (dp1 >= dp0) {
        s = dp1 - dp0;
    }
    else {
        s = dp0 - dp1;
        swap(i, j);
    }
    
    return s;
}


void
TSP_Predictor::
predict(const double* const X, unsigned* y, unsigned m, unsigned n, double w) const
{
	const double* xi = static_cast<const double*>(X + p.i*n);
	const double* xj = static_cast<const double*>(X + p.j*n);

	assert(p.i < m);
	assert(p.j < m);
	
	for (unsigned k=0; k < n; ++k) {
		y[k] = (xi[k] < w*xj[k] ? 1 : 0);
	}
}
        
// Predict the label of one vector.
// x - data vector of m features
unsigned
TSP_Predictor::predict(const double* const x, unsigned m, double w) const
{
	assert(p.i < m);
	assert(p.j < m);

	return (x[p.i] < w*x[p.j] ? 1 : 0);
}
