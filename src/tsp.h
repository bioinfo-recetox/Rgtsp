#ifndef __TSP_H__
#define __TSP_H__

#include "btypes.h"

/**
 * Top Scoring Pairs: a simple decision rule that decides for one class or the
 * other based on the relative order of pairs of features: if F1 < F2, then
 * decide class "1", otherwise decide class "0".
 * For details see Geman et al, Classifying gene expression profiles from
 * pairwise mRNA comparisons, Stat Appl Genet Mol Biol, 2004 3(19)
 */

/**
 * Data representation:
 * -generally, all vector and matrices are stored as a contiguous block of
 * memory;
 * -data matrices (denoted by capital letters) have the samples by columns and
 * features by rows;
 * -labels: only binary classification is handled, with labels coded as 0 and 1
 */

/**
 * TSP: a class that encapsulates the algorithm for scanning all pairwise
 * interactions. The strength of interaction is given by the score() function
 * and a higher score is supposed to indicate a more valuable pair of features.
 * For classification purposes, Geman's TSP algorithm is implemented. Other
 * criteria can be implemented by overriding the virtual method score().
 * 
 * Two methods for selecting the significant pairs are provided:
 *   get_topN()   - returns top N pairs
 *   get_topS()   - returns all pairs with a score above a threshold
 */
class TSP
{
    protected:
        double* X;                     // data matrix, features by rows
        unsigned* y;                   // labels (0/1)
        const unsigned m, n;           // no. of features (m) and samples (n)
        unsigned n0, n1;               // sample size for each class
        double* w0;                    // weights for class '0' or 0 for elements in class '1'
        double* w1;                    // weights for class '1' or 0 for elements in class '0'
        
    public:
        TSP(const double* const _X, const unsigned* const _y,
            const unsigned _m, const unsigned _n);
        
        virtual ~TSP();
        
        // Return top <max_pairs> pairs (with weights, eventually)
        virtual int get_topN(FeatureCollection_dec& final_collection,
                             const unsigned max_pairs=1) const;
        virtual int get_topN(FeatureCollection_dec& final_collection,
                              const double* const weights,
                              const unsigned max_pairs=1) const;
        
        // Return all pairs with a score above <min_score> (with weights, eventually)
        virtual int get_topS(FeatureCollection_dec& final_collection,
                             const double min_score=0.5) const;
        virtual int get_topS(FeatureCollection_dec& final_collection,
                             const double* const weights,
                             const double min_score=0.5) const;
                
    protected:
        // Compute the score for the feature pair (i,j) and exchange i with j
        // if the score would be negative.
        virtual double score(unsigned& i, unsigned& j) const;
        
        // Weighted version of the score function: each sample is allowed to
        // have a weight.
        virtual double wscore(unsigned& i, unsigned& j) const;
        
        virtual const double* feature(unsigned i) const {
            return static_cast<const double*>(X + i*n);
        }
};

/**
 * TSP_Predictor: a simple wrapper class for prediction functions. This class
 * does not own the data (samples for which the labels are to be predicted).
 * For each new data matrix, a new TSP_Predictor object must be created.
 */
class TSP_Predictor {
    private:
        TSP_Predictor() {}   // there will be no default constructor
    protected:
        FeaturePair p;
    public:
        TSP_Predictor(const FeaturePair& _p) : p(_p) {}
        TSP_Predictor(unsigned i, unsigned j, double s=0.0) :
            p(s, i, j) {}
        // Different functions to predict the label(s):
        
        // Predict all the labels for a data set.
        // X - data matrix with m rows (features) and n columns (samples)
        // y - a vector of n elements where the predicted labels will be stored
        // w - a weight for feature j
        virtual void
        predict(const double* const X, unsigned* y, unsigned m, unsigned n,
                double w=1) const;
        
        // Predict the label of one vector.
        // x - data vector of m features
        virtual unsigned predict(const double* const x, unsigned m, double w=1) const;        
};

#endif // __TSP_H__
