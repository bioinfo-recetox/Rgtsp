#ifndef __FPAIR_H__
#define __FPAIR_H__

/*
 * Title      : btypes.h
 * Author     : Vlad Popovici
 * Description: Definition of some basic types.
 */

#include <set>

typedef long long unsigned LU_INT;  // widest unsigned type available

/**
 * FeaturePair: encodes a pair of features and the corresponding score.
 * The class provides methods for ordering pairs of features. All members are
 * public, is the user's responsability to ensure consistency of the objects.
 */
class FeaturePair
{
public:
    double score;
    unsigned i, j;
	
    FeaturePair() : score(-1.0), i(0), j(0) {}
    
    FeaturePair(const FeaturePair& p) : score(p.score), i(p.i), j(p.j) {}
	
    FeaturePair(double _score, unsigned _i, unsigned _j) :
        score(_score), i(_i), j(_j) {}
        
    bool operator== (const FeaturePair& other) const
    { return (score == other.score); }
};

/**
 * CompareFeaturePairs_inc: a simple class for ordering increasingly
 * lists of FeaturePair objects.
 */
class CompareFeaturePairs_inc          // order increasingly the list
{
public:
    bool operator() (const FeaturePair& one, const FeaturePair& two) const
    { return (one.score < two.score); }
};

/**
 * CompareFeaturePairs_dec: a simple class for ordering decreasingly
 * lists of FeaturePair objects.
 */
class CompareFeaturePairs_dec          // order decreasingly the list
{
public:
    bool operator() (const FeaturePair& one, const FeaturePair& two) const
    { return (one.score > two.score); }
};

/**
 * FeatureCollection_inc: an increasingly ordered list of feature pairs.
 */
typedef std::multiset<FeaturePair, CompareFeaturePairs_inc> FeatureCollection_inc;

/**
 * FeatureCollection_dec: an decreasingly ordered list of feature pairs.
 */
typedef std::multiset<FeaturePair, CompareFeaturePairs_dec> FeatureCollection_dec;

#endif // __FPAIR_H__
