#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include "misc/definitions.h"


namespace CSA
{


const std::string SA_EXTENSION = ".sa";


class SuffixArray
{
  public:

//--------------------------------------------------------------------------
//  CONSTRUCTION
//--------------------------------------------------------------------------

    // The input data is expected to contain end markers, but they are not included in
    // the suffix array. Unlike in RLCSA, there is no padding between the sequences.
    // In the last two constructors, the ownership of data is transfered to the suffix array.
    // The third version assumes that the rank array has unique lexicographic ranks for all suffixes.
    SuffixArray(const std::string& base_name, bool print = false);
    SuffixArray(uchar* _data, uint bytes, uint threads = 1);
    SuffixArray(uchar* _data, usint* ra, uint bytes, uint threads = 1);
    ~SuffixArray();

    void writeTo(const std::string& base_name, bool write_data = false) const;

    inline bool isOk() const { return this->ok; }

//--------------------------------------------------------------------------
//  QUERIES
//--------------------------------------------------------------------------

    // Returns the closed range containing the matches.
    pair_type count(const std::string& pattern) const;

    // User must free the array.
    uint* locate(pair_type range) const;

    // Returns SA[index].
    uint locate(uint index) const;

//--------------------------------------------------------------------------
//  LCP
//--------------------------------------------------------------------------

    // This uses the Irreducible LCP algorithm. User must free the array.
    // Note that the array includes positions corresponding to end markers.
    uint* getLCPArray(bool getPLCP) const;

//--------------------------------------------------------------------------
//  REPORTING
//--------------------------------------------------------------------------

    inline uint getSize() const { return this->data_size - this->sequences; }
    inline uint getNumberOfSequences() const { return this->sequences; }

    inline const usint* getRanks() const { return this->ranks; }

    // Returns the size of the data structure.
    usint reportSize(bool print = false) const;

//--------------------------------------------------------------------------
//  INTERNAL VARIABLES
//--------------------------------------------------------------------------

  private:
    bool ok;

    uchar* data;
    uint*  sa;
    uint*  original_sa;
    usint* ranks; // Lexicographic ranks relative to a larger suffix array.
    uint   data_size;
    uint   sequences;

//--------------------------------------------------------------------------
//  INTERNAL STUFF
//--------------------------------------------------------------------------

    uint match(const std::string& pattern, uint index) const;

    // These are not allowed.
    SuffixArray();
    SuffixArray(const SuffixArray&);
    SuffixArray& operator = (const SuffixArray&);
};


} // namespace CSA


#endif // SUFFIXARRAY_H
