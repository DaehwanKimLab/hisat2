#ifndef LCPSAMPLES_H
#define LCPSAMPLES_H

#include <cstdio>
#include <fstream>


#ifdef SUCCINCT_LCP_VECTOR
#include "bits/succinctvector.h"
#else
#include "bits/deltavector.h"
#endif

#include "bits/array.h"
#include "misc/utils.h"


namespace CSA
{


#ifdef SUCCINCT_LCP_VECTOR
typedef SuccinctVector LCPVector;
#else
typedef DeltaVector LCPVector;
#endif


class LCPSamples
{
  public:
    #ifdef SUCCINCT_LCP_VECTOR
    const static usint INDEX_BLOCK_SIZE = 32;
    #else
    const static usint INDEX_BLOCK_SIZE = 16;
    #endif
    const static usint VALUE_BLOCK_SIZE = 32;

    explicit LCPSamples(std::ifstream& sample_file);
    explicit LCPSamples(FILE* sample_file);

    // Input pairs are of form (i, LCP[i]).
    LCPSamples(pair_type* input, usint _size, usint _items, bool report = false, bool free_input = false);

    LCPSamples(LCPVector::Encoder& index_encoder, ArrayEncoder& value_encoder, usint _size);

    ~LCPSamples();

    void writeTo(std::ofstream& sample_file) const;
    void writeTo(FILE* file) const;

    inline bool isSampled(usint sa_index) const
    {
      LCPVector::Iterator iter(*(this->indexes));
      return iter.isSet(sa_index);
    }

    // Behavior is undefined if there is no sample at index.
    inline usint getSampleAt(usint sa_index) const
    {
      LCPVector::Iterator index_iter(*(this->indexes));
      Array::Iterator value_iter(*(this->values));
      return value_iter.getItem(index_iter.rank(sa_index) - 1) - 1;
    }

    inline usint getNumberOfSamples() const { return this->items; }

    usint reportSize() const;

  private:
    usint size, items;

    LCPVector* indexes;
    Array*     values;

    // These are not allowed.
    LCPSamples();
    LCPSamples(const LCPSamples&);
    LCPSamples& operator = (const LCPSamples&);
};


} // namespace CSA


#endif // LCPSAMPLES_H
