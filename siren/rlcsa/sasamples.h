#ifndef SASAMPLES_H
#define SASAMPLES_H

#include <cstdio>
#include <fstream>

#include "sampler.h"
#include "misc/utils.h"
#include "bits/bitbuffer.h"

#include "bits/deltavector.h"

#ifdef SUCCINCT_SA_VECTOR
#include "bits/succinctvector.h"
#endif

namespace CSA
{


#ifdef SUCCINCT_SA_VECTOR
typedef SuccinctVector SAVector;
#else
typedef DeltaVector SAVector;
#endif


class SASamples
{
  public:
    #ifdef SUCCINCT_SA_VECTOR
    const static usint INDEX_BLOCK_SIZE = 32;
    #else
    const static usint INDEX_BLOCK_SIZE = 16;
    #endif

    SASamples(std::ifstream& sample_file, usint sample_rate, bool _weighted);
    SASamples(FILE* sample_file, usint sample_rate, bool _weighted);

    // These assume < 4 GB data.
    SASamples(short_pair* sa, DeltaVector* end_points, usint data_size, usint sample_rate, usint threads);
    SASamples(short_pair* sa, Sampler* sampler, usint threads); // Use the given samples.

    // Use these samples. Assumes regular sampling.
    SASamples(pair_type* sample_pairs, usint data_size, usint sample_rate, usint threads);

    ~SASamples();

    // Destroys contents of index and increment.
    // We assume index and increment have same sample rate.
    // positions must not containt the positions of end of sequence markers.
    // number_of_sequences is subtracted from each position before the value is used.
    SASamples(SASamples& index, SASamples& increment, usint* positions, usint number_of_positions, usint number_of_sequences);

    void writeTo(std::ofstream& sample_file) const;
    void writeTo(FILE* sample_file) const;

    // Returns (i, inverseSA(i)) such that i is the last sampled position up to value.
    // Value is actual 0-based suffix array value.
    // Returns (size, size) if value is too large.
    pair_type inverseSA(usint value) const;

    // Returns the value of ith sample in suffix array order.
    inline usint getSample(usint i) const
    {
      return std::min(this->samples->readItemConst(i) * this->rate, this->size - 1);
    }

    // Returns (ind, sample number) where ind >= index or (size, ???).
    inline pair_type getFirstSampleAfter(usint index) const
    {
      SAVector::Iterator iter(*(this->indexes));
      return iter.valueAfter(index);
    }

    inline bool isSampled(usint index) const
    {
      SAVector::Iterator iter(*(this->indexes));
      return iter.isSet(index);
    }

    inline usint getSampleAt(usint index) const
    {
      SAVector::Iterator iter(*(this->indexes));
      return this->getSample(iter.rank(index) - 1);
    }

    inline usint getSampleRate() const { return this->rate; }
    inline usint getNumberOfSamples() const { return this->items; }

    inline bool isWeighted() const { return this->weighted; }
    inline bool supportsLocate() const { return (this->samples != 0); }
    inline bool supportsDisplay() const { return (this->inverse_samples != 0); }

    usint reportSize() const;

    // Removes structures not necessary for merging.
    void strip();

  private:
    bool weighted;
    usint rate, size, items;

    SAVector*   indexes;
    ReadBuffer* samples;

    SAVector*   inverse_indexes;
    ReadBuffer* inverse_samples;

    void buildInverseSamples();

    // Weighted case.
    void buildSamples(pair_type* sample_pairs, bool inverse, usint threads);

    // Note: contents of original samples are deleted.
    void mergeSamples(SASamples& index, SASamples& increment, usint* positions, usint n, usint skip);

    // These are not allowed.
    SASamples();
    SASamples(const SASamples&);
    SASamples& operator = (const SASamples&);
};


} // namespace CSA


#endif // SASAMPLES_H
