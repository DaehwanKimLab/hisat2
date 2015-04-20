#ifndef ADAPTIVE_SAMPLES_H
#define ADAPTIVE_SAMPLES_H


#include <fstream>

#include "rlcsa.h"
#include "misc/parameters.h"


namespace CSA
{

//--------------------------------------------------------------------------

const parameter_type CANDIDATE_SAMPLES   = parameter_type("CANDIDATE_SAMPLES", 0);
const parameter_type HALF_GREEDY_SAMPLES = parameter_type("HALF_GREEDY_SAMPLES", 0);
const parameter_type SAMPLE_PROMOTE_RATE = parameter_type("SAMPLE_PROMOTE_RATE", 0);
const parameter_type SAMPLE_WINDOW_SIZE  = parameter_type("SAMPLE_WINDOW_SIZE", 65536);

//--------------------------------------------------------------------------

/*
  This uses both POSITIONS and RANGES parts of external module interface.
  FIXME Currently supports only one sequence.
  FIXME The implementation currently assumes regular sampling with Sampler as input.
  FIXME Assumes that key_type is 32 bits.
*/

class AdaptiveSamples
{
  public:
    typedef uint key_type;
    typedef std::pair<key_type, key_type> sample_type;

    AdaptiveSamples(const RLCSA& rlcsa, const std::string& base_name);
    ~AdaptiveSamples();

//--------------------------------------------------------------------------

    inline bool isOk() const { return this->ok; }
    inline bool supportsLocate() const { return true; }
    inline bool supportsDisplay() const { return this->half_greedy; }

    usint getNumberOfSamples() const;
    double getLoad() const;
    double getPrimaryLoad() const;
    double getCandidateLoad() const;

    usint reportSize() const;
    void report() const;

//--------------------------------------------------------------------------

    usint locate(usint i, bool steps);
    usint* locate(pair_type range, bool steps);
    void display(pair_type range, uchar* data);

//--------------------------------------------------------------------------

  private:
    const RLCSA& index;

    sample_type* samples;
    usint        size, items;

    sample_type* candidate_samples;
    usint        c_items;

    SASamples*   regular_samples;

    key_type text_start;

    bool      use_candidates, half_greedy;
    double    promote_probability;
    usint     window_size, pos, sum;
    key_type* window;

    bool ok;
    
    void addSample(sample_type sample);
    void addCandidate(sample_type sample);
    void trySample(sample_type sample, usint offsets);

    key_type hash(key_type key);
    key_type secondaryHash(key_type key);

    // These are not allowed.
    AdaptiveSamples();
    AdaptiveSamples(const AdaptiveSamples&);
    AdaptiveSamples& operator = (const AdaptiveSamples&);
};

//--------------------------------------------------------------------------

};  // namespace CSA


#endif // ADAPTIVE_SAMPLES_H
