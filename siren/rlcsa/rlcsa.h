#ifndef RLCSA_H
#define RLCSA_H

#include <fstream>
#include <vector>

#include "bits/deltavector.h"
#include "bits/rlevector.h"
#include "bits/nibblevector.h"
#include "bits/succinctvector.h"

#include "sasamples.h"
#include "alphabet.h"
#include "lcpsamples.h"
#include "misc/parameters.h"
#include "sampler.h"
#include "suffixarray.h"


namespace CSA
{


const std::string PSI_EXTENSION = ".psi";
const std::string ARRAY_EXTENSION = ".rlcsa.array";
const std::string SA_SAMPLES_EXTENSION = ".rlcsa.sa_samples";
const std::string PARAMETERS_EXTENSION = ".rlcsa.parameters";
const std::string DOCUMENT_EXTENSION = ".rlcsa.docs";
const std::string LCP_SAMPLES_EXTENSION = ".lcp_samples";
const std::string PLCP_EXTENSION = ".plcp";


const parameter_type RLCSA_BLOCK_SIZE  = parameter_type("RLCSA_BLOCK_SIZE", 32);
const parameter_type SAMPLE_RATE       = parameter_type("SAMPLE_RATE", 128);
const parameter_type SUPPORT_LOCATE    = parameter_type("SUPPORT_LOCATE", 1);
const parameter_type SUPPORT_DISPLAY   = parameter_type("SUPPORT_DISPLAY", 1);
const parameter_type WEIGHTED_SAMPLES  = parameter_type("WEIGHTED_SAMPLES", 0);


#ifdef USE_NIBBLE_VECTORS
typedef NibbleVector  PsiVector;
#else
typedef RLEVector     PsiVector;
#endif

#ifdef SUCCINCT_LCP_VECTOR
typedef SuccinctVector PLCPVector;
#else
typedef RLEVector PLCPVector;
#endif


class RLCSA
{
  friend class RLCSABuilder;

  public:

//--------------------------------------------------------------------------
//  CONSTRUCTION
//--------------------------------------------------------------------------

    const static usint ENDPOINT_BLOCK_SIZE = 16;

    explicit RLCSA(const std::string& base_name, bool print = false);

    /*
      Build RLCSA for multiple sequences, treating each \0 as an end marker.
      There must be nonzero characters between the \0s, and the last character must also be \0.
      FIXME Crashes if bytes >= 4 GB.
    */ 
    RLCSA(uchar* data, usint bytes, usint block_size, usint sa_sample_rate, usint threads, bool delete_data);

    /*
      Same as before, but this time we build the suffix array for the ranks array.
    */
    RLCSA(uchar* data, usint* ranks, usint bytes, usint block_size, usint sa_sample_rate, usint threads, bool delete_data);

    /*
      Build an RLCSA for a single sequence, and use the samples given by sampler, if not 0.
      The sequence does not contain an end marker.
      FIXME This cannot be merged. Maybe should enforce it somehow.
    */
    RLCSA(uchar* data, usint bytes, usint block_size, usint sa_sample_rate, usint threads, Sampler* sampler, bool delete_data);

    // Destroys contents of index and increment.
    RLCSA(RLCSA& index, RLCSA& increment, usint* positions, usint block_size, usint threads = 1);
    ~RLCSA();

    void writeTo(const std::string& base_name) const;

    inline bool isOk() const { return this->ok; }

//--------------------------------------------------------------------------
//  QUERIES
//--------------------------------------------------------------------------

    // These queries use SA ranges.

    // Returns the closed range containing the matches.
    pair_type count(const std::string& pattern) const;

    // Used when merging CSAs.
    void reportPositions(uchar* data, usint length, usint* positions) const;

    // Returns SA[range]. User must free the buffer. Latter version uses buffer provided by the user.
    // Direct locate means locating one position at a time.
    // Steps means that the returned values are the number of Psi steps taken, not SA values.
    usint* locate(pair_type range, bool direct = false, bool steps = false) const;
    usint* locate(pair_type range, usint* data, bool direct = false, bool steps = false) const;

    // Returns SA[index].
    usint locate(usint index, bool steps = false) const;

    // Returns T^{sequence}[range]. User must free the buffer.
    // Third version uses buffer provided by the user.
    uchar* display(usint sequence, bool include_end_marker = false) const;
    uchar* display(usint sequence, pair_type range) const;
    uchar* display(usint sequence, pair_type range, uchar* data) const;

    // Displays the intersection of T[position - context, position + len + context - 1]
    // and T^{getSequenceForPosition(position)}.
    // This is intended for displaying an occurrence of a pattern of length 'len'
    // with 'context' extra characters on both sides.
    // The actual length of the returned string is written into result_length.
    uchar* display(usint position, usint len, usint context, usint& result_length) const;

    // Returns the actual length of the prefix. User must provide the buffer.
    usint displayPrefix(usint sequence, usint len, uchar* data) const;

    // Returns at most max_len characters starting from T[SA[index]].
    // User must provide the buffer. Returns the number of characters in buffer.
    usint displayFromPosition(usint index, usint max_len, uchar* data) const;

    // Get the range of SA values for the sequence identified by
    // a sequence number or a SA value.
    pair_type getSequenceRange(usint number) const;
    pair_type getSequenceRangeForPosition(usint value) const;

    // Get the sequence number for given SA value(s).
    // The returned array is the same as the parameter.
    usint getSequenceForPosition(usint value) const;
    usint* getSequenceForPosition(usint* value, usint length) const;

    // Changes SA value to (sequence, offset).
    pair_type getRelativePosition(usint value) const;

    // Returns the BWT of the collection including end of sequence markers.
    uchar* readBWT() const;
    uchar* readBWT(pair_type range) const;

    // Returns the number of equal letter runs in the BWT. Runs consisting of end markers are ignored.
    usint countRuns() const;

    // Return the suffix array for a particular sequence. User must free the suffix array.
    SuffixArray* getSuffixArrayForSequence(usint number) const;

//--------------------------------------------------------------------------
//  SUPPORT FOR EXTERNAL MODULES: POSITIONS
//--------------------------------------------------------------------------

    // The return values of these functions are BWT indexes/ranges.

    inline usint psi(usint sa_index) const
    {
      if(sa_index >= this->data_size)
      {
        return this->data_size + this->number_of_sequences;
      }

      usint c = this->getCharacter(sa_index);
      return this->psiUnsafe(sa_index, c);
    }

    // This version returns a run.
    inline pair_type psi(usint sa_index, usint max_length) const
    {
      if(sa_index >= this->data_size)
      {
        return pair_type(this->data_size + this->number_of_sequences, 0);
      }

      usint c = this->getCharacter(sa_index);
      PsiVector::Iterator iter(*(this->array[c]));
      return iter.selectRun(sa_index - this->alphabet->cumulative(c), max_length);
    }

    inline usint LF(usint sa_index, usint c) const
    {
      if(c >= CHARS)
      {
        return this->data_size + this->number_of_sequences;
      }
      if(this->array[c] == 0)
      {
        if(c < this->alphabet->getFirstChar()) { return this->number_of_sequences - 1; }
        return this->alphabet->cumulative(c) + this->number_of_sequences - 1;
      }
      this->convertToBWTIndex(sa_index);

      PsiVector::Iterator iter(*(this->array[c]));
      return this->LF(sa_index, c, iter);
    }

    inline void convertToSAIndex(usint& bwt_index) const { bwt_index -= this->number_of_sequences; }
    inline void convertToBWTIndex(usint& sa_index) const { sa_index += this->number_of_sequences; }

    inline bool hasImplicitSample(usint bwt_index) const
    {
      return (bwt_index < this->number_of_sequences);
    }

    inline usint getImplicitSample(usint bwt_index) const
    {
      DeltaVector::Iterator iter(*(this->end_points));
      return iter.select(bwt_index) + 1;
    }

    // This is an unsafe function that returns a character value.
    inline usint getCharacter(usint sa_index) const
    {
      return this->alphabet->charAt(sa_index);
    }

//--------------------------------------------------------------------------
//  SUPPORT FOR EXTERNAL MODULES: RANGES
//--------------------------------------------------------------------------

    pair_type getSARange() const;
    pair_type getBWTRange() const;
    pair_type getCharRange(usint c) const;

    void convertToBWTRange(pair_type& sa_range) const;
    void convertToSARange(pair_type& bwt_range) const;
    void convertToSARange(std::vector<pair_type>& bwt_ranges) const;

    // This is an unsafe function that does not check its parameters.
    pair_type LF(pair_type bwt_range, usint c) const;

    // User must free the returned vector.
    std::vector<usint>* locateRange(pair_type range) const;
    std::vector<usint>* locateRanges(std::vector<pair_type>& ranges) const;
    
//--------------------------------------------------------------------------
//  REPORTING
//--------------------------------------------------------------------------

    inline bool supportsLocate() const { return this->support_locate; }
    inline bool supportsDisplay() const { return this->support_display; }
    inline usint getSize() const { return this->data_size; }
    inline usint getTextSize() const { return this->end_points->getSize(); }
    inline usint getNumberOfSequences() const { return this->number_of_sequences; }
    inline usint getBlockSize() const { return this->array[this->alphabet->getFirstChar()]->getBlockSize(); }

    // Returns the size of the data structure.
    usint reportSize(bool print = false) const;

    void printInfo() const;

//--------------------------------------------------------------------------
//  LCP EXPERIMENTS
//--------------------------------------------------------------------------

    // Optimized version:
    //   - Interleaves main loop with computing irreducible values.
    //   - Encodes maximal runs from a true local maximum to a true local minimum.
    PLCPVector* buildPLCP(usint block_size) const; //, usint threads = 1) const;

    // Returns the number of samples. sampled_values will be a pointer to the samples.
    usint sampleLCP(usint sample_rate, pair_type*& sampled_values, bool report = false) const;

    usint lcp(usint sa_index, const LCPSamples& lcp_samples, bool steps = false) const;

    // Calculate LCP[index] directly.
    usint lcpDirect(usint sa_index) const;

    // Writes PLCP[start] to PLCP[stop - 1].
    inline void encodePLCPRun(PLCPVector::Encoder& plcp, usint start, usint stop, usint first_val) const
    {
      plcp.addRun(2 * start + first_val, stop - start);
//      std::cerr << "(" << start << ", " << stop << ", " << first_val << ")" << std::endl;
    }

//--------------------------------------------------------------------------
//  INTERNAL VARIABLES
//--------------------------------------------------------------------------

  private:
    bool  ok;
    usint data_size;

    PsiVector* array[CHARS];
    Alphabet*  alphabet;
    SASamples* sa_samples;

    bool  support_locate, support_display;

    // A sequence starts at the next multiple of sample_rate after the end of previous sequence.
    usint sample_rate;
    usint number_of_sequences;
    DeltaVector* end_points;

//--------------------------------------------------------------------------
//  INTERNAL VERSIONS OF QUERIES
//--------------------------------------------------------------------------

    void  directLocate(pair_type range, usint* data, bool steps) const;
    usint directLocate(usint index, bool steps) const;
    void  locateUnsafe(pair_type range, usint* data, bool steps) const;
    bool  processRun(pair_type run, usint* data, usint* offsets, bool* finished, PsiVector::Iterator** iters, bool steps) const;
    void  displayUnsafe(pair_type range, uchar* data, bool get_ranks = false, usint* ranks = 0) const;

    void locateRange(pair_type range, std::vector<usint>& vec) const;

//--------------------------------------------------------------------------
//  INTERNAL VERSIONS OF BASIC OPERATIONS
//--------------------------------------------------------------------------

    inline usint psi(usint index, PsiVector::Iterator** iters) const
    {
      usint c = this->getCharacter(index);
      return iters[c]->select(index - this->alphabet->cumulative(c));
    }

    inline pair_type psi(usint index, usint max_length, PsiVector::Iterator** iters) const
    {
      usint c = this->getCharacter(index);
      return iters[c]->selectRun(index - this->alphabet->cumulative(c), max_length);
    }

    // Returns psi(index), assuming the suffix of rank index begins with c.
    inline usint psiUnsafe(usint index, usint c) const
    {
      PsiVector::Iterator iter(*(this->array[c]));
      return this->psiUnsafe(index, c, iter);
    }

    // As above, but with a given iterator.
    inline usint psiUnsafe(usint index, usint c, PsiVector::Iterator& iter) const
    {
      return iter.select(index - this->alphabet->cumulative(c));
    }

    // As above, but returns the next value of psi.
    inline usint psiUnsafeNext(PsiVector::Iterator& iter) const
    {
      return iter.selectNext();
    }

    inline usint LF(usint bwt_index, usint c, PsiVector::Iterator& iter) const
    {
      return this->alphabet->cumulative(c) + this->number_of_sequences + iter.rank(bwt_index) - 1;
    }

//--------------------------------------------------------------------------
//  INTERNAL STUFF
//--------------------------------------------------------------------------

    // Creates an array of iterators for every vector in this->array.
    PsiVector::Iterator** getIterators() const;
    void deleteIterators(PsiVector::Iterator** iters) const;

    void mergeEndPoints(RLCSA& index, RLCSA& increment);
    void mergeSamples(RLCSA& index, RLCSA& increment, usint* positions);

    void buildCharIndexes(usint* distribution);
    void buildRLCSA(uchar* data, usint* ranks, usint bytes, usint block_size, usint threads, Sampler* sampler, bool multiple_sequences, bool delete_data);

    // Removes structures not necessary for merging.
    void strip();

    // These are not allowed.
    RLCSA();
    RLCSA(const RLCSA&);
    RLCSA& operator = (const RLCSA&);
};


} // namespace CSA


#endif // RLCSA_H
