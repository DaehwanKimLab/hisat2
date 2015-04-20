#ifndef RLEVECTOR_H
#define RLEVECTOR_H

#include <fstream>

#include "bitvector.h"


namespace CSA
{


/*
  This class is used to construct a RLEVector.
*/

class RLEEncoder : public VectorEncoder
{
  public:
    RLEEncoder(usint block_bytes, usint superblock_size = VectorEncoder::SUPERBLOCK_SIZE);
    ~RLEEncoder();

    void setBit(usint value);
    void setRun(usint start, usint len);

    // These versions try to combine the runs if possible.
    void addBit(usint value);
    void addRun(usint start, usint len);
    void flush(); // Call this when finished.

    inline void RLEncode(usint diff, usint len)
    {
      this->size += diff + len - 1;
      this->items += len;
      this->buffer->writeDeltaCode(diff);
      this->buffer->writeDeltaCode(len);
    }

  protected:
    pair_type run;

    // These are not allowed.
    RLEEncoder();
    RLEEncoder(const RLEEncoder&);
    RLEEncoder& operator = (const RLEEncoder&);
};


/*
  This is a run-length encoded bitvector using delta coding.
*/

class RLEVector : public BitVector
{
  public:
    typedef RLEEncoder Encoder;

    explicit RLEVector(std::ifstream& file);
    explicit RLEVector(FILE* file);
    RLEVector(Encoder& encoder, usint universe_size);
    ~RLEVector();

//--------------------------------------------------------------------------

    usint reportSize() const;

//--------------------------------------------------------------------------

    class Iterator : public BitVector::Iterator
    {
      public:
        explicit Iterator(const RLEVector& par);
        ~Iterator();

        usint rank(usint value, bool at_least = false);

        usint select(usint index);
        usint selectNext();

        pair_type valueBefore(usint value);
        pair_type valueAfter(usint value);
        pair_type nextValue();

        pair_type selectRun(usint index, usint max_length);
        pair_type selectNextRun(usint max_length);

        bool isSet(usint value);

        usint countRuns();

      protected:

        void valueLoop(usint value);

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------

  protected:

    // These are not allowed.
    RLEVector();
    RLEVector(const RLEVector&);
    RLEVector& operator = (const RLEVector&);
};


} // namespace CSA


#endif // RLEVECTOR_H
