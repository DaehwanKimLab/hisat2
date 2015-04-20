#ifndef BITVECTOR_H
#define BITVECTOR_H

#include <cstdio>
#include <fstream>
#include <list>

#include "bitbuffer.h"


namespace CSA
{


/*
  This class provides the core functionality for encoding a bitvector.
*/

class VectorEncoder
{
  public:
    const static usint SUPERBLOCK_SIZE = MEGABYTE;

    // We assume superblock size is divisible by block and sample size.
    VectorEncoder(usint block_bytes, usint superblock_size = SUPERBLOCK_SIZE, bool _use_small_blocks = true);
    ~VectorEncoder();

/*
    This should be implemented in any inherited class.

    // These functions are assumed to be greedy, encoding the 1-bits immediately.
    void setBit(usint value);  // Values must be in increasing order.
    void setRun(usint start, usint len);

    // These versions may combine 1-bits into maximal runs.
    // Use flush() to finish the encoding when using these functions.
    // Do not mix with the greedy versions.
    void addBit(usint value);
    void addRun(usint start, usint len);
    void flush();
*/

    void addNewBlock();
    void setFirstBit(usint value);

    usint size, items, blocks;
    usint block_size, superblock_bytes;
    bool  use_small_blocks;

    WriteBuffer*      buffer;

    std::list<usint*> array_blocks;
    usint*            array;
    usint             blocks_in_superblock, current_blocks;

    std::list<usint*> sample_blocks;
    usint*            samples;
    usint             samples_in_superblock, current_samples;

  protected:

    // These are not allowed.
    VectorEncoder();
    VectorEncoder(const VectorEncoder&);
    VectorEncoder& operator = (const VectorEncoder&);
};


/*
  This class provides the core functionality for a bitvector.
  A bitvector must have at least one 1-bit.
*/

class BitVector
{
  public:
    const static usint INDEX_RATE = 5;

    explicit BitVector(std::ifstream& file);
    explicit BitVector(FILE* file);
    BitVector(VectorEncoder& encoder, usint universe_size);
    explicit BitVector(WriteBuffer& vector);
    ~BitVector();

//--------------------------------------------------------------------------

    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;

    inline usint getSize() const { return this->size; }
    inline usint getNumberOfItems() const { return this->items; }
    inline usint getBlockSize() const { return this->block_size; }

    // This returns only the sizes of the dynamically allocated structures.
    usint reportSize() const;

    usint getCompressedSize() const;

    // Removes structures not necessary for merging.
    void strip();

//--------------------------------------------------------------------------

    class Iterator
    {
      public:
        explicit Iterator(const BitVector& par);
        ~Iterator();

        inline bool hasNext() const
        {
          return (this->sample.first + this->cur < this->parent.items - 1);
        }

      protected:
        const BitVector& parent;

        usint      block;
        pair_type  sample;
        usint      cur, val, run; // cur == 0 is the sample
        usint      block_items;

        ReadBuffer buffer, samples;

        /*
          These functions return the sample corresponding to the block the given
          index/value might be found in. Parameters are assumed to be valid.
        */
        usint sampleForIndex(usint index);
        usint sampleForValue(usint value);

        inline usint getSampledIndex(usint sample_number)
        {
          return this->samples.readItem(2 * sample_number);
        }

        inline usint getSampledValue(usint sample_number)
        {
          return this->samples.readItem(2 * sample_number + 1);
        }

        inline void getSample(usint sample_number)
        {
          this->block = sample_number;
          this->samples.goToItem(2 * sample_number);
          this->sample.first = this->samples.readItem();
          this->sample.second = this->samples.readItem();
          this->cur = 0;
          this->val = this->sample.second;
          this->block_items = this->samples.readItem() - this->sample.first - 1;
          this->buffer.moveBuffer(this->parent.array + (this->block * this->parent.block_size));
        }

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

/*
    These should be implemented in any actual iterator. If the query parameter is out of bounds,
    rank-type queries return the number of items, while select-type queries return the size of
    the vector.

    // rank invalidates the "next" functionality
    // regular:   \sum_{i = 0}^{value} V[i]
    // at_least:  \sum_{i = 0}^{value - 1} V[i] + 1
    usint rank(usint value, bool at_least = false);

    usint select(usint index);      // \min value: \sum_{i = 0}^{value} V[i] = index + 1
    usint selectNext();

    // (\max i <= value: V[i] = 1, rank(i) - 1)
    // Returns (size, items) if not found.
    pair_type valueBefore(usint value);

    pair_type valueAfter(usint value); // (\min i >= value: V[i] = 1, rank(i) - 1)
    pair_type nextValue();

    // These versions of select return (value, length_of_run).
    // max_length is an upper bound for the length of the run returned.
    // V[value] is not included in the length of the run
    // These functions are not greedy: the actual length of the run can be more than reported.
    // This can happen even if max_length was not reached.
    // length_of_run is actually the number of extra items returned past value

    pair_type selectRun(usint index, usint max_length);
    pair_type selectNextRun(usint max_length);

    // isSet invalidates the "next" functionality
    bool isSet(usint value); // V[value]

    // Counts the number of 1-bit runs.
    usint countRuns();
*/

//--------------------------------------------------------------------------

  protected:
    usint size, items;

    const usint* array;
    usint        block_size;
    usint        number_of_blocks;

    /*
      Each sample is (rank(i) - 1, i) where V[i] = 1.
      Number of samples is number_of_blocks + 1.
    */
    ReadBuffer*  samples;
    usint        integer_bits;

    ReadBuffer*  rank_index;
    usint        rank_rate;

    ReadBuffer*  select_index;
    usint        select_rate;

    /*
       These functions build a higher level index for faster rank/select queries.
       The index consists of about (number of samples) / INDEX_RATE pointers.
       The bitvector cannot be used without the index.
    */
    void indexForRank();
    void indexForSelect();

    // These are used in disk storage.
    void writeHeader(std::ofstream& file) const;
    void writeHeader(FILE* file) const;
    void writeArray(std::ofstream& file) const;
    void writeArray(FILE* file) const;
    void readHeader(std::ifstream& file);
    void readHeader(FILE* file);
    void readArray(std::ifstream& file);
    void readArray(FILE* file);

    void copyArray(VectorEncoder& encoder, bool use_directly = false);

    // These are not allowed.
    BitVector();
    BitVector(const BitVector&);
    BitVector& operator = (const BitVector&);
};


} // namespace CSA


#endif // BITVECTOR_H
