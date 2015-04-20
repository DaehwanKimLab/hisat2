#ifndef _RLCSA_ARRAY_H
#define _RLCSA_ARRAY_H

#include <cstdio>
#include <fstream>
#include <list>

#include "../misc/definitions.h"
#include "bitbuffer.h"


namespace CSA
{

/*
  Delta-coded array for storing positive integers.
  Do not try to encode value 0!
*/


class ArrayEncoder
{
  public:
    const static usint SUPERBLOCK_SIZE = MEGABYTE;

    // We assume that superblock size is divisible by block and sample size.
    ArrayEncoder(usint block_bytes, usint superblock_size = SUPERBLOCK_SIZE);
    ~ArrayEncoder();

    void writeItem(usint value);
    inline void addItem(usint value) { this->writeItem(value); }

    usint items, blocks;
    usint block_size, superblock_bytes;

    WriteBuffer*      buffer;

    std::list<usint*> array_blocks;
    usint*            array;
    usint             blocks_in_superblock, current_blocks;

    std::list<usint*> sample_blocks;
    usint*            samples;
    usint             samples_in_superblock, current_samples;

  protected:

    // These are not allowed.
    ArrayEncoder();
    ArrayEncoder(const ArrayEncoder&);
    ArrayEncoder& operator = (const ArrayEncoder&);
};


class Array
{
  public:
    const static usint INDEX_RATE = 5;
    typedef ArrayEncoder Encoder;

    explicit Array(std::ifstream& file);
    explicit Array(FILE* file);
    explicit Array(Encoder& encoder);
    ~Array();

//--------------------------------------------------------------------------

    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;

    inline usint getSize() const { return this->items; }
    inline usint getBlockSize() const { return this->block_size; }

    usint reportSize() const;

//--------------------------------------------------------------------------

    class Iterator
    {
      public:
        Iterator(const Array& par);
        ~Iterator();

        // Returns 0 if the index is invalid.
        usint readItem(usint index);
        inline usint getItem(usint index) { return this->readItem(index); }
        usint nextItem();

        inline bool hasNext() const
        {
          return (this->sample + this->cur < this->parent.getSize());
        }

      private:
        const Array& parent;

        usint block, sample, cur, block_items;

        ReadBuffer buffer, samples;

        /*
          This function returns the sample corresponding to the block the given
          index might be found in. Parameter is assumed to be valid.
        */
        usint sampleForIndex(usint index);

        inline void getSample(usint sample_number)
        {
          this->block = sample_number;
          this->sample = this->samples.readItem(sample_number);
          this->cur = 0;
          this->block_items = this->samples.readItem() - this->sample;
          this->buffer.moveBuffer(this->parent.array + (this->block * this->parent.block_size));
        }

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------

  protected:
    usint          items;

    ReadBuffer*    buffer;
    usint*         array;
    bool           delete_array;  // Do we own the array?
    usint          block_size;
    usint          number_of_blocks;

    /*
      Each sample tells how many items are contained in the previous blocks.
    */
    ReadBuffer*    samples;
    usint          integer_bits;

    ReadBuffer*    item_index;
    usint          index_rate;

    /*
      Build a higher level index for faster queries.
      The index consists of about (number of samples) / INDEX_RATE pointers.
      The array cannot be used without the index.
    */
    void buildIndex();

    // These are not allowed.
    Array();
    Array(const Array&);
    Array& operator = (const Array&);
};


} // namespace CSA


#endif // _RLCSA_ARRAY_H
