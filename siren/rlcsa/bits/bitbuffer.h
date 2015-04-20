#ifndef BITBUFFER_H
#define BITBUFFER_H

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "../misc/definitions.h"


namespace CSA
{


class ReadBuffer
{
  public:
    ReadBuffer(std::ifstream& file, usint words);
    ReadBuffer(FILE* file, usint words);
    ReadBuffer(std::ifstream& file, usint _items, usint item_size);
    ReadBuffer(FILE* file, usint _items, usint item_size);

    // These versions do not delete the data when deleted.
    ReadBuffer(const usint* buffer, usint words);
    ReadBuffer(const usint* buffer, usint _items, usint item_size);
    ReadBuffer(const ReadBuffer& original);

    ~ReadBuffer();

//--------------------------------------------------------------------------

    void claimData(); // Make the buffer own the current data.
    const usint* rawData(); // Get access to the raw data.

    void writeTo(std::ofstream& file) const;
    inline void writeBuffer(std::ofstream& file) const { this->writeTo(file); }
    void writeTo(FILE* file) const;
    inline void writeBuffer(FILE* file) const { this->writeTo(file); }

    // Write the buffer in an exportable format.
    // If the buffer contains 1-bit items, each bit becomes '0' or '1'.
    // For other item sizes, each item becomes a word.
    void exportTo(const std::string& filename) const;

    // The buffer will no longer own the data.
    void moveBuffer(const usint* buffer);

    usint reportSize() const;

    inline usint getSize() const { return this->size; } // In words.
    inline usint getItemSize() const { return this->item_bits; }
    inline usint getNumberOfItems() const { return this->items; }
    inline bool ownsData() const { return this->free_buffer; }

//--------------------------------------------------------------------------

    inline void reset()
    {
      this->pos = 0;
      this->bits = WORD_BITS;
      this->current = 0;
    }

    inline void skipBits(usint count)
    {
      if(count < this->bits)
      {
        this->bits -= count;
        return;
      }

      count -= this->bits;
      this->pos += 1 + count / WORD_BITS;
      this->bits = WORD_BITS - count % WORD_BITS;
    }

//--------------------------------------------------------------------------

    inline usint bitsLeft() const
    {
      return this->bits + WORD_BITS * (this->size - this->pos - 1);
    }

    // Returns nonzero if bit is 1
    inline usint isSet(usint index) const
    {
      return this->data[index / WORD_BITS] & ((usint)1 << (WORD_BITS - index % WORD_BITS - 1));
    }

    // Returns nonzero if bit is 1
    inline usint readBit()
    {
      this->bits--;
      usint bit = this->data[this->pos] & ((usint)1 << this->bits);

      if(this->bits == 0) { this->pos++; this->bits = WORD_BITS; }

      return bit;
    }

    inline usint readBits(usint count)
    {
      usint value = 0;

      if(count >= this->bits)
      {
        count -= this->bits;
        value |= HIGHER(GET(this->data[this->pos], this->bits), count);
        this->pos++; this->bits = WORD_BITS;
      }
      if(count > 0)
      {
        this->bits -= count;
        value |= GET(LOWER(this->data[this->pos], this->bits), count);
      }

      return value;
    }

//--------------------------------------------------------------------------

    /*
      These operations work on 4-bit nibbles.
      Do not use these with the bit-level operations.
    */

    inline usint readNibble()
    {
      this->bits -= 4;
      usint value = GET(LOWER(this->data[this->pos], this->bits), 4);

      if(this->bits == 0) { this->pos++; this->bits = WORD_BITS; }

      return value;
    }

    // Nibble code for positive integers.
    inline usint readNibbleCode()
    {
      usint temp, value = 0, shift = 0;
      do
      {
        temp = this->readNibble();
        value |= (temp & 0x7) << shift;
        shift += 3;
      }
      while((temp & 0x8) == 0);

      return value + 1;
    }

    // This version reads the code only if value <= limit.
    inline usint readNibbleCode(usint limit)
    {
       usint _pos = this->pos, _bits = this->bits;
       usint value = this->readNibbleCode();
       if(value > limit)
       {
         this->pos = _pos; this->bits = _bits;
         return 0;
       }
       return value;
    }

//--------------------------------------------------------------------------

    /*
      These operations work on fixed-size items. No sanity checks are made
      for parameter values.
    */

    inline void goToItem(usint item)
    {
      usint b = item * this->item_bits;
      this->pos = b / WORD_BITS;
      this->bits = WORD_BITS - b % WORD_BITS;
      this->current = item;
    }

    inline usint nextItem()
    {
      this->current++;
      return this->readBits(this->item_bits);
    }

    inline usint readItem() { return this->nextItem(); }

    inline usint readItem(usint item)
    {
      this->goToItem(item);
      return this->nextItem();
    }

    inline usint readFirstItem()
    {
      return this->readItem(0);
    }

    inline usint readItemConst(usint item) const
    {
      usint b = item * this->item_bits;
      usint p = b / WORD_BITS;
      b = WORD_BITS - b % WORD_BITS;

      usint c = this->item_bits;
      usint value = 0;

      if(c >= b)
      {
        c -= b;
        value |= HIGHER(GET(this->data[p], b), c);
        p++; b = WORD_BITS;
      }
      if(c > 0)
      {
        b -= c;
        value |= GET(LOWER(this->data[p], b), c);
      }

      return value;
    }

    inline bool hasNextItem() const
    {
      return (this->current < this->items);
    }

    inline void skipItem()
    {
      this->skipBits(this->item_bits);
      this->current++;
    }

//--------------------------------------------------------------------------

    /*
      Delta coding for positive integers
    */

    inline usint readDeltaCode()
    {
      usint len = 0;
      while(this->readBit() == 0) { len++; }

      usint temp = (((usint)1 << len) | this->readBits(len)) - 1;
      temp = ((usint)1 << temp) | this->readBits(temp);
      return temp;
    }

    // This version reads the code only if value <= limit.
    inline usint readDeltaCode(usint limit)
    {
       usint _pos = this->pos, _bits = this->bits;
       usint value = this->readDeltaCode();
       if(value > limit)
       {
         this->pos = _pos; this->bits = _bits;
         return 0;
       }
       return value;
    }

//--------------------------------------------------------------------------

  private:
    const usint* data;
    usint size, item_bits, items;
    bool  free_buffer;

    // Iterator data
    usint pos, bits, current;

    inline static usint bitsToWords(usint _bits) { return (_bits + WORD_BITS - 1) / WORD_BITS; }

    // These are not allowed.
    ReadBuffer();
    ReadBuffer& operator = (const ReadBuffer&);
};


//--------------------------------------------------------------------------


class WriteBuffer
{
  public:
    explicit WriteBuffer(usint words);
    WriteBuffer(usint _items, usint item_size);

    // These versions do not delete the data when deleted.
    WriteBuffer(usint* buffer, usint words);
    WriteBuffer(usint* buffer, usint _items, usint item_size);

    ~WriteBuffer();

//--------------------------------------------------------------------------

    // Get access to the raw data.
    usint* rawData();
    usint* stealData(); // The buffer no longer owns the data.

    // These transfer the ownership of the data to the ReadBuffer.
    // Resizing does not affect the allocated size.
    // The last version does not check whether the items fit in the buffer.
    ReadBuffer* getReadBuffer();
    ReadBuffer* getReadBuffer(usint words);
    ReadBuffer* getReadBuffer(usint _items, usint item_size);

    void writeTo(std::ofstream& file) const;
    inline void writeBuffer(std::ofstream& file) const { this->writeTo(file); }
    void writeTo(FILE* file) const;
    inline void writeBuffer(FILE* file) const { this->writeTo(file); }

    // The buffer will no longer own the data.
    void moveBuffer(usint* buffer);

    usint reportSize() const;

    inline usint getSize() const { return this->size; } // In words.
    inline usint getItemSize() const { return this->item_bits; }
    inline usint getNumberOfItems() const { return this->items; }
    inline bool ownsData() const { return this->free_buffer; }

//--------------------------------------------------------------------------

    inline void reset()
    {
      this->pos = 0;
      this->bits = WORD_BITS;
      this->current = 0;
    }

    inline void skipBits(usint count)
    {
      if(count < this->bits)
      {
        this->bits -= count;
        return;
      }

      count -= this->bits;
      this->pos += 1 + count / WORD_BITS;
      this->bits = WORD_BITS - count % WORD_BITS;
    }

//--------------------------------------------------------------------------

    inline usint bitsLeft() const
    {
      return this->bits + WORD_BITS * (this->size - this->pos - 1);
    }

    inline void writeBits(usint value, usint count)
    {
      if(count >= this->bits)
      {
        count -= this->bits;
        this->data[this->pos] &= HIGHER(WORD_MAX, this->bits);  // Clear low-order bits.
        this->data[this->pos] |= GET(LOWER(value, count), this->bits);
        this->pos++; this->bits = WORD_BITS;
      }
      if(count > 0)
      {
        this->bits -= count;
        // Clear low-order bits.
        this->data[this->pos] &= ~HIGHER(LOWER(WORD_MAX, WORD_BITS - count), this->bits);
        this->data[this->pos] |= HIGHER(GET(value, count), this->bits);
      }
    }

//--------------------------------------------------------------------------

    // Returns nonzero if bit is 1
    inline usint isSet(usint index) const
    {
      return this->data[index / WORD_BITS] & ((usint)1 << (WORD_BITS - index % WORD_BITS - 1));
    }

    inline void setBit(usint index)
    {
      this->data[index / WORD_BITS] |= (usint)1 << (WORD_BITS - index % WORD_BITS - 1);
    }

    inline void unsetBit(usint index)
    {
      this->data[index / WORD_BITS] &= ~((usint)1 << (WORD_BITS - index % WORD_BITS - 1));
    }

//--------------------------------------------------------------------------

    /*
      These operations work on fixed-size items.
    */

    inline void goToItem(usint item)
    {
      usint b = item * this->item_bits;
      this->pos = b / WORD_BITS;
      this->bits = WORD_BITS - b % WORD_BITS;
      this->current = item;
    }

    inline bool hasNextItem() const
    {
      return (this->current < this->items);
    }

    inline void writeItem(usint item)
    {
      this->writeBits(item, this->item_bits);
      this->current++;
    }

    inline void skipItem()
    {
      this->skipBits(this->item_bits);
      this->current++;
    }

//--------------------------------------------------------------------------

    /*
      Nibble coding for positive integers.
    */

    inline usint nibbleCodeLength(usint value) const
    {
      usint b = 0;
      value--;

      do
      {
        b += 4;
        value >>= 3;
      }
      while(value > 0);

      return b;
    }

    // Something breaks very badly if value > 15.
    inline void writeNibble(usint value)
    {
      this->bits -= 4;
      this->data[this->pos] |= HIGHER(value, this->bits);
      if(this->bits == 0) { this->pos++; this->bits = WORD_BITS; }
    }

    // It is assumed that there is enough space for the code.
    inline void writeNibbleCode(usint value)
    {
      value--;
      while(value > 0x7)
      {
        this->writeNibble(value & 0x7);
        value >>= 3;
      }
      this->writeNibble(value | 0x8);
    }

//--------------------------------------------------------------------------

    /*
      Delta coding for positive integers
    */

    inline bool canDeltaCode(usint value) const
    {
      return this->deltaCodeLength(value) <= this->bitsLeft();
    }

    inline usint deltaCodeLength(usint value) const
    {
      usint len = length(value);
      usint llen = length(len);
      return (len + llen + llen - 2);
    }

    // This version returns false if there is no space left for the encoding.
    inline bool writeDeltaCode(usint value)
    {
      usint len = length(value);
      usint llen = length(len);

      if(len + llen + llen - 2 > this->bitsLeft()) { return false; }

      // this->writeBits(0, llen - 1); // Now included in the next writeBits()
      this->writeBits(len, llen + llen - 1);
      this->writeBits(value, len - 1);
      return true;
    }

    // This version assumes the code fits into the buffer.
    inline void writeDeltaCodeDirect(usint value)
    {
      usint len = length(value);
      usint llen = length(len);

      // this->writeBits(0, llen - 1); // Now included in the next writeBits()
      this->writeBits(len, llen + llen - 1);
      this->writeBits(value, len - 1);
    }

    // We assume the code fits into usint:
    //  32-bit:  value < 2^24
    //  64-bit:  value < 2^54
    inline void writeDeltaCodeFast(usint value)
    {
      usint len = length(value);

      value ^= ((usint)1 << (len - 1));
      this->writeBits((len << (len - 1)) | value, len + 2 * length(len) - 2);
    }

//--------------------------------------------------------------------------

  private:
    usint* data;
    usint size, item_bits, items;
    bool free_buffer;

    // Iterator data
    usint pos, bits, current;

    inline static usint bitsToWords(usint _bits) { return (_bits + WORD_BITS - 1) / WORD_BITS; }

    // These are not allowed.
    WriteBuffer();
    WriteBuffer(const WriteBuffer&);
    WriteBuffer& operator = (const WriteBuffer&);
};


//--------------------------------------------------------------------------


} // namespace CSA


#endif // BITBUFFER_H
