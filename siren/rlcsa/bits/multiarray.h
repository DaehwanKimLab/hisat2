#ifndef MULTIARRAY_H
#define MULTIARRAY_H

#include <fstream>
#include <list>

#include "array.h"
#include "bitbuffer.h"
#include "succinctvector.h"


namespace CSA
{

/*
  A single array split into multiple logical arrays. Can store either fixed-width non-negative integers
  or Delta-coded positive integers.
*/


class MultiArray
{
  public:
    static MultiArray* createFixed(usint items, usint bits);
    static MultiArray* createDelta(usint block_size);
    static MultiArray* readFrom(std::ifstream& file);
    static MultiArray* readFrom(FILE* file);

    virtual ~MultiArray();

    virtual bool writeTo(std::ofstream& file) const;
    virtual bool writeTo(FILE* file) const;
    virtual usint reportSize() const; // Does not include sizeof(*this).

//--------------------------------------------------------------------------

    inline bool isWritable() const { return (this->status == writable); }
    inline bool isOk() const { return (this->status == ok); }

    inline usint getSize() const { return this->items; }
    inline usint getNumberOfArrays() const { return this->arrays; }

    enum encoding_type { fixed_width, delta_coded };
    inline encoding_type getEncoding() const { return this->encoding; }

//--------------------------------------------------------------------------

    virtual bool isFull() const = 0;
    virtual void writeItem(usint item) = 0; // No safety checks!

    bool nextArray(); // Returns false if already at the beginning of an array.
    bool finishWriting(); // Fails if we are at the beginning of an empty array.

//--------------------------------------------------------------------------

    class Iterator
    {
      public:
        virtual ~Iterator();

        virtual void goToItem(usint array, usint index) = 0;
        virtual usint readItem(usint array, usint index) = 0;
        virtual usint nextItem() = 0;
        inline bool hasNext() const { return (this->pos < this->parent.items); }

        // End should be one position past the last item we are interested in.
        void setEnd(usint array, usint index);
        inline bool atEnd() const { return (this->pos >= this->limit); }

      protected:
        const MultiArray& parent;

        SuccinctVector::Iterator array_iter;

        usint pos, limit;

        explicit Iterator(const MultiArray& par);

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

    // User must delete the iterator. This is ugly.
    virtual Iterator* getIterator() const = 0;

//--------------------------------------------------------------------------

  protected:
    enum status_type { error, writable, ok };

    const static usint DELTA_FLAG = 0x01;
    const static usint ITEM_BITS_MASK = 0x00FF0000;
    const static usint ITEM_BITS_SHIFT = 16;

    const static usint VECTOR_BLOCK_SIZE = 32;

    status_type   status;
    encoding_type encoding;
    usint         items, arrays, current_array_begin;

    SuccinctVector::Encoder* array_encoder;
    SuccinctVector*          array_borders;

    MultiArray(std::ifstream& file, usint flags);
    MultiArray(FILE* file, usint flags);
    MultiArray(usint _items, usint _bits);  // Fixed width.
    explicit MultiArray(usint block_size);  // Delta coded.

    virtual usint getFlags() const = 0;

    virtual void encodeItems() = 0;

    // These are not allowed.
    MultiArray();
    MultiArray(const MultiArray&);
    MultiArray& operator = (const MultiArray&);
};

//--------------------------------------------------------------------------

class FixedMultiArray : public MultiArray
{
  public:
    FixedMultiArray(std::ifstream& file, usint flags);
    FixedMultiArray(FILE* file, usint flags);
    FixedMultiArray(usint _items, usint bits);
    virtual ~FixedMultiArray();

    virtual bool writeTo(std::ofstream& file) const;
    virtual bool writeTo(FILE* file) const;
    virtual usint reportSize() const;

    virtual bool isFull() const;
    virtual void writeItem(usint item);

    class Iterator : public MultiArray::Iterator
    {
      public:
        explicit Iterator(const FixedMultiArray& par);
        virtual ~Iterator();

        virtual void goToItem(usint array, usint index);
        virtual usint readItem(usint array, usint index);
        virtual usint nextItem();

      protected:
        ReadBuffer fixed_iter;

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

    virtual MultiArray::Iterator* getIterator() const { return new Iterator(*this); }

  protected:
    virtual usint getFlags() const;
    virtual void encodeItems();

  private:
    WriteBuffer* fixed_buffer;
    ReadBuffer*  fixed_items;

    // These are not allowed.
    FixedMultiArray();
    FixedMultiArray(const FixedMultiArray&);
    FixedMultiArray& operator = (const FixedMultiArray&);
};

//--------------------------------------------------------------------------

class DeltaMultiArray : public MultiArray
{
  public:
    DeltaMultiArray(std::ifstream& file, usint flags);
    DeltaMultiArray(FILE* file, usint flags);
    explicit DeltaMultiArray(usint block_size);
    virtual ~DeltaMultiArray();

    virtual bool writeTo(std::ofstream& file) const;
    virtual bool writeTo(FILE* file) const;
    virtual usint reportSize() const;

    virtual bool isFull() const;
    virtual void writeItem(usint item);

    class Iterator : public MultiArray::Iterator
    {
      public:
        explicit Iterator(const DeltaMultiArray& par);
        virtual ~Iterator();

        virtual void goToItem(usint array, usint index);
        virtual usint readItem(usint array, usint index);
        virtual usint nextItem();

      protected:
        Array::Iterator delta_iter;

        pair_type buffer;

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

    virtual MultiArray::Iterator* getIterator() const { return new Iterator(*this); }

  protected:
    virtual usint getFlags() const;
    virtual void encodeItems();

  private:
    Array::Encoder* delta_buffer;
    Array*          delta_items;

    // These are not allowed.
    DeltaMultiArray();
    DeltaMultiArray(const DeltaMultiArray&);
    DeltaMultiArray& operator = (const DeltaMultiArray&);
};

} // namespace CSA


#endif // MULTIARRAY_H
