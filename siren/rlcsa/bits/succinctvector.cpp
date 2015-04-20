#include <cstring>

#include "succinctvector.h"
#include "../misc/utils.h"


namespace CSA
{


SuccinctVector::SuccinctVector(std::ifstream& file) :
  BitVector()
{
  this->readHeader(file);
  this->readArray(file);

  this->integer_bits = length(this->size);
  this->indexForRank();
  this->indexForSelect();
}

SuccinctVector::SuccinctVector(FILE* file) :
  BitVector()
{
  if(file == 0) { return; }

  this->readHeader(file);
  this->readArray(file);

  this->integer_bits = length(this->size);
  this->indexForRank();
  this->indexForSelect();
}

SuccinctVector::SuccinctVector(Encoder& encoder, usint universe_size) :
  BitVector()
{
  if(encoder.items == 0)
  {
    std::cerr << "SuccinctVector: Cannot create a bitvector with no 1-bits!" << std::endl;
    return;
  }

  this->size = universe_size;
  this->items = encoder.items;
  this->block_size = encoder.block_size;

  this->number_of_blocks = (this->size + this->block_size * WORD_BITS - 1) / (this->block_size * WORD_BITS);
  this->copyArray(encoder);

  this->integer_bits = length(this->size);
  this->indexForRank();
  this->indexForSelect();
}

SuccinctVector::SuccinctVector(ReadBuffer& buffer, usint block_bytes, usint universe_size) :
  BitVector()
{
  this->initializeUsing(buffer.rawData(), block_bytes, universe_size);
}

SuccinctVector::SuccinctVector(WriteBuffer& buffer, usint block_bytes, usint universe_size) :
  BitVector()
{
  if(buffer.ownsData())
  {
    this->initializeFrom(buffer.stealData(), block_bytes, universe_size);
  }
  else
  {
    this->initializeUsing(buffer.rawData(), block_bytes, universe_size);
  }
}

SuccinctVector::~SuccinctVector()
{
}

//--------------------------------------------------------------------------

void
SuccinctVector::writeTo(std::ofstream& file) const
{
  this->writeHeader(file);
  this->writeArray(file);
}

void
SuccinctVector::writeTo(FILE* file) const
{
  this->writeHeader(file);
  this->writeArray(file);
}

usint
SuccinctVector::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

void
SuccinctVector::initializeUsing(const usint* buffer, usint block_bytes, usint universe_size)
{
  this->size = universe_size;
  this->block_size = BYTES_TO_WORDS(block_bytes);
  this->number_of_blocks = (this->size + this->block_size * WORD_BITS - 1) / (this->block_size * WORD_BITS);
  this->integer_bits = length(this->size);

  usint* buf = new usint[this->block_size * this->number_of_blocks];
  memcpy(buf, buffer, sizeof(usint) * this->block_size * this->number_of_blocks);
  this->array = buf;

  this->indexForRank();
  if(this->items == 0)
  {
    std::cerr << "SuccinctVector: Cannot create a bitvector with no 1-bits!" << std::endl;
    return;
  }
  this->indexForSelect();
}

void
SuccinctVector::initializeFrom(usint* buffer, usint block_bytes, usint universe_size)
{
  this->size = universe_size;
  this->block_size = BYTES_TO_WORDS(block_bytes);
  this->number_of_blocks = (this->size + this->block_size * WORD_BITS - 1) / (this->block_size * WORD_BITS);
  this->array = buffer;
  this->integer_bits = length(this->size);

  this->indexForRank();
  if(this->items == 0)
  {
    std::cerr << "SuccinctVector: Cannot create a bitvector with no 1-bits!" << std::endl;
    return;
  }
  this->indexForSelect();
}

void
SuccinctVector::indexForRank()
{
  delete this->rank_index; this->rank_index = 0;

  WriteBuffer buffer(this->number_of_blocks + 1, this->integer_bits);
  const usint* data = this->array;

  buffer.writeItem(0);
  usint bitcount = 0;
  for(usint block = 0; block < this->number_of_blocks; block++)
  {
    for(usint word = 0; word < this->block_size; word++, ++data) { bitcount += popcount(*data); }
    buffer.writeItem(bitcount);
  }

  this->items = bitcount;
  this->rank_index = buffer.getReadBuffer();
}

void
SuccinctVector::indexForSelect()
{
  delete this->select_index; this->select_index = 0;

  this->select_rate = (this->items + this->number_of_blocks - 1) / this->number_of_blocks;
  usint index_samples = (this->items + this->select_rate - 1) / this->select_rate + 1;
  WriteBuffer buffer(index_samples, this->integer_bits);

  // pointer is the largest block number such that rank_index[pointer] <= current.
  usint pointer = 0;
  usint next_rank = this->rank_index->readItem(1);
  for(usint current = 0; current < this->items; current += this->select_rate)
  {
    while(next_rank <= current)
    {
      pointer++; next_rank = this->rank_index->readItem();
    }
    buffer.writeItem(pointer);
  }
  buffer.writeItem(this->number_of_blocks - 1);

  this->select_index = buffer.getReadBuffer();
}

//--------------------------------------------------------------------------

SuccinctVector::Iterator::Iterator(const SuccinctVector& par) :
  parent(par),
  cur(0)
{
}

SuccinctVector::Iterator::~Iterator()
{
}

usint
SuccinctVector::Iterator::rank(usint value, bool at_least)
{
  const SuccinctVector& par = (const SuccinctVector&)(this->parent);

  if(value >= par.size) { return par.items; }

  usint block = value / (par.block_size * WORD_BITS);
  usint result = par.rank_index->readItemConst(block);
  const usint* data = par.array + block * par.block_size;

  // Value is now the number of bits to be read within the block.
  // Add popcounts for entire words as long as possible.
  value = value - block * par.block_size * WORD_BITS + 1;
  while(value >= WORD_BITS)
  {
    result += popcount(*data);
    ++data; value -= WORD_BITS;
  }

  // Add popcount for the last bits.
  if(value > 0)
  {
    result += popcount(GET(LOWER(*data, WORD_BITS - value), value));
  }

  return result;
}

usint
SuccinctVector::Iterator::select(usint index)
{
  const SuccinctVector& par = (const SuccinctVector&)(this->parent);
  ReadBuffer rank_index(*(par.rank_index));

  if(index >= par.items) { return par.size; }

  // The range of blocks that contains the 1-bit of rank (index + 1).
  usint low = par.select_index->readItemConst(index / par.select_rate);
  usint high = par.select_index->readItemConst(index / par.select_rate + 1);

  // Use binary search if the range is too long.
  while(high - low >= SuccinctVector::SHORT_RANGE)
  {
    usint mid = low + (high - low) / 2;
    usint val = rank_index.readItem(mid);
    if(val > index) { high = mid - 1; }
    else            { low = mid; }
  }

  // Find the correct block by linear search.
  for(usint limit = rank_index.readItem(low + 1);
      limit <= index;
      limit = rank_index.readItem(), low++);

  // Finally scan the block to find the 1-bit of rank (index + 1).
  // result is the number of bits scanned.
  // cur is the rank of the last 1-bit in the scanned range.
  const usint* data = par.array + low * par.block_size;
  usint result = low * par.block_size * WORD_BITS;
  this->cur = rank_index.readItem(low);
  for(usint count = popcount(*data); this->cur + count <= index; count = popcount(*data))
  {
    ++data; result += WORD_BITS; this->cur += count;
  }
  for(usint bit = (usint)1 << (WORD_BITS - 1); this->cur <= index; bit >>= 1)
  {
    result++; if(*data & bit) { this->cur++; }
  }

  this->cur--;
  return result - 1;
}

usint
SuccinctVector::Iterator::selectNext()
{
  return this->select(this->cur + 1);
}

pair_type
SuccinctVector::Iterator::selectRun(usint index, usint max_length)
{
  return pair_type(this->select(index), 0);
}

pair_type
SuccinctVector::Iterator::selectNextRun(usint max_length)
{
  return pair_type(this->select(this->cur + 1), 0);
}

pair_type
SuccinctVector::Iterator::valueBefore(usint value)
{
  const SuccinctVector& par = (const SuccinctVector&)(this->parent);
  if(value >= par.size) { return pair_type(par.size, par.items); }

	usint temp = this->rank(value);
	if(temp == 0)         { return pair_type(par.size, par.items); }
	return pair_type(this->select(temp - 1), temp - 1);
}

pair_type
SuccinctVector::Iterator::valueAfter(usint value)
{
  const SuccinctVector& par = (const SuccinctVector&)(this->parent);
  if(value >= par.size) { return pair_type(par.size, par.items); }

	if(value == 0) { return pair_type(this->select(0), 0); }
	usint temp = this->rank(value - 1);
	return pair_type(this->select(temp), temp);
}

pair_type
SuccinctVector::Iterator::nextValue()
{
	usint temp = this->cur + 1;
	return pair_type(this->select(temp), temp);
}

bool
SuccinctVector::Iterator::isSet(usint value)
{
  const SuccinctVector& par = (const SuccinctVector&)(this->parent);
  if(value >= par.size) { return false; }
  return (par.array[value / WORD_BITS] & ((usint)1 << (WORD_BITS - value % WORD_BITS - 1)));
}

usint
SuccinctVector::Iterator::countRuns()
{
  usint runs = 0;
  usint prev = this->select(0);
  while(this->hasNext())
  {
    usint temp = this->selectNext();
    if(temp != prev + 1) { runs++; }
    prev = temp;
  }
  return runs;
}

//--------------------------------------------------------------------------

SuccinctEncoder::SuccinctEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size, false)
{
}

SuccinctEncoder::~SuccinctEncoder()
{
}

void
SuccinctEncoder::setBit(usint value)
{
  // We can set any bit in the current superblock as well as those bits after it.
  if(value < this->array_blocks.size() * this->superblock_bytes * CHAR_BIT) { return; }
  this->size = std::max(this->size, value + 1);

  value -= this->array_blocks.size() * this->superblock_bytes * CHAR_BIT;
  while(value >= this->superblock_bytes * CHAR_BIT)
  {
    this->addNewBlock();
    value -= this->superblock_bytes * CHAR_BIT;
  }

  if(!(this->buffer->isSet(value))) { this->items++; }
  this->buffer->setBit(value);
}

void
SuccinctEncoder::setRun(usint start, usint len)
{
  for(usint i = start; i < start + len; i++) { this->setBit(i); }
}

void
SuccinctEncoder::addBit(usint value)
{
  this->setBit(value);
}

void
SuccinctEncoder::addRun(usint start, usint len)
{
  this->setRun(start, len);
}

void
SuccinctEncoder::flush()
{
}

} // namespace CSA
