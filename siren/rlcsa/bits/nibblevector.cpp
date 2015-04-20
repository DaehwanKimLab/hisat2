#include <cstdlib>

#include "nibblevector.h"
#include "../misc/utils.h"


namespace CSA
{


NibbleVector::NibbleVector(std::ifstream& file) :
  BitVector(file)
{
}

NibbleVector::NibbleVector(FILE* file) :
  BitVector(file)
{
}

NibbleVector::NibbleVector(Encoder& encoder, usint universe_size) :
  BitVector(encoder, universe_size)
{
}

NibbleVector::~NibbleVector()
{
}

//--------------------------------------------------------------------------

usint
NibbleVector::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

NibbleVector::Iterator::Iterator(const NibbleVector& par) :
  BitVector::Iterator(par),
  use_rle(false)
{
}

NibbleVector::Iterator::~Iterator()
{
}

// FIXME gap encoding for all operations

usint
NibbleVector::Iterator::rank(usint value, bool at_least)
{
  const NibbleVector& par = (const NibbleVector&)(this->parent);

  if(value >= par.size) { return par.items; }

  this->valueLoop(value);

  usint idx = this->sample.first + this->cur + 1;
  if(!at_least && this->val > value)
  {
    idx--;
  }
  if(at_least && this->val < value)
  {
    this->getSample(this->block + 1);
    idx = this->sample.first + this->cur + 1;
  }
  return idx;
}

usint
NibbleVector::Iterator::select(usint index)
{
  const NibbleVector& par = (const NibbleVector&)(this->parent);

  if(index >= par.items) { return par.size; }
  this->getSample(this->sampleForIndex(index));

  usint lim = index - this->sample.first;
  while(this->cur < lim)
  {
    this->val += this->buffer.readNibbleCode();
    usint temp = this->buffer.readNibbleCode();
    this->val += temp - 1;
    this->cur += temp;
  }
  if(this->cur > lim)
  {
    this->run = this->cur - lim;
    this->cur -= this->run;
    this->val -= this->run;
  }

  return this->val;
}

usint
NibbleVector::Iterator::selectNext()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return this->val;
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer.readNibbleCode();
    this->run = this->buffer.readNibbleCode() - 1;
  }

  return this->val;
}

pair_type
NibbleVector::Iterator::valueBefore(usint value)
{
  const NibbleVector& par = (const NibbleVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->getSample(this->sampleForValue(value));
  if(this->val > value) { return pair_type(par.size, par.items); }
  this->run = 0;

  while(this->cur < this->block_items && this->val < value)
  {
    usint temp = this->buffer.readNibbleCode(value - this->val);
    if(temp == 0) { break; }
    this->val += temp;

    temp = this->buffer.readNibbleCode();
    this->cur += temp;
    this->val += temp - 1;
  }
  if(this->val > value)
  {
    this->run = this->val - value;
    this->val = value;
    this->cur -= this->run;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
NibbleVector::Iterator::valueAfter(usint value)
{
  const NibbleVector& par = (const NibbleVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->valueLoop(value);

  if(this->val < value)
  {
    this->getSample(this->block + 1);
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
NibbleVector::Iterator::nextValue()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return pair_type(this->val, this->sample.first);
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer.readNibbleCode();
    this->run = this->buffer.readNibbleCode() - 1;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
NibbleVector::Iterator::selectRun(usint index, usint max_length)
{
  usint value = this->select(index);

  usint len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

pair_type
NibbleVector::Iterator::selectNextRun(usint max_length)
{
  usint value = this->selectNext();

  usint len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

bool
NibbleVector::Iterator::isSet(usint value)
{
  const NibbleVector& par = (const NibbleVector&)(this->parent);

  if(value >= par.size) { return false; }

  this->valueLoop(value);

  return (this->val == value);
}

usint
NibbleVector::Iterator::countRuns()
{
  const NibbleVector& par = (const NibbleVector&)(this->parent);

  if(par.items == 0) { return 0; }

  usint runs = 1;
  pair_type res = this->selectRun(0, par.items);
  usint last = res.first + res.second;

  while(last < par.size)
  {
    res = this->selectNextRun(par.items);
    if(res.first < par.size && res.first > last + 1) { runs++; }
    last = res.first + res.second;
  }

  return runs;
}

// FIXME for gap encoding
void
NibbleVector::Iterator::valueLoop(usint value)
{
  this->getSample(this->sampleForValue(value));

  if(this->val >= value) { return; }
  while(this->cur < this->block_items)
  {
    this->val += this->buffer.readNibbleCode();
    this->cur++;
    this->run = this->buffer.readNibbleCode() - 1;
    if(this->val >= value) { break; }

    this->cur += this->run;
    this->val += this->run;
    if(this->val >= value)
    {
      this->run = this->val - value;
      this->val = value;
      this->cur -= this->run;
      break;
    }
    this->run = 0;
  }
}

//--------------------------------------------------------------------------

NibbleEncoder::NibbleEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size),
  run(EMPTY_PAIR)
{
}

NibbleEncoder::~NibbleEncoder()
{
}

void
NibbleEncoder::setBit(usint value)
{
  this->setRun(value, 1);
}

// FIXME for gap encoding
void
NibbleEncoder::setRun(usint start, usint len)
{
  if(this->items == 0)
  {
    this->setFirstBit(start);
    if(len > 1)
    {
      this->nibbleEncode(1, len - 1);
    }
    return;
  }
  if(start < this->size || len == 0) { return; }

  // Write as much into the buffer as possible.
  usint diff = start + 1 - this->size;
  usint free_bits = this->buffer->bitsLeft();
  usint code_bits = this->buffer->nibbleCodeLength(diff);
  if(free_bits >= code_bits + 4) // At least a part of the run fits into the block.
  {
    free_bits -= code_bits;
    usint run_bits = this->buffer->nibbleCodeLength(len);
    if(run_bits <= free_bits)
    {
      this->nibbleEncode(diff, len);
      return;
    }

    // Encode as much as possible and let the rest spill.
    usint llen = (usint)1 << (3 * (free_bits / 4));
    this->nibbleEncode(diff, llen);
    len -= llen;

    // A new sample will be added.
    this->size++;
    this->items++;
  }
  else
  {
    this->size = start + 1;
    this->items++;
  }

  // Didn't fit into the block. A new sample & block required.
  this->addNewBlock();
  if(len > 1)
  {
    this->nibbleEncode(1, len - 1);
  }
}

void
NibbleEncoder::addBit(usint value)
{
  this->addRun(value, 1);
}

void
NibbleEncoder::addRun(usint start, usint len)
{
  if(this->run.second == 0)
  {
    this->run = pair_type(start, len);
  }
  else if(start == this->run.first + this->run.second)
  {
    this->run.second += len;
  }
  else
  {
    this->setRun(this->run.first, this->run.second);
    this->run = pair_type(start, len);
  }
}

void
NibbleEncoder::flush()
{
  this->setRun(this->run.first, this->run.second);
  this->run.second = 0;
}


} // namespace CSA
