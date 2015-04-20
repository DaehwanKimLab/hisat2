#include <cstdlib>

#include "deltavector.h"


namespace CSA
{


DeltaVector::DeltaVector(std::ifstream& file) :
  BitVector(file)
{
}

DeltaVector::DeltaVector(FILE* file) :
  BitVector(file)
{
}

DeltaVector::DeltaVector(Encoder& encoder, usint universe_size) :
  BitVector(encoder, universe_size)
{
}

DeltaVector::~DeltaVector()
{
}

//--------------------------------------------------------------------------

usint
DeltaVector::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

DeltaVector::Iterator::Iterator(const DeltaVector& par) :
  BitVector::Iterator(par)
{
}

DeltaVector::Iterator::~Iterator()
{
}

usint
DeltaVector::Iterator::rank(usint value, bool at_least)
{
  const DeltaVector& par = (const DeltaVector&)(this->parent);

  if(value >= par.size) { return par.items; }
  this->getSample(this->sampleForValue(value));

  while(this->cur < this->block_items && this->val < value)
  {
    this->val += this->buffer.readDeltaCode();
    this->cur++;
  }

  usint idx = this->sample.first + this->cur + 1;
  if(!at_least && this->val > value) { idx--; }
  if(at_least && this->val < value)
  {
    this->getSample(this->block + 1);
    idx = this->sample.first + this->cur + 1;
  }
  return idx;
}

usint
DeltaVector::Iterator::select(usint index)
{
  const DeltaVector& par = (const DeltaVector&)(this->parent);

  if(index >= par.items) { return par.size; }
  this->getSample(this->sampleForIndex(index));

  usint lim = index - this->sample.first;
  for(; this->cur < lim; this->cur++)
  {
    this->val += this->buffer.readDeltaCode();
  }

  return this->val;
}

usint
DeltaVector::Iterator::selectNext()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return this->val;
  }

  this->cur++;
  this->val += this->buffer.readDeltaCode();
  return this->val;
}

pair_type
DeltaVector::Iterator::valueBefore(usint value)
{
  const DeltaVector& par = (const DeltaVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->getSample(this->sampleForValue(value));
  if(this->val > value) { return pair_type(par.size, par.items); }

  while(this->cur < this->block_items && this->val < value)
  {
    usint temp = this->buffer.readDeltaCode(value - this->val);
    if(temp == 0) { break; }
    this->val += temp;
    this->cur++;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
DeltaVector::Iterator::valueAfter(usint value)
{
  const DeltaVector& par = (const DeltaVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }
  this->getSample(this->sampleForValue(value));

  while(this->cur < this->block_items && this->val < value)
  {
    this->val += this->buffer.readDeltaCode();
    this->cur++;
  }
  if(this->val < value)
  {
    this->getSample(this->block + 1);
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
DeltaVector::Iterator::nextValue()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    return pair_type(this->val, this->sample.first);
  }

  this->cur++;
  this->val += this->buffer.readDeltaCode();
  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
DeltaVector::Iterator::selectRun(usint index, usint max_length)
{
  return pair_type(this->select(index), 0);
}

pair_type
DeltaVector::Iterator::selectNextRun(usint max_length)
{
  return pair_type(this->selectNext(), 0);
}

bool
DeltaVector::Iterator::isSet(usint value)
{
  const DeltaVector& par = (const DeltaVector&)(this->parent);

  if(value >= par.size) { return false; }
  this->getSample(this->sampleForValue(value));

  while(this->cur < this->block_items && this->val < value)
  {
    this->val += this->buffer.readDeltaCode();
    this->cur++;
  }

  return (this->val == value);
}

//--------------------------------------------------------------------------

DeltaEncoder::DeltaEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size)
{
}

DeltaEncoder::~DeltaEncoder()
{
}

void
DeltaEncoder::setBit(usint value)
{
  if(this->items == 0)
  {
    this->setFirstBit(value);
    return;
  }
  if(value < this->size) { return; }

  usint diff = value + 1 - this->size;
  this->size = value + 1;
  this->items++;
  if(this->buffer->writeDeltaCode(diff)) { return; }

  // Didn't fit into the block. A new sample & block required.
  this->addNewBlock();
}

void
DeltaEncoder::setRun(usint start, usint len)
{
  for(usint i = start; i < start + len; i++) { this->setBit(i); }
}

void
DeltaEncoder::addBit(usint value)
{
  this->setBit(value);
}

void
DeltaEncoder::addRun(usint start, usint len)
{
  this->setRun(start, len);
}

void
DeltaEncoder::flush()
{
}

} // namespace CSA
