#include <cstring>
#include <cstdlib>

#include "bitvector.h"


namespace CSA
{


BitVector::BitVector(std::ifstream& file) :
  rank_index(0), select_index(0)
{
  this->readHeader(file);
  this->readArray(file);

  this->integer_bits = length(this->size);
  this->samples = new ReadBuffer(file, 2 * (this->number_of_blocks + 1), this->integer_bits);

  this->indexForRank();
  this->indexForSelect();
}

BitVector::BitVector(FILE* file) :
  rank_index(0), select_index(0)
{
  this->readHeader(file);
  this->readArray(file);

  this->integer_bits = length(this->size);
  this->samples = new ReadBuffer(file, 2 * (this->number_of_blocks + 1), this->integer_bits);

  this->indexForRank();
  this->indexForSelect();
}

BitVector::BitVector(VectorEncoder& encoder, usint universe_size) :
  size(universe_size), items(encoder.items),
  block_size(encoder.block_size),
  number_of_blocks(encoder.blocks),
  rank_index(0), select_index(0)
{
  if(this->items == 0)
  {
    std::cerr << "BitVector: Cannot create a bitvector with no 1-bits!" << std::endl;
    return;
  }
  this->copyArray(encoder);

  this->integer_bits = length(this->size);
  WriteBuffer sample_buffer(2 * (this->number_of_blocks + 1), this->integer_bits);

  // Compress the samples.
  for(std::list<usint*>::iterator iter = encoder.sample_blocks.begin(); iter != encoder.sample_blocks.end(); iter++)
  {
    usint* buf = *iter;
    for(usint i = 0; i < 2 * encoder.samples_in_superblock; i++)
    {
      sample_buffer.writeItem(buf[i]);
    }
  }
  for(usint i = 0; i < 2 * encoder.current_samples; i++)
  {
    sample_buffer.writeItem(encoder.samples[i]);
  }
  sample_buffer.writeItem(this->items);
  sample_buffer.writeItem(this->size);

  this->samples = sample_buffer.getReadBuffer();

  this->indexForRank();
  this->indexForSelect();
}

BitVector::BitVector() :
  array(0), samples(0), rank_index(0), select_index(0)
{
}

BitVector::~BitVector()
{
  delete[] this->array;
  delete this->samples;
  delete this->rank_index;
  delete this->select_index;
}

//--------------------------------------------------------------------------

void
BitVector::writeTo(std::ofstream& file) const
{
  this->writeHeader(file);
  this->writeArray(file);
  this->samples->writeBuffer(file);
}

void
BitVector::writeTo(FILE* file) const
{
  this->writeHeader(file);
  this->writeArray(file);
  this->samples->writeBuffer(file);
}

void
BitVector::writeHeader(std::ofstream& file) const
{
  file.write((char*)&(this->size), sizeof(this->size));
  file.write((char*)&(this->items), sizeof(this->items));
  file.write((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.write((char*)&(this->block_size), sizeof(this->block_size));
}

void
BitVector::writeHeader(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(&(this->size), sizeof(this->size), 1, file);
  std::fwrite(&(this->items), sizeof(this->items), 1, file);
  std::fwrite(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file);
  std::fwrite(&(this->block_size), sizeof(this->block_size), 1, file);
}

void
BitVector::writeArray(std::ofstream& file) const
{
  file.write((char*)(this->array), this->block_size * this->number_of_blocks * sizeof(usint));
}

void
BitVector::writeArray(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(this->array, this->block_size * sizeof(usint), this->number_of_blocks, file);
}

void
BitVector::readHeader(std::ifstream& file)
{
  file.read((char*)&(this->size), sizeof(this->size));
  file.read((char*)&(this->items), sizeof(this->items));
  file.read((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.read((char*)&(this->block_size), sizeof(this->block_size));
}

void
BitVector::readHeader(FILE* file)
{
  if(file == 0) { return; }
  if(!std::fread(&(this->size), sizeof(this->size), 1, file)) { return; }
  if(!std::fread(&(this->items), sizeof(this->items), 1, file)) { return; }
  if(!std::fread(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file)) { return; }
  if(!std::fread(&(this->block_size), sizeof(this->block_size), 1, file)) { return; }
}

void
BitVector::readArray(std::ifstream& file)
{
  usint* array_buffer = new usint[this->block_size * this->number_of_blocks];
  file.read((char*)(array_buffer), this->block_size * this->number_of_blocks * sizeof(usint));
  this->array = array_buffer;
}

void
BitVector::readArray(FILE* file)
{
  if(file == 0) { return; }
  usint* array_buffer = new usint[this->block_size * this->number_of_blocks];
  if(!std::fread(array_buffer, this->block_size * sizeof(usint), this->number_of_blocks, file)) { return; }
  this->array = array_buffer;
}

//--------------------------------------------------------------------------

void
BitVector::copyArray(VectorEncoder& encoder, bool use_directly)
{
  if(use_directly && encoder.array_blocks.size() == 0)
  {
    this->array = encoder.array; encoder.array = 0;
    return;
  }

  usint* array_buffer = new usint[this->block_size * this->number_of_blocks];
  usint pos = 0, total_words = this->block_size * this->number_of_blocks;

  for(std::list<usint*>::iterator iter = encoder.array_blocks.begin(); iter != encoder.array_blocks.end(); iter++)
  {
    memcpy(array_buffer + pos, *iter, encoder.superblock_bytes);
    pos += encoder.superblock_bytes / sizeof(usint);
  }

  usint remainder = std::min(total_words - pos, (usint)(encoder.superblock_bytes / sizeof(usint)));
  memcpy(array_buffer + pos, encoder.array, remainder * sizeof(usint));
  pos += remainder;
  if(pos < total_words) { memset(array_buffer + pos, 0, (total_words - pos) * sizeof(usint)); }

  this->array = array_buffer;
}

//--------------------------------------------------------------------------

usint
BitVector::reportSize() const
{
  // We assume the reportSize() of derived classes includes any class variables of BitVector.
  usint bytes = this->block_size * this->number_of_blocks * sizeof(usint);
  if(this->samples != 0) { bytes += this->samples->reportSize(); }
  if(this->rank_index != 0) { bytes += this->rank_index->reportSize(); }
  if(this->select_index != 0) { bytes += this->select_index->reportSize(); }
  return bytes;
}

usint
BitVector::getCompressedSize() const
{
  return this->block_size * this->number_of_blocks * sizeof(usint);
}

//--------------------------------------------------------------------------

void
BitVector::strip()
{
  delete this->rank_index; this->rank_index = 0;
}

//--------------------------------------------------------------------------

void
BitVector::indexForRank()
{
  delete this->rank_index;

  usint value_samples = (this->number_of_blocks + BitVector::INDEX_RATE - 1) / BitVector::INDEX_RATE;
  this->rank_rate = (this->size + value_samples - 1) / value_samples;
  value_samples = (this->size + this->rank_rate - 1) / this->rank_rate + 1;
  WriteBuffer index_buffer(value_samples, length(this->number_of_blocks - 1));

  // current is value, pointer is sample number.
  usint current = 0, pointer = 0;
  this->samples->goToItem(2);
  while(this->samples->hasNextItem())
  {
    this->samples->skipItem();
    usint limit = this->samples->readItem();  // Next sampled value.
    while(current < limit)
    {
      index_buffer.writeItem(pointer);
      current += this->rank_rate;
    }
    pointer++;
  }
  index_buffer.writeItem(this->number_of_blocks - 1);

  this->rank_index = index_buffer.getReadBuffer();
}

void
BitVector::indexForSelect()
{
  delete this->select_index;

  usint index_samples = (this->number_of_blocks + BitVector::INDEX_RATE - 1) / BitVector::INDEX_RATE;
  this->select_rate = (this->items + index_samples - 1) / index_samples;
  index_samples = (this->items + this->select_rate - 1) / this->select_rate + 1;
  WriteBuffer index_buffer(index_samples, length(this->number_of_blocks - 1));

  // current is index, pointer is sample number.
  usint current = 0, pointer = 0;
  this->samples->goToItem(2);
  while(this->samples->hasNextItem())
  {
    usint limit = this->samples->readItem();  // Next sampled index.
    this->samples->skipItem();
    while(current < limit)
    {
      index_buffer.writeItem(pointer);
      current += this->select_rate;
    }
    pointer++;
  }
  index_buffer.writeItem(this->number_of_blocks - 1);

  this->select_index = index_buffer.getReadBuffer();
}

//--------------------------------------------------------------------------

BitVector::Iterator::Iterator(const BitVector& par) :
  parent(par),
  buffer(par.array, par.block_size),
  samples(*(par.samples))
{
}

BitVector::Iterator::~Iterator()
{
}

usint
BitVector::Iterator::sampleForIndex(usint index)
{
  usint low = this->parent.select_index->readItemConst(index / this->parent.select_rate);
  usint high = this->parent.number_of_blocks - 1;

  this->samples.goToItem(2 * low + 2);
  for(; low < high; low++)
  {
    if(this->samples.readItem() > index) { return low; }
    this->samples.skipItem();
  }

  return low;
}

usint
BitVector::Iterator::sampleForValue(usint value)
{
  usint low = this->parent.rank_index->readItemConst(value / this->parent.rank_rate);
  usint high = this->parent.number_of_blocks - 1;

  this->samples.goToItem(2 * low + 3);
  for(; low < high; low++)
  {
    if(this->samples.readItem() > value) { return low; }
    this->samples.skipItem();
  }

  return low;
}

//--------------------------------------------------------------------------

VectorEncoder::VectorEncoder(usint block_bytes, usint superblock_size, bool _use_small_blocks) :
  size(0), items(0), blocks(0),
  block_size(BYTES_TO_WORDS(block_bytes)),
  superblock_bytes(superblock_size),
  use_small_blocks(_use_small_blocks),
  blocks_in_superblock(1), current_blocks(0),
  samples(0), samples_in_superblock(0), current_samples(0)
{
  this->array = new usint[this->superblock_bytes / sizeof(usint)];
  memset(this->array, 0, this->superblock_bytes);
  if(use_small_blocks)
  {
    this->blocks_in_superblock = this->superblock_bytes / (sizeof(usint) * this->block_size);
    this->buffer = new WriteBuffer(this->array, this->block_size);
    this->samples = new usint[this->superblock_bytes / sizeof(usint)];
    this->samples_in_superblock = this->superblock_bytes / (2 * sizeof(usint));
  }
  else
  {
    this->buffer = new WriteBuffer(this->array, BYTES_TO_WORDS(this->superblock_bytes));
    this->blocks = this->current_blocks = 1;
  }
}

VectorEncoder::~VectorEncoder()
{
  delete[] this->array;

  delete this->buffer;
  for(std::list<usint*>::iterator iter = this->array_blocks.begin(); iter != this->array_blocks.end(); iter++)
  {
    delete[] *iter;
  }

  delete[] this->samples;
  for(std::list<usint*>::iterator iter = this->sample_blocks.begin(); iter != this->sample_blocks.end(); iter++)
  {
    delete[] *iter;
  }
}

void
VectorEncoder::addNewBlock()
{
  // Do we need a new superblock for the block?
  this->blocks++; this->current_blocks++;
  if(this->current_blocks > this->blocks_in_superblock)
  {
    this->array_blocks.push_back(this->array);
    this->array = new usint[this->superblock_bytes / sizeof(usint)];
    memset(this->array, 0, this->superblock_bytes);
    this->current_blocks = 1;
  }
  this->buffer->moveBuffer(this->array + (this->block_size * (this->current_blocks - 1)));

  // Do we need a new superblock for the sample?
  if(this->use_small_blocks)
  {
    this->current_samples++;
    if(this->current_samples > this->samples_in_superblock)
    {
      this->sample_blocks.push_back(this->samples);
      this->samples = new usint[this->superblock_bytes / sizeof(usint)];
      this->current_samples = 1;
    }
    this->samples[2 * this->current_samples - 2] = this->items - 1;
    this->samples[2 * this->current_samples - 1] = this->size - 1;
  }
}

void
VectorEncoder::setFirstBit(usint value)
{
  this->samples[0] = 0;
  this->samples[1] = value;

  this->size = value + 1;
  this->items = 1;
  this->blocks = 1;

  this->current_blocks = 1;
  this->current_samples = 1;
}


} // namespace CSA
