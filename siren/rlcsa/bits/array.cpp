#include <cstring>
#include <cstdlib>

#include "array.h"


namespace CSA
{


Array::Array(std::ifstream& file) :
  delete_array(true), item_index(0)
{
  file.read((char*)&(this->items), sizeof(this->items));
  file.read((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.read((char*)&(this->block_size), sizeof(this->block_size));

  this->array = new usint[this->block_size * this->number_of_blocks];
  file.read((char*)(this->array), this->block_size * this->number_of_blocks * sizeof(usint));
  this->buffer = new ReadBuffer(this->array, this->block_size);

  this->integer_bits = length(this->items);
  this->samples = new ReadBuffer(file, this->number_of_blocks + 1, this->integer_bits);

  this->buildIndex();
}

Array::Array(FILE* file) :
  delete_array(true), item_index(0)
{
  if(!std::fread(&(this->items), sizeof(this->items), 1, file)) { return; }
  if(!std::fread(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file)) { return; }
  if(!std::fread(&(this->block_size), sizeof(this->block_size), 1, file)) { return; }

  this->array = new usint[this->block_size * this->number_of_blocks];
  if(!std::fread(this->array, this->block_size * sizeof(usint), this->number_of_blocks, file)) { return; }
  this->buffer = new ReadBuffer(this->array, this->block_size);

  this->integer_bits = length(this->items);
  this->samples = new ReadBuffer(file, this->number_of_blocks + 1, this->integer_bits);

  this->buildIndex();
}

Array::Array(Encoder& encoder) :
  items(encoder.items),
  delete_array(true),
  block_size(encoder.block_size),
  number_of_blocks(encoder.blocks),
  item_index(0)
{
  usint* array_buffer = new usint[this->block_size * this->number_of_blocks];
  this->array = array_buffer;
  this->buffer = new ReadBuffer(this->array, this->block_size);

  this->integer_bits = length(this->items);
  WriteBuffer sample_buffer(this->number_of_blocks + 1, this->integer_bits);

  // Copy & linearize the array.
  usint pos = 0;
  for(std::list<usint*>::iterator iter = encoder.array_blocks.begin(); iter != encoder.array_blocks.end(); iter++)
  {
    memcpy(array_buffer + pos, *iter, encoder.superblock_bytes);
    pos += encoder.block_size * encoder.blocks_in_superblock;
  }
  memcpy(array_buffer + pos, encoder.array, encoder.current_blocks * encoder.block_size * sizeof(usint));

  // Compress the samples.
  for(std::list<usint*>::iterator iter = encoder.sample_blocks.begin(); iter != encoder.sample_blocks.end(); iter++)
  {
    usint* buf = *iter;
    for(usint i = 0; i < encoder.samples_in_superblock; i++)
    {
      sample_buffer.writeItem(buf[i]);
    }
  }
  for(usint i = 0; i < encoder.current_samples; i++)
  {
    sample_buffer.writeItem(encoder.samples[i]);
  }
  sample_buffer.writeItem(this->items);

  this->samples = sample_buffer.getReadBuffer();

  this->buildIndex();
}

Array::~Array()
{
  if(this->delete_array) { delete[] this->array; }
  delete this->buffer;
  delete this->samples;
  delete this->item_index;
}

//--------------------------------------------------------------------------

void
Array::writeTo(std::ofstream& file) const
{
  file.write((char*)&(this->items), sizeof(this->items));
  file.write((char*)&(this->number_of_blocks), sizeof(this->number_of_blocks));
  file.write((char*)&(this->block_size), sizeof(this->block_size));
  file.write((char*)(this->array), this->block_size * this->number_of_blocks * sizeof(usint));
  this->samples->writeBuffer(file);
}

void
Array::writeTo(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(&(this->items), sizeof(this->items), 1, file);
  std::fwrite(&(this->number_of_blocks), sizeof(this->number_of_blocks), 1, file);
  std::fwrite(&(this->block_size), sizeof(this->block_size), 1, file);
  std::fwrite(this->array, this->block_size * sizeof(usint), this->number_of_blocks, file);
  this->samples->writeTo(file);
}

usint
Array::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += this->buffer->reportSize();
  bytes += this->block_size * this->number_of_blocks * sizeof(usint);
  bytes += this->samples->reportSize();
  bytes += this->item_index->reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

void
Array::buildIndex()
{
  delete this->item_index;

  usint index_samples = (this->number_of_blocks + Array::INDEX_RATE - 1) / Array::INDEX_RATE;
  this->index_rate = (this->items + index_samples - 1) / index_samples;
  index_samples = (this->items + this->index_rate - 1) / this->index_rate + 1;
  WriteBuffer index_buffer(index_samples, length(this->number_of_blocks - 1));

  usint current = 0, pointer = 0;
  this->samples->goToItem(1);
  while(this->samples->hasNextItem())
  {
    usint limit = this->samples->readItem();
    while(current < limit)
    {
      index_buffer.writeItem(pointer);
      current += this->index_rate;
    }
    pointer++;
  }
  index_buffer.writeItem(this->number_of_blocks - 1);

  this->item_index = index_buffer.getReadBuffer();
}

//--------------------------------------------------------------------------

Array::Iterator::Iterator(const Array& par) :
  parent(par),
  buffer(*(par.buffer)),
  samples(*(par.samples))
{
}

Array::Iterator::~Iterator()
{
}

usint
Array::Iterator::readItem(usint index)
{
  if(index >= this->parent.items) { return 0; }
  this->getSample(this->sampleForIndex(index));

  usint value = 0;
  while(this->sample + this->cur <= index)
  {
    value = this->buffer.readDeltaCode();
    this->cur++;
  }

  return value;
}

usint
Array::Iterator::nextItem()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
  }

  this->cur++;
  return this->buffer.readDeltaCode();
}

usint
Array::Iterator::sampleForIndex(usint index)
{
  usint low = this->parent.item_index->readItemConst(index / this->parent.index_rate);
  usint high = this->parent.number_of_blocks - 1;

  this->samples.goToItem(low + 1);
  for(; low < high; low++)
  {
    if(this->samples.readItem() > index) { return low; }
  }

  return low;
}

//--------------------------------------------------------------------------

ArrayEncoder::ArrayEncoder(usint block_bytes, usint superblock_size) :
  items(0), blocks(1),
  block_size(BYTES_TO_WORDS(block_bytes)),
  superblock_bytes(superblock_size)
{
  this->array = new usint[this->superblock_bytes / sizeof(usint)];
  memset(this->array, 0, this->superblock_bytes);
  this->blocks_in_superblock = this->superblock_bytes / (sizeof(usint) * this->block_size);
  this->current_blocks = 1;

  this->samples = new usint[this->superblock_bytes / sizeof(usint)];
  this->samples_in_superblock = this->superblock_bytes / sizeof(usint);
  this->current_samples = 1;
  this->samples[0] = 0;

  this->buffer = new WriteBuffer(this->array, this->block_size);
}

ArrayEncoder::~ArrayEncoder()
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
ArrayEncoder::writeItem(usint value)
{
  this->items++;
  if(this->buffer->writeDeltaCode(value)) { return; }

  // Didn't fit into the block. A new block & sample required.
  this->blocks++;
  this->current_blocks++;
  this->current_samples++;

  // Do we need a new superblock for the block?
  if(this->current_blocks > this->blocks_in_superblock)
  {
    this->array_blocks.push_back(this->array);
    this->array = new usint[this->superblock_bytes / sizeof(usint)];
    memset(this->array, 0, this->superblock_bytes);
    this->current_blocks = 1;
  }
  this->buffer->moveBuffer(this->array + (this->block_size * (this->current_blocks - 1)));

  // Do we need a new superblock for the sample?
  if(this->current_samples > this->samples_in_superblock)
  {
    this->sample_blocks.push_back(this->samples);
    this->samples = new usint[this->superblock_bytes / sizeof(usint)];
    this->current_samples = 1;
  }
  this->samples[this->current_samples - 1] = this->items - 1;

  // Finally, write the item.
  this->buffer->writeDeltaCode(value);
}


} // namespace CSA
