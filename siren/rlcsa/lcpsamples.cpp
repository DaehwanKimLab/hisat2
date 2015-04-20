#include <algorithm>
#include <fstream>
#include <iostream>

#include "lcpsamples.h"
#include "misc/utils.h"


namespace CSA
{


LCPSamples::LCPSamples(std::ifstream& sample_file)
{
  this->indexes = new LCPVector(sample_file);
  this->values = new Array(sample_file);

  this->size = indexes->getSize();
  this->items = indexes->getNumberOfItems();
}

LCPSamples::LCPSamples(FILE* sample_file)
{
  if(sample_file == 0) { return; }

  this->indexes = new LCPVector(sample_file);
  this->values = new Array(sample_file);

  this->size = indexes->getSize();
  this->items = indexes->getNumberOfItems();
}

LCPSamples::LCPSamples(pair_type* input, usint _size, usint _items, bool report, bool free_input) :
  size(_size),
  items(_items)
{
  LCPVector::Encoder index_encoder(LCPSamples::INDEX_BLOCK_SIZE);
  ArrayEncoder value_encoder(LCPSamples::VALUE_BLOCK_SIZE);

  usint max_sample = 0;
  for(usint i = 0; i < this->items; i++)
  {
    index_encoder.setBit(input[i].first);
    value_encoder.addItem(input[i].second + 1);
    max_sample = std::max(max_sample, input[i].second);
  }

  this->indexes = new LCPVector(index_encoder, this->size);
  this->values = new Array(value_encoder);

  if(free_input) { delete[] input; }

  if(report)
  {
    usint max_bits = length(max_sample);
    double total_size = max_bits * (double)(this->items) / (CHAR_BIT * (double)MEGABYTE);
    std::cout << "Maximum sample value: " << max_sample << " (" << max_bits << " bits)" << std::endl;
    std::cout << "Raw samples: " << total_size << " MB" << std::endl;
    std::cout << "Encoded samples: " << (this->values->reportSize() / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << std::endl;
  }
}

LCPSamples::LCPSamples(LCPVector::Encoder& index_encoder, ArrayEncoder& value_encoder, usint _size) :
  size(_size)
{
  this->indexes = new LCPVector(index_encoder, this->size);
  this->values = new Array(value_encoder);
  this->items = this->indexes->getNumberOfItems();
}

LCPSamples::~LCPSamples()
{
  delete this->indexes;
  delete this->values;
}

void
LCPSamples::writeTo(std::ofstream& sample_file) const
{
  this->indexes->writeTo(sample_file);
  this->values->writeTo(sample_file);
}

void
LCPSamples::writeTo(FILE* sample_file) const
{
  this->indexes->writeTo(sample_file);
  this->values->writeTo(sample_file);
}

//--------------------------------------------------------------------------

usint
LCPSamples::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += this->indexes->reportSize();
  bytes += this->values->reportSize();
  return bytes;
}

//--------------------------------------------------------------------------


} // namespace CSA
