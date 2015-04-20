#include <cstdlib>
#include <cstring>

#include "bitbuffer.h"


namespace CSA
{

//--------------------------------------------------------------------------

ReadBuffer::ReadBuffer(std::ifstream& file, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(true)
{
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  file.read((char*)buffer, this->size * sizeof(usint));
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(FILE* file, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(true)
{
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  if(file != 0)
  {
    if(!std::fread(buffer, this->size * sizeof(usint), 1, file)) { return; }
  }
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(std::ifstream& file, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(true)
{
  this->size = bitsToWords(this->items * this->item_bits);
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  file.read((char*)buffer, this->size * sizeof(usint));
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(FILE* file, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(true)
{
  this->size = bitsToWords(this->items * this->item_bits);
  usint* buffer = new usint[this->size];
  memset(buffer, 0, this->size * sizeof(usint));
  if(file != 0)
  {
    if(!std::fread(buffer, this->size * sizeof(usint), 1, file)) { return; }
  }
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(const usint* buffer, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(false)
{
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(const usint* buffer, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(false)
{
  this->size = bitsToWords(this->items * this->item_bits);
  this->data = buffer;
  this->reset();
}

ReadBuffer::ReadBuffer(const ReadBuffer& original) :
  data(original.data),
  size(original.size),
  item_bits(original.item_bits),
  items(original.items),
  free_buffer(false)
{
  this->reset();
}

ReadBuffer::~ReadBuffer()
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
}

//--------------------------------------------------------------------------

void
ReadBuffer::claimData()
{
  this->free_buffer = true;
}

const usint*
ReadBuffer::rawData()
{
  return this->data;
}

void
ReadBuffer::writeTo(std::ofstream& file) const
{
  file.write((const char*)this->data, this->size * sizeof(usint));
}

void
ReadBuffer::writeTo(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(this->data, this->size * sizeof(usint), 1, file);
}

void
ReadBuffer::moveBuffer(const usint* buffer)
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
  this->free_buffer = false;

  this->data = buffer;
  this->reset();
}

usint
ReadBuffer::reportSize() const
{
  usint bytes = sizeof(*this);
  if(this->free_buffer) { bytes += this->size * sizeof(usint); }
  return bytes;
}

void
ReadBuffer::exportTo(const std::string& filename) const
{
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "ReadBuffer::exportTo(): Cannot open output file " << filename << std::endl;
    return;
  }

  if(this->items == 0 || this->item_bits == 0) { this->writeTo(output); }
  else if(this->item_bits == 1)
  {
    for(usint i = 0; i < this->items; i++)
    {
      char temp = (this->isSet(i) ? '1' : '0');
      output.write(&temp, sizeof(temp));
    }
  }
  else
  {
    for(usint i = 0; i < this->items; i++)
    {
      usint temp = this->readItemConst(i);
      output.write((char*)(&temp), sizeof(temp));
    }
  }

  output.close();
}

//--------------------------------------------------------------------------

WriteBuffer::WriteBuffer(usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(true)
{
  this->data = new usint[words];
  memset(this->data, 0, this->size * sizeof(usint));
  this->reset();
}

WriteBuffer::WriteBuffer(usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(true)
{
  this->size = bitsToWords(this->items * this->item_bits);
  this->data = new usint[this->size];
  memset(this->data, 0, this->size * sizeof(usint));
  this->reset();
}

WriteBuffer::WriteBuffer(usint* buffer, usint words) :
  size(words),
  item_bits(1),
  items(0),
  free_buffer(false)
{
  this->data = buffer;
  this->reset();
}

WriteBuffer::WriteBuffer(usint* buffer, usint _items, usint item_size) :
  item_bits(item_size),
  items(_items),
  free_buffer(false)
{
  this->size = bitsToWords(this->items * this->item_bits);
  this->data = buffer;
  this->reset();
}

WriteBuffer::~WriteBuffer()
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
}

//--------------------------------------------------------------------------

usint*
WriteBuffer::rawData()
{
  return this->data;
}

usint*
WriteBuffer::stealData()
{
  this->free_buffer = false;
  return this->data;
}

ReadBuffer*
WriteBuffer::getReadBuffer()
{
  if(this->items > 0) { return this->getReadBuffer(this->items, this->item_bits); }
  else { return this->getReadBuffer(this->size); }
}

ReadBuffer*
WriteBuffer::getReadBuffer(usint words)
{
  words = std::min(words, this->size);

  ReadBuffer* buffer = new ReadBuffer(this->data, words);
  if(this->free_buffer)
  {
    buffer->claimData();
    this->free_buffer = false;
  }
  return buffer;
}

ReadBuffer*
WriteBuffer::getReadBuffer(usint _items, usint item_size)
{
  ReadBuffer* buffer = new ReadBuffer(this->data, _items, item_size);
  if(this->free_buffer)
  {
    buffer->claimData();
    this->free_buffer = false;
  }
  return buffer;
}

void
WriteBuffer::writeTo(std::ofstream& file) const
{
  file.write((char*)this->data, this->size * sizeof(usint));
}

void
WriteBuffer::writeTo(FILE* file) const
{
  if(file == 0) { return; }
  std::fwrite(this->data, this->size * sizeof(usint), 1, file);
}

void
WriteBuffer::moveBuffer(usint* buffer)
{
  if(this->free_buffer)
  {
    delete[] this->data;
  }
  this->free_buffer = false;

  this->data = buffer;
  this->reset();
}

usint
WriteBuffer::reportSize() const
{
  usint bytes = sizeof(*this);
  if(this->free_buffer) { bytes += this->size * sizeof(usint); }
  return bytes;
}

//--------------------------------------------------------------------------

} // namespace CSA
