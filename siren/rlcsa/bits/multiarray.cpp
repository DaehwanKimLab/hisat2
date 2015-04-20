#include <cstdlib>

#include "multiarray.h"


namespace CSA
{

//--------------------------------------------------------------------------

MultiArray*
MultiArray::createFixed(usint items, usint bits)
{
  return new FixedMultiArray(items, bits);
}

MultiArray*
MultiArray::createDelta(usint block_size)
{
  return new DeltaMultiArray(block_size);
}

MultiArray*
MultiArray::readFrom(std::ifstream& file)
{
  usint flags = 0;
  file.read((char*)&flags, sizeof(flags));
  if(flags & DELTA_FLAG) { return new DeltaMultiArray(file, flags); }
  else { return new FixedMultiArray(file, flags); }
}

MultiArray*
MultiArray::readFrom(FILE* file)
{
  if(file == 0) { return 0; }
  usint flags = 0;
  if(!std::fread(&flags, sizeof(flags), 1, file)) { return 0; }
  if(flags & DELTA_FLAG) { return new DeltaMultiArray(file, flags); }
  else { return new FixedMultiArray(file, flags); }
}

MultiArray::~MultiArray()
{
  delete this->array_encoder; this->array_encoder = 0;
  delete this->array_borders; this->array_borders = 0;
}

MultiArray::MultiArray(std::ifstream& file, usint flags) :
  status(error), encoding(fixed_width),
  items(0), arrays(0), current_array_begin(0),
  array_encoder(0), array_borders(0)
{
  if(flags & DELTA_FLAG) { this->encoding = delta_coded; }

  this->array_borders = new SuccinctVector(file);
  this->items = this->array_borders->getSize();
  this->arrays = this->array_borders->getNumberOfItems();
}

MultiArray::MultiArray(FILE* file, usint flags) :
  status(error), encoding(fixed_width),
  items(0), arrays(0), current_array_begin(0),
  array_encoder(0), array_borders(0)
{
  if(flags & DELTA_FLAG) { this->encoding = delta_coded; }

  this->array_borders = new SuccinctVector(file);
  this->items = this->array_borders->getSize();
  this->arrays = this->array_borders->getNumberOfItems();
}

MultiArray::MultiArray(usint _items, usint _item_bits) :
  status(error), encoding(fixed_width), items(0), arrays(1), current_array_begin(0),
  array_encoder(new SuccinctVector::Encoder(VECTOR_BLOCK_SIZE)), array_borders(0)
{
  this->array_encoder->setBit(0);
}

MultiArray::MultiArray(usint block_size) :
  status(error), encoding(delta_coded), items(0), arrays(1), current_array_begin(0),
  array_encoder(new SuccinctVector::Encoder(VECTOR_BLOCK_SIZE)), array_borders(0)
{
  this->array_encoder->setBit(0);
}

bool
MultiArray::writeTo(std::ofstream& file) const
{
  if(!(this->isOk()))
  {
    std::cerr << "MultiArray::writeTo(): Cannot write MultiArray with status ";
    if(this->status == error) { std::cerr << "error"; }
    else if(this->status == writable) { std::cerr << "writable"; }
    std::cerr << std::endl;
    return false;
  }

  usint flags = this->getFlags();
  file.write((char*)&flags, sizeof(flags));
  this->array_borders->writeTo(file);

  return true;
}

bool
MultiArray::writeTo(FILE* file) const
{
  if(!(this->isOk()))
  {
    std::cerr << "MultiArray::writeTo(): Cannot write MultiArray with status ";
    if(this->status == error) { std::cerr << "error"; }
    else if(this->status == writable) { std::cerr << "writable"; }
    std::cerr << std::endl;
    return false;
  }
  if(file == 0) { return false; }

  usint flags = this->getFlags();
  std::fwrite(&flags, sizeof(flags), 1, file);
  this->array_borders->writeTo(file);

  return true;
}

usint
MultiArray::reportSize() const
{
  usint bytes = 0;
  if(this->array_borders != 0) { bytes += this->array_borders->reportSize(); }
  return bytes;
}

//--------------------------------------------------------------------------

bool
MultiArray::nextArray()
{
  if(this->isFull() || this->items <= this->current_array_begin) { return false; }

  this->array_encoder->setBit(this->items);
  this->arrays++; this->current_array_begin = this->items;
  return true;
}

bool
MultiArray::finishWriting()
{
  if(!(this->isWritable())) { return false; }
  if(this->items <= this->current_array_begin)
  {
    std::cerr << "MultiArray::finishWriting(): Buffer ends with an empty array!" << std::endl;
    return false;
  }

  this->array_borders = new SuccinctVector(*(this->array_encoder), this->items);
  delete this->array_encoder; this->array_encoder = 0;

  this->encodeItems();
  this->status = ok;
  return true;
}

//--------------------------------------------------------------------------

FixedMultiArray::FixedMultiArray(std::ifstream& file, usint flags) :
  MultiArray(file, flags),
  fixed_buffer(0), fixed_items(0)
{
  usint item_bits = (flags & ITEM_BITS_MASK) >> ITEM_BITS_SHIFT;
  this->fixed_items = new ReadBuffer(file, this->items, item_bits);
  this->status = ok;
}

FixedMultiArray::FixedMultiArray(FILE* file, usint flags) :
  MultiArray(file, flags),
  fixed_buffer(0), fixed_items(0)
{
  usint item_bits = (flags & ITEM_BITS_MASK) >> ITEM_BITS_SHIFT;
  this->fixed_items = new ReadBuffer(file, this->items, item_bits);
  this->status = ok;
}

FixedMultiArray::FixedMultiArray(usint _items, usint bits) :
  MultiArray(items, bits),
  fixed_buffer(new WriteBuffer(_items, bits)), fixed_items(0)
{
  this->status = writable;
}

FixedMultiArray::~FixedMultiArray()
{
  delete this->fixed_buffer; this->fixed_buffer = 0;
  delete this->fixed_items; this->fixed_items = 0;
}

bool
FixedMultiArray::writeTo(std::ofstream& file) const
{
  if(!(MultiArray::writeTo(file))) { return false; }
  this->fixed_items->writeBuffer(file);
  return true;
}

bool
FixedMultiArray::writeTo(FILE* file) const
{
  if(!(MultiArray::writeTo(file))) { return false; }
  this->fixed_items->writeBuffer(file);
  return true;
}

usint
FixedMultiArray::reportSize() const
{
  usint bytes = sizeof(*this) + MultiArray::reportSize();
  if(this->fixed_items != 0) { bytes += this->fixed_items->reportSize(); }
  return bytes;
}

bool
FixedMultiArray::isFull() const
{
  if(!(this->isWritable())) { return true; }
  return !(this->fixed_buffer->hasNextItem());
}

void
FixedMultiArray::writeItem(usint item)
{
  this->fixed_buffer->writeItem(item);
  this->items++;
}

usint
FixedMultiArray::getFlags() const
{
  return this->fixed_items->getItemSize() << ITEM_BITS_SHIFT;
}

void
FixedMultiArray::encodeItems()
{
  this->fixed_items = this->fixed_buffer->getReadBuffer();
  delete this->fixed_buffer; this->fixed_buffer = 0;
}

//--------------------------------------------------------------------------

DeltaMultiArray::DeltaMultiArray(std::ifstream& file, usint flags) :
  MultiArray(file, flags),
  delta_buffer(0), delta_items(new Array(file))
{
  this->status = ok;
}

DeltaMultiArray::DeltaMultiArray(FILE* file, usint flags) :
  MultiArray(file, flags),
  delta_buffer(0), delta_items(new Array(file))
{
  this->status = ok;
}

DeltaMultiArray::DeltaMultiArray(usint block_size) :
  MultiArray(block_size),
  delta_buffer(new Array::Encoder(block_size)), delta_items(0)
{
  this->status = writable;
}

DeltaMultiArray::~DeltaMultiArray()
{
  delete this->delta_buffer; this->delta_buffer = 0;
  delete this->delta_items; this->delta_items = 0;
}

bool
DeltaMultiArray::writeTo(std::ofstream& file) const
{
  if(!(MultiArray::writeTo(file))) { return false; }
  this->delta_items->writeTo(file);
  return true;
}

bool
DeltaMultiArray::writeTo(FILE* file) const
{
  if(!(MultiArray::writeTo(file))) { return false; }
  this->delta_items->writeTo(file);
  return true;
}

usint
DeltaMultiArray::reportSize() const
{
  usint bytes = sizeof(*this) + MultiArray::reportSize();
  if(this->delta_items != 0) { bytes += this->delta_items->reportSize(); }
  return bytes;
}

bool
DeltaMultiArray::isFull() const
{
  return !(this->isWritable());
}

void
DeltaMultiArray::writeItem(usint item)
{
  this->delta_buffer->writeItem(item);
  this->items++;
}

usint
DeltaMultiArray::getFlags() const
{
  return DELTA_FLAG;
}

void
DeltaMultiArray::encodeItems()
{
  this->delta_items = new Array(*(this->delta_buffer));
  delete this->delta_buffer; this->delta_buffer = 0;
}

//--------------------------------------------------------------------------

MultiArray::Iterator::Iterator(const MultiArray& par) :
  parent(par),
  array_iter(*(par.array_borders)),
  pos(0), limit(0)
{
}

FixedMultiArray::Iterator::Iterator(const FixedMultiArray& par) :
  MultiArray::Iterator(par),
  fixed_iter(*(par.fixed_items))
{
}

DeltaMultiArray::Iterator::Iterator(const DeltaMultiArray& par) :
  MultiArray::Iterator(par),
  delta_iter(*(par.delta_items)),
  buffer(par.items, 0)
{
}

MultiArray::Iterator::~Iterator()
{
}

FixedMultiArray::Iterator::~Iterator()
{
}

DeltaMultiArray::Iterator::~Iterator()
{
}

void
FixedMultiArray::Iterator::goToItem(usint array, usint index)
{
  this->pos = std::min(this->array_iter.select(array) + index, this->parent.getSize());
  this->fixed_iter.goToItem(this->pos);
}

void
DeltaMultiArray::Iterator::goToItem(usint array, usint index)
{
  this->pos = std::min(this->array_iter.select(array) + index, this->parent.getSize());
  this->buffer = pair_type(this->pos, this->delta_iter.readItem(this->pos));
}

usint
FixedMultiArray::Iterator::readItem(usint array, usint index)
{
  this->goToItem(array, index);
  usint value = this->fixed_iter.readItem(this->pos); this->pos++;
  return value;
}

usint
DeltaMultiArray::Iterator::readItem(usint array, usint index)
{
  this->goToItem(array, index); this->pos++;
  return this->buffer.second;
}

usint
FixedMultiArray::Iterator::nextItem()
{
  usint value = this->fixed_iter.readItem(); this->pos++;
  return value;
}

usint
DeltaMultiArray::Iterator::nextItem()
{
  usint value = (this->buffer.first == this->pos ? this->buffer.second : this->delta_iter.nextItem());
  this->pos++;
  return value;
}

void
MultiArray::Iterator::setEnd(usint array, usint index)
{
  this->limit = std::min(this->array_iter.select(array) + index, this->parent.items);
}

//--------------------------------------------------------------------------

} // namespace CSA
