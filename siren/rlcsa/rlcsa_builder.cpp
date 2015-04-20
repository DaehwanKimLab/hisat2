#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "rlcsa_builder.h"
#include "misc/utils.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif


namespace CSA
{


RLCSABuilder::RLCSABuilder(usint _block_size, usint _sample_rate, usint _buffer_size, usint _threads) :
  block_size(_block_size), sample_rate(_sample_rate), buffer_size(_buffer_size),
  threads(_threads),
  buffer(0)
{
  this->reset();
  this->build_time = this->search_time = this->sort_time = this->merge_time = 0.0;
}

RLCSABuilder::~RLCSABuilder()
{
  delete this->index;
  delete[] this->buffer;
}

//--------------------------------------------------------------------------

void
RLCSABuilder::insertSequence(char* sequence, usint length, bool delete_sequence)
{
  if(sequence == 0 || length == 0 || !this->ok)
  {
    if(delete_sequence) { delete[] sequence; }
    return;
  }

  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(this->threads);
  #endif

  if(this->buffer == 0)
  {
    this->addLongSequence((uchar*)sequence, length, delete_sequence);
    return;
  }

  if(this->buffer_size - this->chars > length)
  {
    memcpy(this->buffer + this->chars, sequence, length);
    if(delete_sequence) { delete[] sequence; }
    this->chars += length;
    this->buffer[this->chars] = 0;
    this->chars++;
  }
  else
  {
    if(this->chars > 0) { this->flush(); }
    if(length >= this->buffer_size - 1)
    {
      this->addLongSequence((uchar*)sequence, length, delete_sequence);
    }
    else
    {
      memcpy(this->buffer + this->chars, sequence, length);
      if(delete_sequence) { delete[] sequence; }
      this->chars += length;
      this->buffer[this->chars] = 0;
      this->chars++;
    }
  }
}

void
RLCSABuilder::insertFromFile(const std::string& base_name)
{
  if(!this->ok) { return; }

  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(this->threads);
  #endif

  this->flush();

  std::ifstream input(base_name.c_str(), std::ios_base::binary);
  if(!input) { return; }
  RLCSA* increment = new RLCSA(base_name);
  usint data_size = increment->getSize() + increment->getNumberOfSequences();
  uchar* data = new uchar[data_size];
  input.read((char*)data, data_size);
  input.close();

  this->addRLCSA(increment, data, data_size, true);
}

void
RLCSABuilder::insertCollection(const std::string& base_name)
{
  if(!this->ok) { return; }

  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(this->threads);
  #endif

  this->flush();

  std::ifstream input(base_name.c_str(), std::ios_base::binary);
  if(!input) { return; }

  usint data_size = fileSize(input);
  uchar* data = new uchar[data_size];
  input.read((char*)data, data_size);
  input.close();

  if(this->index == 0)
  {
    this->addCollection(data, data_size, true);
    return;
  }

  std::vector<usint> end_markers;
  usint* ranks = this->getRanks(data, data_size, end_markers);

  double start = readTimer();
  #pragma omp parallel for schedule(static)
  for(usint i = 0; i < data_size; i++) { ranks[i] += end_markers.size() + 1 + data[i]; }
  for(usint i = 0; i < end_markers.size(); i++)
  {
    ranks[end_markers[i]] = this->index->getNumberOfSequences() + i;
  }
  RLCSA* increment = new RLCSA(data, ranks, data_size, this->block_size, this->sample_rate, this->threads, false);
  #pragma omp parallel for schedule(static)
  for(usint i = 0; i < data_size; i++) { ranks[i] -= data[i]; }
  delete[] data;
  double mark = readTimer();
  this->build_time += mark - start;

  parallelSort(ranks, ranks + data_size);
  #pragma omp parallel for schedule(static)
  for(usint i = end_markers.size(); i < data_size; i++) { ranks[i] += i - end_markers.size(); }
  this->sort_time += readTimer() - mark;

  this->mergeRLCSA(increment, ranks, data_size);
}

//--------------------------------------------------------------------------

RLCSA*
RLCSABuilder::getRLCSA()
{
  if(this->chars > 0) { this->flush(); }

  RLCSA* temp = this->index;
  this->reset();

  return temp;
}

char*
RLCSABuilder::getBWT(usint& length)
{
  this->flush();

  if(this->index == 0 || !(this->ok))
  {
    length = 0;
    return 0;
  }

  length = this->index->getSize() + this->index->getNumberOfSequences();
  return (char*)(this->index->readBWT());
}

bool
RLCSABuilder::isOk()
{
  return this->ok;
}

double
RLCSABuilder::getBuildTime()
{
  return this->build_time;
}

double
RLCSABuilder::getSearchTime()
{
  return this->search_time;
}

double
RLCSABuilder::getSortTime()
{
  return this->sort_time;
}

double
RLCSABuilder::getMergeTime()
{
  return this->merge_time;
}

//--------------------------------------------------------------------------

void
RLCSABuilder::flush()
{
  if(this->buffer == 0 || this->chars == 0) { return; }

  double start = readTimer();
  RLCSA* temp = new RLCSA(this->buffer, this->chars, this->block_size, this->sample_rate, this->threads, (this->index == 0));
  this->build_time += readTimer() - start;
  this->addRLCSA(temp, this->buffer, this->chars, (this->index != 0));

  this->chars = 0;
  this->buffer = new uchar[this->buffer_size];
}

void
RLCSABuilder::reset()
{
  this->index = 0;

  if(this->buffer_size != 0)
  {
    delete[] this->buffer;
    this->buffer = new uchar[this->buffer_size];
  }
  this->chars = 0;

  this->ok = true;
}

//--------------------------------------------------------------------------

void
RLCSABuilder::addRLCSA(RLCSA* increment, uchar* sequence, usint length, bool delete_sequence)
{
  if(this->index == 0)
  {
    if(delete_sequence) { delete[] sequence; }
    this->setRLCSA(increment);
    return;
  }

  std::vector<usint> end_markers;
  usint* ranks = this->getRanks(sequence, length, end_markers);
  if(delete_sequence) { delete[] sequence; }

  double mark = readTimer();
  parallelSort(ranks, ranks + length);
  #pragma omp parallel for schedule(static)
  for(usint i = 0; i < length; i++) { ranks[i] += i + 1; }
  this->sort_time += readTimer() - mark;

  this->mergeRLCSA(increment, ranks, length);
}

void
RLCSABuilder::setRLCSA(RLCSA* new_index)
{
  this->index = new_index;
  this->ok &= this->index->isOk();
}

void
RLCSABuilder::mergeRLCSA(RLCSA* increment, usint* ranks, usint length)
{
  double mark = readTimer();

  RLCSA* merged = new RLCSA(*(this->index), *increment, ranks, this->block_size, this->threads);
  delete[] ranks;
  delete this->index;
  delete increment;
  this->index = merged;

  this->merge_time += readTimer() - mark;

  this->ok &= this->index->isOk();
}

//--------------------------------------------------------------------------

usint*
RLCSABuilder::getRanks(uchar* sequence, usint length, std::vector<usint>& end_markers)
{
  double start = readTimer();

  usint sequences = 0;
  for(usint i = 0; i < length; i++)
  {
    if(sequence[i] == 0) { end_markers.push_back(i); sequences++; }
  }

  usint* ranks = new usint[length];
  #ifdef MULTITHREAD_SUPPORT
  usint chunk = std::max((usint)1, sequences / (8 * this->threads));
  #endif
  #pragma omp parallel for schedule(dynamic, chunk)
  for(usint i = 0; i < sequences; i++)
  {
    usint begin = (i > 0 ? end_markers[i - 1] + 1 : 0);
    this->index->reportPositions(sequence + begin, end_markers[i] - begin, ranks + begin);
  }

  this->index->strip();

  this->search_time += readTimer() - start;
  return ranks;
}

//--------------------------------------------------------------------------

void
RLCSABuilder::addLongSequence(uchar* sequence, usint length, bool delete_sequence)
{
  double start = readTimer();
  RLCSA* temp = new RLCSA(sequence, length, this->block_size, this->sample_rate, this->threads, 0, false);
  this->build_time += readTimer() - start;
  this->addRLCSA(temp, sequence, length, delete_sequence); 
}

void
RLCSABuilder::addCollection(uchar* sequence, usint length, bool delete_sequence)
{
  double start = readTimer();
  RLCSA* temp = new RLCSA(sequence, length, this->block_size, this->sample_rate, this->threads, delete_sequence);
  this->build_time += readTimer() - start;
  this->setRLCSA(temp);
}

//--------------------------------------------------------------------------

} // namespace CSA
