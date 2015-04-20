#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>

#include "rlcsa.h"
#include "misc/utils.h"
#include "bits/vectors.h"


#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif



namespace CSA
{


RLCSA::RLCSA(const std::string& base_name, bool print) :
  ok(false),
  alphabet(0),
  sa_samples(0), support_locate(false), support_display(false),
  end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }

  std::string array_name = base_name + ARRAY_EXTENSION;
  std::ifstream array_file(array_name.c_str(), std::ios_base::binary);
  if(!array_file)
  {
    std::cerr << "RLCSA: Error opening Psi array file!" << std::endl;
    return;
  }

  usint distribution[CHARS];
  array_file.read((char*)distribution, CHARS * sizeof(usint));
  this->alphabet = new Alphabet(distribution); this->data_size = this->alphabet->getDataSize();

  Parameters parameters;
  parameters.read(base_name + PARAMETERS_EXTENSION);
  for(usint c = 0; c < CHARS; c++)
  {
    if(this->alphabet->hasChar(c)) { this->array[c] = new PsiVector(array_file); }
  }

  this->end_points = new DeltaVector(array_file);
  this->number_of_sequences = this->end_points->getNumberOfItems();

  array_file.read((char*)&(this->sample_rate), sizeof(this->sample_rate));
  array_file.close();

  if(parameters.get(SUPPORT_LOCATE) || parameters.get(SUPPORT_DISPLAY))
  {
    std::string sa_sample_name = base_name + SA_SAMPLES_EXTENSION;
    std::ifstream sa_sample_file(sa_sample_name.c_str(), std::ios_base::binary);
    if(!sa_sample_file)
    {
      std::cerr << "RLCSA: Error opening suffix array sample file!" << std::endl;
      return;
    }
    bool weighted = parameters.get(WEIGHTED_SAMPLES);
    this->sa_samples = new SASamples(sa_sample_file, this->sample_rate, weighted);
    sa_sample_file.close();

    this->support_locate = this->sa_samples->supportsLocate();
    this->support_display = this->sa_samples->supportsDisplay();
  }

  if(print) { parameters.print(); }

  this->ok = true;
}

RLCSA::RLCSA(uchar* data, usint bytes, usint block_size, usint sa_sample_rate, usint threads, bool delete_data) :
  ok(false),
  alphabet(0),
  sa_samples(0), support_locate(false), support_display(false),
  sample_rate(sa_sample_rate), end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }
 
  if(!data || bytes == 0)
  {
    std::cerr << "RLCSA: No input data given!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }
  if(block_size < 2 * sizeof(usint) || block_size % sizeof(usint) != 0)
  {
    std::cerr << "RLCSA: Block size must be a multiple of " << sizeof(usint) << " bytes!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }

  this->buildRLCSA(data, 0, bytes, block_size, threads, 0, true, delete_data);
}

RLCSA::RLCSA(uchar* data, usint* ranks, usint bytes, usint block_size, usint sa_sample_rate, usint threads, bool delete_data) :
  ok(false),
  sa_samples(0), support_locate(false), support_display(false),
  sample_rate(sa_sample_rate), end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }
 
  if(!data || !ranks || bytes == 0)
  {
    std::cerr << "RLCSA: No input data given!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }
  if(block_size < 2 * sizeof(usint) || block_size % sizeof(usint) != 0)
  {
    std::cerr << "RLCSA: Block size must be a multiple of " << sizeof(usint) << " bytes!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }

  this->buildRLCSA(data, ranks, bytes, block_size, threads, 0, true, delete_data);
}

RLCSA::RLCSA(uchar* data, usint bytes, usint block_size, usint sa_sample_rate, usint threads, Sampler* sampler, bool delete_data) :
  ok(false),
  alphabet(0),
  sa_samples(0), support_locate(false), support_display(false),
  sample_rate(sa_sample_rate),
  end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }

  if(!data || bytes == 0)
  {
    std::cerr << "RLCSA: No input data given!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }
  if(block_size < 2 * sizeof(usint) || block_size % sizeof(usint) != 0)
  {
    std::cerr << "RLCSA: Block size must be a multiple of " << sizeof(usint) << " bytes!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }
  if(sampler != 0 && sampler->getStatus() != Sampler::SAMPLED)
  {
    std::cerr << "RLCSA: No samples given!" << std::endl;
    if(delete_data) { delete[] data; }
    return;
  }

  this->buildRLCSA(data, 0, bytes, block_size, threads, sampler, false, delete_data);
}

RLCSA::RLCSA(RLCSA& index, RLCSA& increment, usint* positions, usint block_size, usint threads) :
  ok(false),
  alphabet(0),
  sa_samples(0), support_locate(false), support_display(false),
  end_points(0)
{
  for(usint c = 0; c < CHARS; c++) { this->array[c] = 0; }

  if(!index.isOk() || !increment.isOk())
  {
    return; // Fail silently. Actual error has already been reported.
  }
  if(positions == 0)
  {
    std::cerr << "RLCSA: Positions for insertions not available!" << std::endl;
    return;
  }
  if(index.sample_rate != increment.sample_rate)
  {
    std::cerr << "RLCSA: Cannot combine indexes with different sample rates!" << std::endl;
    return;
  }

  index.strip();
  increment.strip();

  // Build character tables etc.
  usint distribution[CHARS];
  for(usint c = 0; c < CHARS; c++)
  {
    distribution[c] = index.alphabet->countOf(c) + increment.alphabet->countOf(c);
  }
  this->alphabet = new Alphabet(distribution); this->data_size = this->alphabet->getDataSize();
  this->sample_rate = index.sample_rate;
  this->number_of_sequences = index.number_of_sequences + increment.number_of_sequences;


  // Merge end points, SA samples, and Psi.
  usint psi_size = this->data_size + this->number_of_sequences;
  bool should_be_ok = true;

  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif
  #pragma omp parallel for schedule(dynamic, 1)
  for(int c = -2; c < (int)CHARS; c++)
  {
    if(c == -2)      { this->mergeEndPoints(index, increment); }
    else if(c == -1) { this->mergeSamples(index, increment, positions);  }
    else if(this->alphabet->hasChar(c) != 0)
    {
      this->array[c] = mergeVectors<PsiVector>(index.array[c], increment.array[c], positions, increment.data_size + increment.number_of_sequences, psi_size, block_size);
      index.array[c] = 0;
      increment.array[c] = 0;

      if(this->array[c] == 0)
      {
        std::cerr << "RLCSA: Merge failed for vectors " << c << std::endl;
        should_be_ok = false;
      }
    }
  }

  this->ok = should_be_ok;
}

RLCSA::~RLCSA()
{
  for(usint c = 0; c < CHARS; c++) { delete this->array[c]; this->array[c] = 0; }
  delete this->alphabet; this->alphabet = 0;
  delete this->sa_samples; this->sa_samples = 0;
  delete this->end_points; this->end_points = 0;
}

//--------------------------------------------------------------------------

void
RLCSA::writeTo(const std::string& base_name) const
{
  std::string array_name = base_name + ARRAY_EXTENSION;
  std::ofstream array_file(array_name.c_str(), std::ios_base::binary);
  if(!array_file)
  {
    std::cerr << "RLCSA: Error creating Psi array file!" << std::endl;
    return;
  }

  this->alphabet->writeTo(array_file);
  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c] != 0)
    {
      this->array[c]->writeTo(array_file);
    }
  }

  this->end_points->writeTo(array_file);
  array_file.write((char*)&(this->sample_rate), sizeof(this->sample_rate));
  array_file.close();

  if(this->sa_samples != 0)
  {
    std::string sa_sample_name = base_name + SA_SAMPLES_EXTENSION;
    std::ofstream sa_sample_file(sa_sample_name.c_str(), std::ios_base::binary);
    if(!sa_sample_file)
    {
      std::cerr << "RLCSA: Error creating suffix array sample file!" << std::endl;
      return;
    }

    this->sa_samples->writeTo(sa_sample_file);
    sa_sample_file.close();
  }

  Parameters parameters;
  parameters.set(RLCSA_BLOCK_SIZE.first, this->getBlockSize() * sizeof(usint));
  parameters.set(SAMPLE_RATE.first, this->sample_rate);
  parameters.set(SUPPORT_LOCATE.first, this->support_locate);
  parameters.set(SUPPORT_DISPLAY.first, this->support_display);
  if(this->sa_samples != 0 && this->sa_samples->isWeighted())
  {
    parameters.set(WEIGHTED_SAMPLES.first, 1);
  }
  else { parameters.set(WEIGHTED_SAMPLES); }
  parameters.write(base_name + PARAMETERS_EXTENSION);
}

//--------------------------------------------------------------------------

pair_type
RLCSA::count(const std::string& pattern) const
{
  if(pattern.length() == 0) { return this->getSARange(); }

  std::string::const_reverse_iterator iter = pattern.rbegin();
  pair_type index_range = this->getCharRange((uchar)*iter);
  if(isEmpty(index_range)) { return index_range; }

  for(++iter; iter != pattern.rend(); ++iter)
  {
    index_range = this->LF(index_range, (uchar)*iter);
    if(isEmpty(index_range)) { return EMPTY_PAIR; }
  }
  this->convertToSARange(index_range);

  return index_range;
}

//--------------------------------------------------------------------------

void
RLCSA::reportPositions(uchar* data, usint length, usint* positions) const
{
  if(data == 0 || length == 0 || positions == 0) { return; }

  PsiVector::Iterator** iters = this->getIterators();

  usint current = this->number_of_sequences - 1;
  positions[length] = current; // "immediately after current"
  for(sint i = (sint)(length - 1); i >= 0; i--)
  {
    usint c = (usint)data[i];
    if(this->array[c] != 0)
    {
      current = this->LF(current, c, *(iters[c]));
    }
    else
    {
      if(c < this->alphabet->getFirstChar()) // No previous characters either.
      {
        current = this->number_of_sequences - 1;
      }
      else
      {
        current = this->alphabet->cumulative(c) - 1 + this->number_of_sequences;
      }
    }
    positions[i] = current; // "immediately after current"
  }

  this->deleteIterators(iters);
}

//--------------------------------------------------------------------------

usint*
RLCSA::locate(pair_type range, bool direct, bool steps) const
{
  if(!(this->support_locate) || isEmpty(range) || range.second >= this->data_size) { return 0; }

  usint* data = new usint[length(range)];
  if(direct) { this->directLocate(range, data, steps); }
  else       { this->locateUnsafe(range, data, steps); }

  return data;
}

usint*
RLCSA::locate(pair_type range, usint* data, bool direct, bool steps) const
{
  if(!(this->support_locate) || isEmpty(range) || range.second >= this->data_size || data == 0) { return 0; }

  if(direct) { this->directLocate(range, data, steps); }
  else       { this->locateUnsafe(range, data, steps); }

  return data;
}

usint
RLCSA::locate(usint index, bool steps) const
{
  if(!(this->support_locate) || index >= this->data_size) { return (steps ? 0 : this->data_size); }

  return this->directLocate(index + this->number_of_sequences, steps);
}

void
RLCSA::directLocate(pair_type range, usint* data, bool steps) const
{
  this->convertToBWTRange(range);
  for(usint i = 0, j = range.first; j <= range.second; i++, j++)
  {
    data[i] = this->directLocate(j, steps);
  }
}

usint
RLCSA::directLocate(usint index, bool steps) const
{
  usint offset = 0;
  while(true)
  {
    if(this->hasImplicitSample(index))
    {
      if(steps) { return offset; }
      else      { return this->getImplicitSample(index) - offset; }
    }
    index -= this->number_of_sequences;
    if(this->sa_samples->isSampled(index))
    {
      return (steps ? offset : this->sa_samples->getSampleAt(index) - offset);
    }
    index = this->psi(index);
    offset++;
  }
}

void
RLCSA::locateUnsafe(pair_type range, usint* data, bool steps) const
{
  this->convertToBWTRange(range);
  usint items = length(range);
  usint* offsets = new usint[items];
  bool* finished = new bool[items];  // FIXME This could be more space efficient...

  PsiVector::Iterator** iters = this->getIterators();

  for(usint i = 0, j = range.first; i < items; i++, j++)
  {
    data[i] = j;
    offsets[i] = 0;
    finished[i] = false;
  }

  bool found = false;
  while(!found)
  {
    found = true;
    pair_type run = EMPTY_PAIR;
    for(usint i = 0; i < items; i++)
    {
      if(finished[i])
      {
        continue; // The run might continue after this.
      }
      else if(isEmpty(run))
      {
        run = pair_type(i, i);
      }
      else if(data[i] - data[run.first] == i - run.first)
      {
        run.second = i;
      }
      else
      {
        found &= this->processRun(run, data, offsets, finished, iters, steps);
        run = pair_type(i, i);
      }
    }
    if(!isEmpty(run)) { found &= this->processRun(run, data, offsets, finished, iters, steps); }
  }

  this->deleteIterators(iters);
  delete[] offsets;
  delete[] finished;
}

bool
RLCSA::processRun(pair_type run, usint* data, usint* offsets, bool* finished, PsiVector::Iterator** iters, bool steps) const
{
  bool found = true;
  usint run_start = 0, run_left = 0;
  pair_type next_sample = pair_type(0, 0);

  for(usint i = run.first; i <= run.second; i++)
  {
    if(finished[i])
    {
      if(run_left > 0) { run_left--; }
      continue;
    }
    if(data[i] < this->number_of_sequences) // Implicit sample here.
    {
      DeltaVector::Iterator iter(*(this->end_points));
      data[i] = (steps ? offsets[i] : iter.select(data[i]) + 1 - offsets[i]);
      finished[i] = true;
      if(run_left > 0) { run_left--; }
      continue;
    }
    if(next_sample.first < data[i]) // Need another sample.
    {
      next_sample = this->sa_samples->getFirstSampleAfter(data[i] - this->number_of_sequences);
      next_sample.first += this->number_of_sequences;
    }
    if(data[i] < next_sample.first) // No sample found for current position.
    {
      if(run_left > 0)
      {
        data[i] = data[run_start] + i - run_start;
        run_left--;
      }
      else
      {
        pair_type value = this->psi(data[i] - this->number_of_sequences, run.second - i, iters);
        data[i] = value.first;
        run_left = value.second;
        run_start = i;
      }
      offsets[i]++;
      found = false;
    }
    else  // Sampled position found.
    {
      data[i] = (steps ? offsets[i] : this->sa_samples->getSample(next_sample.second) - offsets[i]);
      finished[i] = true;
      if(run_left > 0) { run_left--; }
    }
  }
  return found;
}

//--------------------------------------------------------------------------

uchar*
RLCSA::display(usint sequence, bool include_end_marker) const
{
  if(!(this->support_display)) { return 0; }

  pair_type seq_range = this->getSequenceRange(sequence);
  if(isEmpty(seq_range)) { return 0; }

  uchar* data = new uchar[length(seq_range) + include_end_marker];
  this->displayUnsafe(seq_range, data);
  if(include_end_marker) { data[length(seq_range)] = 0; }

  return data;
}

uchar*
RLCSA::display(usint sequence, pair_type range) const
{
  if(!(this->support_display) || isEmpty(range)) { return 0; }

  pair_type seq_range = this->getSequenceRange(sequence);
  if(isEmpty(seq_range)) { return 0; }

  range.first += seq_range.first; range.second += seq_range.first;
  if(range.second > seq_range.second) { return 0; }

  uchar* data = new uchar[length(range)];
  this->displayUnsafe(range, data);

  return data;
}

uchar*
RLCSA::display(usint sequence, pair_type range, uchar* data) const
{
  if(!(this->support_display) || isEmpty(range) || data == 0) { return 0; }

  pair_type seq_range = this->getSequenceRange(sequence);
  if(isEmpty(seq_range)) { return 0; }

  range.first += seq_range.first; range.second += seq_range.first;
  if(range.second > seq_range.second) { return 0; }

  this->displayUnsafe(range, data);
  return data;
}

uchar*
RLCSA::display(usint position, usint len, usint context, usint& result_length) const
{
  if(!(this->support_display)) { return 0; }

  pair_type range = this->getSequenceRangeForPosition(position);
  if(isEmpty(range)) { return 0; }

  range.first = position - std::min(context, position - range.first);
  range.second = std::min(range.second, position + len + context - 1);
  result_length = length(range);
  if(isEmpty(range)) { return 0; }

  uchar* data = new uchar[length(range)];
  this->displayUnsafe(range, data);

  return data;
}

usint
RLCSA::displayPrefix(usint sequence, usint len, uchar* data) const
{
  if(!(this->support_display) || len == 0 || data == 0) { return 0; }

  pair_type seq_range = this->getSequenceRange(sequence);
  if(isEmpty(seq_range)) { return 0; }

  pair_type range(seq_range.first, std::min(seq_range.second, seq_range.first + len - 1));

  this->displayUnsafe(range, data);
  return length(range);
}

usint
RLCSA::displayFromPosition(usint index, usint max_len, uchar* data) const
{
  if(max_len == 0 || data == 0 || index >= this->data_size) { return 0; }

  for(usint i = 0; i < max_len; i++)
  {
    data[i] = this->getCharacter(index);
    index = this->psiUnsafe(index, data[i]);
    if(index < this->number_of_sequences) { return i + 1; }
    index -= this->number_of_sequences;
  }

  return max_len;
}

void
RLCSA::displayUnsafe(pair_type range, uchar* data, bool get_ranks, usint* ranks) const
{
  pair_type res = this->sa_samples->inverseSA(range.first);
  usint i = res.first, pos = res.second;

  if(length(range) >= 1024)
  {
    PsiVector::Iterator** iters = this->getIterators();
    for(; i < range.first; i++)
    {
      pos = this->psi(pos, iters) - this->number_of_sequences;
    }
    for(; i <= range.second; i++)
    {
      usint c = this->getCharacter(pos);
      data[i - range.first] = c;
      if(get_ranks) { ranks[i - range.first] = pos + this->number_of_sequences; }
      pos = this->psiUnsafe(pos, c, *(iters[c])) - this->number_of_sequences;
    }
    this->deleteIterators(iters);
  }
  else
  {
    for(; i < range.first; i++)
    {
      pos = this->psi(pos) - this->number_of_sequences;
    }
    for(; i <= range.second; i++)
    {
      usint c = this->getCharacter(pos);
      data[i - range.first] = c;
      if(get_ranks) { ranks[i - range.first] = pos + this->number_of_sequences; }
      pos = this->psiUnsafe(pos, c) - this->number_of_sequences;
    }
  }
}

//--------------------------------------------------------------------------

pair_type
RLCSA::getSARange() const
{
  return pair_type(0, this->data_size - 1);
}

pair_type
RLCSA::getBWTRange() const
{
  return pair_type(this->number_of_sequences, this->number_of_sequences + this->data_size - 1);
}

pair_type
RLCSA::getCharRange(usint c) const
{
  if(c >= CHARS) { return EMPTY_PAIR; }
  pair_type index_range = this->alphabet->getRange(c);
  this->convertToBWTRange(index_range);
  return index_range;
}

void
RLCSA::convertToSARange(pair_type& bwt_range) const
{
  bwt_range.first -= this->number_of_sequences;
  bwt_range.second -= this->number_of_sequences;
}

void
RLCSA::convertToBWTRange(pair_type& sa_range) const
{
  sa_range.first += this->number_of_sequences;
  sa_range.second += this->number_of_sequences;
}

void
RLCSA::convertToSARange(std::vector<pair_type>& bwt_ranges) const
{
  for(std::vector<pair_type>::iterator iter = bwt_ranges.begin(); iter != bwt_ranges.end(); ++iter)
  {
    this->convertToSARange(*iter);
  }
}

pair_type
RLCSA::LF(pair_type range, usint c) const
{
  if(c >= CHARS || this->array[c] == 0) { return EMPTY_PAIR; }
  PsiVector::Iterator iter(*(this->array[c]));

  usint start = this->alphabet->cumulative(c) + this->number_of_sequences - 1;
  range.first = start + iter.rank(range.first, true);
  range.second = start + iter.rank(range.second);

  return range;
}

std::vector<usint>*
RLCSA::locateRange(pair_type range) const
{
  std::vector<usint>* results = new std::vector<usint>;
  if(!(this->support_locate)) { return results; }

  this->locateRange(range, *results);
  removeDuplicates(results, false);

  return results;
}

void
RLCSA::locateRange(pair_type range, std::vector<usint>& vec) const
{
  if(isEmpty(range) || range.second >= this->data_size) { return; }

  usint* data = new usint[length(range)];
  this->locateUnsafe(range, data, false);
  for(usint i = 0; i < length(range); i++) { vec.push_back(data[i]); }
  delete[] data;
}

std::vector<usint>*
RLCSA::locateRanges(std::vector<pair_type>& ranges) const
{
  std::vector<usint>* results = new std::vector<usint>;
  if(!(this->support_locate)) { return results; }

  for(std::vector<pair_type>::iterator iter = ranges.begin(); iter != ranges.end(); ++iter)
  {
    this->locateRange(*iter, *results);
  }
  removeDuplicates(results, false);

  return results;
}

//--------------------------------------------------------------------------

pair_type
RLCSA::getSequenceRange(usint number) const
{
  if(number >= this->number_of_sequences) { return EMPTY_PAIR; }

  pair_type result;
  DeltaVector::Iterator iter(*(this->end_points));
  if(number == 0)
  {
    result.first = 0;
    result.second = iter.select(number);
  }
  else
  {
    result.first = nextMultipleOf(this->sample_rate, iter.select(number - 1));
    result.second = iter.selectNext();
  }

  return result;
}

usint
RLCSA::getSequenceForPosition(usint value) const
{
  if(value == 0) { return 0; }
  DeltaVector::Iterator iter(*(this->end_points));
  return iter.rank(value - 1);
}

pair_type
RLCSA::getSequenceRangeForPosition(usint value) const
{
  return this->getSequenceRange(this->getSequenceForPosition(value));
}

usint*
RLCSA::getSequenceForPosition(usint* values, usint len) const
{
  if(values == 0) { return 0; }

  DeltaVector::Iterator iter(*(this->end_points));
  for(usint i = 0; i < len; i++)
  {
    if(values[i] > 0) { values[i] = iter.rank(values[i] - 1); }
  }

  return values;
}

pair_type
RLCSA::getRelativePosition(usint value) const
{
  DeltaVector::Iterator iter(*(this->end_points));
  pair_type result(0, value);

  if(value > 0) { result.first = iter.rank(value - 1); }
  if(result.first > 0)
  {
    result.second -= nextMultipleOf(this->sample_rate, iter.select(result.first - 1));
  }

  return result;
}

//--------------------------------------------------------------------------

uchar*
RLCSA::readBWT() const
{
  return this->readBWT(pair_type(0, this->data_size + this->number_of_sequences - 1));
}

uchar*
RLCSA::readBWT(pair_type range) const
{
  if(isEmpty(range) || range.second >= this->data_size + this->number_of_sequences) { return 0; }

  usint n = length(range);

  uchar* bwt = new uchar[n];
  memset(bwt, 0, n);

  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c] != 0)
    {
      PsiVector::Iterator iter(*(this->array[c]));
      usint pos = iter.valueAfter(range.first).first;
      while(pos <= range.second)
      {
        bwt[pos - range.first] = c;
        pos = iter.selectNext();
      }
    }
  }

  return bwt;
}

usint
RLCSA::countRuns() const
{
  usint runs = 0;
  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c] != 0)
    {
      PsiVector::Iterator iter(*(this->array[c]));
      runs += iter.countRuns();
    }
  }

  return runs;
}

//--------------------------------------------------------------------------

SuffixArray*
RLCSA::getSuffixArrayForSequence(usint number) const
{
  if(!this->supportsDisplay()) { return 0; }

  pair_type seq_range = this->getSequenceRange(number);
  if(isEmpty(seq_range)) { return 0; }

  uchar* data = new uchar[length(seq_range) + 1];
  usint* ranks = new usint[length(seq_range) + 1];
  this->displayUnsafe(seq_range, data, true, ranks);
  data[length(seq_range)] = 0; ranks[length(seq_range)] = number;

  return new SuffixArray(data, ranks, length(seq_range) + 1, 1);
}

//--------------------------------------------------------------------------

PsiVector::Iterator**
RLCSA::getIterators() const
{
  PsiVector::Iterator** iters = new PsiVector::Iterator*[CHARS];
  for(usint i = 0; i < CHARS; i++)
  {
    if(this->array[i] == 0) { iters[i] = 0; }
    else                    { iters[i] = new PsiVector::Iterator(*(this->array[i])); }
  }
  return iters;
}

void
RLCSA::deleteIterators(PsiVector::Iterator** iters) const
{
  if(iters == 0) { return; }
  for(usint i = 0; i < CHARS; i++) { delete iters[i]; }
  delete[] iters;
}

//--------------------------------------------------------------------------

usint
RLCSA::reportSize(bool print) const
{
  usint bytes = 0, temp = 0, bwt = 0;

  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c])
    {
      bytes += this->array[c]->reportSize();
      bwt += this->array[c]->getCompressedSize();
    }
  }
  bytes += sizeof(*this) + this->alphabet->reportSize() + this->end_points->reportSize();
  if(print)
  {
    std::cout << "RLCSA:           " << (bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "  BWT only:      " << (bwt / (double)MEGABYTE) << " MB" << std::endl;
  }

  if(this->support_locate || this->support_display)
  {
    temp = this->sa_samples->reportSize();
    if(print) { std::cout << "SA samples:      " << (temp / (double)MEGABYTE) << " MB" << std::endl; }
    bytes += temp;
  }

  if(print)
  {
    std::cout << "Total size:      " << (bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << std::endl;
  }

  return bytes;
}

void
RLCSA::printInfo() const
{
  double megabytes = this->data_size / (double)MEGABYTE;

  std::cout << "Sequences:       " << this->number_of_sequences << std::endl;
  std::cout << "Original size:   " << megabytes << " MB" << std::endl;
  std::cout << "Block size:      " << (this->getBlockSize() * sizeof(usint)) << " bytes" << std::endl;
  if(this->support_locate || this->support_display)
  {
    std::cout << "Sample rate:     " << this->sample_rate;
    if(this->sa_samples->isWeighted()) { std::cout << " (weighted)"; }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

//--------------------------------------------------------------------------

PLCPVector*
RLCSA::buildPLCP(usint block_size) const
{
  if(block_size < 2 * sizeof(usint) || block_size % sizeof(usint) != 0)
  {
    std::cerr << "PLCP: Block size must be a multiple of " << sizeof(usint) << " bytes!" << std::endl;
    return 0;
  }

  PsiVector::Iterator** iters = this->getIterators();

  PLCPVector::Encoder plcp(block_size);
  std::list<pair_type> matches;
  pair_type prev_range = EMPTY_PAIR;
  for(usint j = 0; j < this->number_of_sequences; j++)
  {
    // Encode the padding as a single run.
    pair_type seq_range = this->getSequenceRange(j);
    if(j > 0 && prev_range.second + 1 < seq_range.first)
    {
      this->encodePLCPRun(plcp, prev_range.second + 1, seq_range.first, seq_range.first - prev_range.second - 2);
    }
    prev_range = seq_range;

    usint maximal = seq_range.first;
    usint x = this->sa_samples->inverseSA(seq_range.first).second, next_x;

    // Invariant: x == inverseSA(i)
    for(usint i = seq_range.first; i <= seq_range.second; i++, x = next_x)
    {
      usint c = this->getCharacter(x);

      // T[i,n] is lexicographically the first suffix beginning with c.
      if(x == this->alphabet->cumulative(c))
      {
        if(!matches.empty())
        {
          maximal = (*(matches.begin())).first;
          matches.clear();
        }
        this->encodePLCPRun(plcp, maximal, i + 1, i - maximal);
        maximal = i + 1;
        next_x = this->psiUnsafe(x, c, *(iters[c])) - this->number_of_sequences;
        continue;
      }

      // Process the previous left matches we are still following.
      usint low = this->alphabet->cumulative(c) + this->number_of_sequences;
      usint start, stop; start = stop = maximal;
      for(std::list<pair_type>::iterator iter = matches.begin(); iter != matches.end();)
      {
        pair_type match = *iter;
        if(match.second < low)  // These no longer match the current character.
        {
          if(match.first < start) { start = match.first; }  // Form a single run from then.
          iter = matches.erase(iter);
        }
        else
        {
          if(match.first < stop) { stop = match.first; }  // End of the run to be encoded.
          (*iter).second = this->psiUnsafe(match.second - this->number_of_sequences, c, *(iters[c]));
          ++iter;
        }
      }
      if(start < stop) { this->encodePLCPRun(plcp, start, stop, i - start); }

      // If PLCP[i] is minimal, we add the left match to the list.
      // We add pairs of type (j, inverseSA(j) + n_of_s) as inverseSA is < 0 for end markers.
      usint next_y = this->psiUnsafe(x - 1, c, *(iters[c]));
      next_x = this->psiUnsafeNext(*(iters[c])) - this->number_of_sequences;
      if(next_y != next_x + this->number_of_sequences - 1)
      {
        matches.push_back(pair_type(maximal, next_y));
        maximal = i + 1;
      }
    }

    if(!matches.empty())
    {
      maximal = (*(matches.begin())).first;
      matches.clear();
    }
    if(maximal <= seq_range.second)
    {
      this->encodePLCPRun(plcp, maximal, seq_range.second + 1, seq_range.second + 1 - maximal);
    }
  }

  this->deleteIterators(iters);

  plcp.flush();
  return new PLCPVector(plcp, 2 * (this->end_points->getSize() + 1));
}



/*PLCPVector*
RLCSA::buildPLCP(usint block_size, usint threads) const
{
  if(block_size < 2 * sizeof(usint) || block_size % sizeof(usint) != 0)
  {
    std::cerr << "PLCP: Block size must be a multiple of " << sizeof(usint) << " bytes!" << std::endl;
    return 0;
  }

  #ifdef MULTITHREAD_SUPPORT
  threads = std::max(threads, (usint)1);
  threads = std::min(threads, this->getNumberOfSequences());
  omp_set_num_threads(threads);
  #else
  threads = 1;
  #endif

  PLCPVector::Encoder* plcp_encoders[threads];
  for(usint i = 0; i < threads; i++) { plcp_encoders[i] = new PLCPVector::Encoder(block_size); }

  #pragma omp parallel for schedule(dynamic, 1)
  for(usint j = 0; j < this->number_of_sequences; j++)
  {
    PsiVector::Iterator** iters = this->getIterators();
    std::list<pair_type> matches;
    PLCPVector::Encoder& plcp = *(plcp_encoders[omp_get_thread_num()]);

    // Encode the padding as a single run.
    pair_type seq_range = this->getSequenceRange(j);
    pair_type prev_range = (j == 0 ? EMPTY_PAIR : this->getSequenceRange(j - 1));
    if(j > 0 && prev_range.second + 1 < seq_range.first)
    {
      this->encodePLCPRun(plcp, prev_range.second + 1, seq_range.first, seq_range.first - prev_range.second - 2);
    }

    usint maximal = seq_range.first;
    usint x = this->sa_samples->inverseSA(seq_range.first).second, next_x;

    // Invariant: x == inverseSA(i)
    for(usint i = seq_range.first; i <= seq_range.second; i++, x = next_x)
    {
      usint c = this->getCharacter(x);

      // T[i,n] is lexicographically the first suffix beginning with c.
      if(x == this->alphabet->cumulative(c))
      {
        if(!matches.empty())
        {
          maximal = (*(matches.begin())).first;
          matches.clear();
        }
        this->encodePLCPRun(plcp, maximal, i + 1, i - maximal);
        maximal = i + 1;
        next_x = this->psiUnsafe(x, c, *(iters[c])) - this->number_of_sequences;
        continue;
      }

      // Process the previous left matches we are still following.
      usint low = this->alphabet->cumulative(c) + this->number_of_sequences;
      usint start, stop; start = stop = maximal;
      for(std::list<pair_type>::iterator iter = matches.begin(); iter != matches.end();)
      {
        pair_type match = *iter;
        if(match.second < low)  // These no longer match the current character.
        {
          if(match.first < start) { start = match.first; }  // Form a single run from then.
          iter = matches.erase(iter);
        }
        else
        {
          if(match.first < stop) { stop = match.first; }  // End of the run to be encoded.
          (*iter).second = this->psiUnsafe(match.second - this->number_of_sequences, c, *(iters[c]));
          ++iter;
        }
      }
      if(start < stop) { this->encodePLCPRun(plcp, start, stop, i - start); }

      // If PLCP[i] is minimal, we add the left match to the list.
      // We add pairs of type (j, inverseSA(j) + n_of_s) as inverseSA is < 0 for end markers.
      usint next_y = this->psiUnsafe(x - 1, c, *(iters[c]));
      next_x = this->psiUnsafeNext(*(iters[c])) - this->number_of_sequences;
      if(next_y != next_x + this->number_of_sequences - 1)
      {
        matches.push_back(pair_type(maximal, next_y));
        maximal = i + 1;
      }
    }

    if(!matches.empty())
    {
      maximal = (*(matches.begin())).first;
      matches.clear();
    }
    if(maximal <= seq_range.second)
    {
      this->encodePLCPRun(plcp, maximal, seq_range.second + 1, seq_range.second + 1 - maximal);
    }
    this->deleteIterators(iters);
  }


  // Produce vectors from each of the threads and merge them to get the final vector.
  PLCPVector* temp_vecs[threads];
  for(usint i = 0; i < threads; i++)
  {
    plcp_encoders[i]->flush();
    temp_vecs[i] = new PLCPVector(*(plcp_encoders[i]), 2 * (this->end_points->getSize() + 1));
    std::cout << "Vector " << i << ": " << temp_vecs[i]->getNumberOfItems() << " / " << temp_vecs[i]->getSize() << std::endl;
    delete plcp_encoders[i]; plcp_encoders[i] = 0;
  }
  if(threads == 1) { return temp_vecs[0]; }

  PLCPVector::Encoder plcp(block_size);
  PLCPVector::Iterator* iters[threads];
  pair_type next_run[threads];
  usint min_vec = 0;
  usint m = 0;
  for(usint i = 0; i < threads; i++)
  {
    iters[i] = new PLCPVector::Iterator(*(temp_vecs[i]));
    next_run[i] = iters[i]->selectRun(0, this->data_size + this->number_of_sequences); m += next_run[i].second + 1;
    if(next_run[i].first < next_run[min_vec].first) { min_vec = i; }
  }
  usint n = 0;
  while(next_run[min_vec].first < this->data_size + this->number_of_sequences)
  {
    plcp.addRun(next_run[min_vec].first, next_run[min_vec].second + 1); n += next_run[min_vec].second + 1;
    next_run[min_vec] = iters[min_vec]->selectNextRun(this->data_size + this->number_of_sequences); m += next_run[min_vec].second + 1;
    for(usint i = 0; i < threads; i++)
    {
      if(next_run[i].first < next_run[min_vec].first) { min_vec = i; }
    }
  }
  std::cout << n << " / " << m << std::endl;

  for(usint i = 0; i < threads; i++)
  {
    delete iters[i]; iters[i] = 0;
    delete temp_vecs[i]; temp_vecs[i] = 0;
  }

  plcp.flush();
  return new PLCPVector(plcp, 2 * (this->end_points->getSize() + 1));
}*/

usint
RLCSA::sampleLCP(usint sample_rate, pair_type*& sampled_values, bool report) const
{
  if(sample_rate == 0)
  {
    sample_rate = this->data_size + 1;
  }

  PsiVector::Iterator** iters = this->getIterators();

  usint runs = this->countRuns();
  usint samples = 0, minimal_sum = 0, nonstrict_minimal_sum = 0, minimal_samples = 0;
  usint max_samples = runs + (this->data_size - runs) / sample_rate;
  sampled_values = new pair_type[max_samples];

  std::list<Triple> matches;
  for(usint j = 0; j < this->number_of_sequences; j++)
  {
    pair_type seq_range = this->getSequenceRange(j);
    usint first_sample = samples; // First minimal sample of the current sequence.
    usint i, x = this->sa_samples->inverseSA(seq_range.first).second, next_x;

    // Invariant: x == inverseSA(i)
    for(i = seq_range.first; i <= seq_range.second; i++, x = next_x)
    {
      usint c = this->getCharacter(x);

      // T[i,n] is lexicographically the first suffix beginning with c.
      if(x == this->alphabet->cumulative(c))
      {
        if(!matches.empty())
        {
          for(std::list<Triple>::iterator iter = matches.begin(); iter != matches.end(); ++iter)
          {
            nonstrict_minimal_sum += i - (*iter).first;
          }
          matches.clear();
        }

        // Minimal sample: SA[x] = i, LCP[x] = 0
        sampled_values[samples] = pair_type(x, 0);
        samples++; minimal_samples++;
        next_x = this->psiUnsafe(x, c, *(iters[c])) - this->number_of_sequences;
        continue;
      }

      // Process the previous left matches we are still following.
      usint low = this->alphabet->cumulative(c) + this->number_of_sequences;
      pair_type potential_sample(this->data_size, this->data_size);
      for(std::list<Triple>::iterator iter = matches.begin(); iter != matches.end();)
      {
        Triple match = *iter;
        if(match.second < low)  // These no longer match the current character.
        {
          if(potential_sample.first != this->data_size) { nonstrict_minimal_sum += potential_sample.second; }

          // Potential minimal sample: SA[match.third] = match.first, LCP[match.third] = i - match.first
          potential_sample = pair_type(match.third, i - match.first);
          iter = matches.erase(iter);
        }
        else
        {
          (*iter).second = this->psiUnsafe(match.second - this->number_of_sequences, c, *(iters[c]));
          ++iter;
        }
      }
      if(potential_sample.first != this->data_size)
      {
        // Last potential sample is minimal.
        sampled_values[samples] = potential_sample;
        samples++; minimal_sum += potential_sample.second; minimal_samples++;
      }

      // If PLCP[i] is minimal, we add the left match to the list.
      // We add triples of type (j, inverseSA(j) + n_of_s, x) as inverseSA is < 0 for end markers.
      // x is the original minimal position.
      usint next_y = this->psiUnsafe(x - 1, c, *(iters[c]));
      next_x = this->psiUnsafeNext(*(iters[c]));
      if(next_y != next_x - 1 || next_x < this->number_of_sequences)
      {
        matches.push_back(Triple(i, next_y, x));
      }
      next_x -= this->number_of_sequences;
    }

    // Are we still following something?
    if(!matches.empty())
    {
      pair_type potential_sample(this->data_size, this->data_size);
      for(std::list<Triple>::iterator iter = matches.begin(); iter != matches.end(); ++iter)
      {
        Triple match = *iter;
        if(potential_sample.first != this->data_size) { nonstrict_minimal_sum += potential_sample.second; }

        // Potential minimal sample: SA[match.third] = match.first, LCP[match.third] = i - match.first
        potential_sample = pair_type(match.third, i - match.first);
      }
      matches.clear();
      if(potential_sample.first != this->data_size)
      {
        // Last potential sample is minimal.
        sampled_values[samples] = potential_sample;
        samples++; minimal_sum += potential_sample.second; minimal_samples++;
      }
    }

    // Add the non-minimal samples.
    if(sample_rate <= this->data_size)
    {
      usint last_sample = samples - 1;
      i = seq_range.first; x = this->sa_samples->inverseSA(seq_range.first).second;
      for(usint current_sample = first_sample; current_sample <= last_sample; current_sample++)
      {
        // Find the next minimal sample and add nonminimal samples if needed.
        pair_type first_nonminimal(i + sample_rate - 1, samples);
        pair_type next_nonminimal = first_nonminimal;
        while(x != sampled_values[current_sample].first)
        {
          if(i == next_nonminimal.first)
          {
            sampled_values[samples] = pair_type(x, 0); samples++;
            next_nonminimal.first += sample_rate; next_nonminimal.second++;
          }
          i++; x = this->psi(x, iters) - this->number_of_sequences;
        }

        // Reduce the nonminimal samples to the current minimal sample.
        for(next_nonminimal = first_nonminimal; next_nonminimal.second < samples; next_nonminimal.first += sample_rate, next_nonminimal.second++)
        {
          sampled_values[next_nonminimal.second].second = sampled_values[current_sample].second + i - next_nonminimal.first;
        }
        i++; x = this->psi(x, iters) - this->number_of_sequences;
      }
    }
  }

  std::sort(sampled_values, sampled_values + samples);

  if(report)
  {
    std::cout << "Samples: " << samples << " (total) / " << minimal_samples << " (minimal)" << std::endl;
    std::cout << "Upper bounds: " << max_samples << " (total) / " << runs << " (minimal)" << std::endl;
    std::cout << "Sum of minimal samples: " << (minimal_sum + nonstrict_minimal_sum) << " (total) / " << minimal_sum << " (strict)" << std::endl;
    std::cout << std::endl;
  }

  this->deleteIterators(iters);

  return samples;
}

usint
RLCSA::lcp(usint index, const LCPSamples& lcp_samples, bool steps) const
{
  if(index >= this->data_size) { return 0; }

  usint offset = 0;
  while(true)
  {
    if(lcp_samples.isSampled(index))
    {
      return (steps ? offset : lcp_samples.getSampleAt(index) + offset);
    }
    index = this->psi(index); this->convertToSAIndex(index);
    offset++;
  }
}

usint
RLCSA::lcpDirect(usint index) const
{
  if(index >= this->data_size || index == 0)
  {
    return 0;
  }

  usint match = index - 1;
  usint value = 0;

  usint c = this->getCharacter(index);
  usint low = this->alphabet->cumulative(c);
  while(match >= low)
  {
    PsiVector::Iterator iter(*(this->array[c]));
    match = this->psiUnsafe(match, c, iter);
    index = this->psiUnsafe(index, c, iter);
    value++;

    if(match < this->number_of_sequences || index < this->number_of_sequences)
    {
      break;
    }
    this->convertToSAIndex(match); this->convertToSAIndex(index);

    c = this->getCharacter(index);
    low = this->alphabet->cumulative(c);
  }

  return value;
}

//--------------------------------------------------------------------------

void
RLCSA::mergeEndPoints(RLCSA& index, RLCSA& increment)
{
  DeltaEncoder* endings = new DeltaEncoder(RLCSA::ENDPOINT_BLOCK_SIZE);

  DeltaVector::Iterator index_iter(*(index.end_points));
  DeltaVector::Iterator increment_iter(*(increment.end_points));

  endings->setBit(index_iter.select(0));
  for(usint i = 1; i < index.number_of_sequences; i++)
  {
    endings->setBit(index_iter.selectNext());
  }
  usint sum = index.end_points->getSize();
  delete index.end_points; index.end_points = 0;

  endings->setBit(sum + increment_iter.select(0));
  for(usint i = 1; i < increment.number_of_sequences; i++)
  {
    endings->setBit(sum + increment_iter.selectNext());
  }
  sum += increment.end_points->getSize();
  delete increment.end_points; increment.end_points = 0;

  this->end_points = new DeltaVector(*endings, sum);
  delete endings;
}


void
RLCSA::mergeSamples(RLCSA& index, RLCSA& increment, usint* positions)
{
  if(index.sa_samples == 0 || increment.sa_samples == 0) { return; }

  positions += increment.number_of_sequences;
  this->sa_samples = new SASamples(*(index.sa_samples), *(increment.sa_samples), positions, increment.data_size, this->number_of_sequences);

  this->support_locate = this->sa_samples->supportsLocate();
  this->support_display = this->sa_samples->supportsDisplay();
}

//--------------------------------------------------------------------------

void
RLCSA::strip()
{
  for(usint c = 0; c < CHARS; c++)
  {
    if(this->array[c] != 0) { this->array[c]->strip(); }
  }
  if(this->sa_samples != 0) { this->sa_samples->strip(); }
  this->end_points->strip();
}

//--------------------------------------------------------------------------

void
RLCSA::buildRLCSA(uchar* data, usint* ranks, usint bytes, usint block_size, usint threads, Sampler* sampler, bool multiple_sequences, bool delete_data)
{
  threads = std::max(threads, (usint)1);
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif


  // Do we store SA samples?
  bool sample_sa = true;
  if(this->sample_rate == 0)
  {
    sample_sa = false; this->sample_rate = 1;
  }


  // Determine the number of sequences and mark their end points.
  DeltaEncoder endings(RLCSA::ENDPOINT_BLOCK_SIZE);
  if(multiple_sequences)
  {
    this->number_of_sequences = 0;
    usint marker = 0;
    usint padding = 0, chars_encountered = 0;

    for(usint i = 0; i < bytes; i++)
    {
      if(data[i] == 0)
      {
        if(i == marker) { break; }  // Empty sequence.
        this->number_of_sequences++;
        marker = i + 1;
        usint pos = chars_encountered + padding - 1;
        endings.setBit(pos);
        padding = ((pos + this->sample_rate) / this->sample_rate) * this->sample_rate - chars_encountered;
      }
      else { chars_encountered++; }
    }

    if(this->number_of_sequences == 0 || marker != bytes)
    {
      std::cerr << "RLCSA: Collection must consist of 0-terminated nonempty sequences!" << std::endl;
      if(delete_data) { delete[] data; }
      return;
    }
    this->end_points = new DeltaVector(endings, chars_encountered + padding);
  }
  else
  {
    this->number_of_sequences = 1;
    DeltaEncoder endings(RLCSA::ENDPOINT_BLOCK_SIZE, RLCSA::ENDPOINT_BLOCK_SIZE);
    endings.setBit(bytes - 1);
    this->end_points = new DeltaVector(endings, bytes);
  }


  // Build character tables etc.
  usint distribution[CHARS];
  for(usint c = 0; c < CHARS; c++) { distribution[c] = 0; }
  for(usint i = 0; i < bytes; i++) { distribution[(usint)data[i]]++; }
  if(multiple_sequences) { distribution[0] = 0; } // \0 is an end marker
  this->alphabet = new Alphabet(distribution); this->data_size = this->alphabet->getDataSize();


  // Build suffix array.
  short_pair* sa = 0;
  if(multiple_sequences)
  {
    if(ranks == 0) { sa = simpleSuffixSort(data, bytes, this->number_of_sequences, threads); }
    else           { sa = simpleSuffixSort(ranks, bytes, threads); }
  }
  else
  {
    bytes++;
    sa = new short_pair[bytes];
    #pragma omp parallel for schedule(static)
    for(usint i = 0; i < bytes - 1; i++) { sa[i].second = data[i] + 1; }
    sa[bytes - 1].second = 0;
    simpleSuffixSort(sa, bytes, threads);
  }
  if(delete_data) { delete[] data; }


  // Sample SA.
  if(sample_sa)
  {
    if(sampler != 0) { this->sa_samples = new SASamples(sa, sampler, threads); }
    else             { this->sa_samples = new SASamples(sa, this->end_points, this->data_size, this->sample_rate, threads); }
    this->support_locate = this->sa_samples->supportsLocate();
    this->support_display = this->sa_samples->supportsDisplay();
  }


  // Build Psi.
  #pragma omp parallel for schedule(static)
  for(usint i = 0; i < bytes; i++) { sa[i].first = sa[(sa[i].first + 1) % bytes].second; }


  // Build RLCSA.
  #pragma omp parallel for schedule(dynamic, 1)
  for(usint c = 0; c < CHARS; c++)
  {
    if(!(this->alphabet->hasChar(c))) { this->array[c] = 0; continue; }

    short_pair* curr = sa + this->alphabet->cumulative(c) + this->number_of_sequences;
    short_pair* limit = curr + this->alphabet->countOf(c);
    PsiVector::Encoder encoder(block_size);
    pair_type run((*curr).first, 1); ++curr;

    for(; curr < limit; ++curr)
    {
      if((*curr).first == run.first + run.second) { run.second++; }
      else
      {
        encoder.addRun(run.first, run.second);
        run = pair_type((*curr).first, 1);
      }
    }
    encoder.addRun(run.first, run.second);
    encoder.flush();

    this->array[c] = new PsiVector(encoder, this->data_size + this->number_of_sequences);
  }
  delete[] sa;


  this->ok = true;
}

//--------------------------------------------------------------------------

} // namespace CSA
