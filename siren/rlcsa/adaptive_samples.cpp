#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "adaptive_samples.h"
#include "misc/utils.h"



namespace CSA
{

//--------------------------------------------------------------------------

AdaptiveSamples::AdaptiveSamples(const RLCSA& rlcsa, const std::string& base_name) :
  index(rlcsa),
  samples(0), size(0), items(0),
  candidate_samples(0), c_items(0),
  regular_samples(0),
  text_start(0),
  promote_probability(0.0),
  pos(0), sum(0), window(0),
  ok(false)
{
  Parameters parameters;
  parameters.read(base_name + PARAMETERS_EXTENSION);
  this->use_candidates = parameters.get(CANDIDATE_SAMPLES);
  this->half_greedy = parameters.get(HALF_GREEDY_SAMPLES);
  this->window_size = parameters.get(SAMPLE_WINDOW_SIZE);

  if(parameters.get(SAMPLE_PROMOTE_RATE) > 0)
  {
    this->promote_probability = 1.0 / parameters.get(SAMPLE_PROMOTE_RATE);
  }

  std::string file_name = base_name + SA_SAMPLES_EXTENSION;
  std::ifstream sample_file(file_name.c_str(), std::ios_base::binary);
  if(!sample_file)
  {
    std::cerr << "AdaptiveSamples: Cannot open sample file!" << std::endl;
    return;
  }

  // Determine the number of adaptive, candidate, and regular samples.
  usint total_samples = 0;
  usint sample_rate = parameters.get(SAMPLE_RATE);
  sample_file.read((char*)&(total_samples), sizeof(total_samples)); // Actually text size.
  sample_file.read((char*)&(total_samples), sizeof(total_samples));
  this->size = (this->half_greedy ? total_samples / 2 : total_samples);
  if(this->use_candidates) { this->size /= 2; }

  // Initialize the samples.
  this->samples = new sample_type[this->size];
  for(key_type i = 0; i < this->size; i++)
  {
    this->samples[i] = sample_type(this->index.getSize(), this->index.getSize());
  }
  if(this->use_candidates)
  {
    this->candidate_samples = new sample_type[this->size];
    for(key_type i = 0; i < this->size; i++) { this->candidate_samples[i] = this->samples[i]; }
  }

  // Insert the initial samples.
  pair_type buffer;
  pair_type* regulars = 0;
  usint reg_samples = 0;
  if(this->half_greedy) { regulars = new pair_type[total_samples / 2]; }
  for(usint i = 0; i < total_samples; i++)
  {
    sample_file.read((char*)&buffer, sizeof(pair_type));
    if(this->half_greedy && buffer.second % (2 * sample_rate) == 0)
    {
      regulars[reg_samples++] = buffer;
    }
    else
    {
      this->addSample(sample_type(buffer.first, buffer.second));
    }
    if(buffer.second == 0) { this->text_start = buffer.first; }
  }
  if(this->half_greedy)
  {
    this->regular_samples = new SASamples(regulars, this->index.getSize(), 2 * sample_rate, 1);
    delete[] regulars; regulars = 0;
  }

  // Initialize the window.
  if(this->promote_probability > 0.0)
  {
    this->window = new key_type[this->window_size];
    for(usint i = 0; i < this->window_size; i++) { this->window[i] = sample_rate / 2; }
    this->sum = this->window_size * sample_rate / 2;
  }

  sample_file.close();
  this->ok = true;
}

AdaptiveSamples::~AdaptiveSamples()
{
  delete[] samples; samples = 0;
  delete[] candidate_samples; candidate_samples = 0;
  delete   regular_samples; regular_samples = 0;
  delete[] window; window = 0;
}

//--------------------------------------------------------------------------

usint
AdaptiveSamples::getNumberOfSamples() const
{
  usint temp = this->items + this->c_items;
  if(this->half_greedy) { temp += this->regular_samples->getNumberOfSamples(); }
  return temp;
}

double
AdaptiveSamples::getLoad() const
{
  double temp = this->items + this->c_items;
  temp /= (this->use_candidates ? 2 * this->size : this->size);
  return temp;
}

double
AdaptiveSamples::getPrimaryLoad() const
{
  return this->items / (double)(this->size);
}

double
AdaptiveSamples::getCandidateLoad() const
{
  return this->c_items / (double)(this->size);
}

usint
AdaptiveSamples::reportSize() const
{
  usint temp = sizeof(*this) + this->size * sizeof(sample_type);
  if(this->use_candidates) { temp += this->size * sizeof(sample_type); }
  if(this->half_greedy) { temp += this->regular_samples->reportSize(); }
  if(this->promote_probability > 0.0) { temp += this->window_size * sizeof(key_type); }
  return temp;
}

void
AdaptiveSamples::report() const
{
  std::cout << "Adaptive samples:" << std::endl;
  std::cout << "Size:     " << this->reportSize() / (double)MEGABYTE << " MB" << std::endl;
  std::cout << "Load:     " << this->getLoad() << " (" << this->getPrimaryLoad() << " / " << this->getCandidateLoad() << ")" << std::endl;

  std::cout << "Options: ";
  if(this->use_candidates) { std::cout << " candidates"; }
  if(this->half_greedy) { std::cout << " half-greedy"; }
  if(this->promote_probability > 0.0)
  {
    std::cout << " random(" << this->promote_probability << ", " << this->window_size << ")";
  }
  std::cout << std::endl;

  std::cout << std::endl;
}

//--------------------------------------------------------------------------

usint
AdaptiveSamples::locate(usint i, bool steps)
{
  if(i >= this->index.getSize()) { return (steps ? 0 : this->index.getSize()); }

  usint offset = 0;
  sample_type sample(i, 0);
  this->index.convertToBWTIndex(i);

  while(true)
  {
    // Priority 1: Implicit samples.
    if(this->index.hasImplicitSample(i))
    {
      sample.second = this->index.getImplicitSample(i) - offset;
      break;
    }
    this->index.convertToSAIndex(i);

    // Priority 2: Regular samples.
    if(this->half_greedy && this->regular_samples->isSampled(i))
    {
      sample.second = this->regular_samples->getSampleAt(i) - offset;
      break;
    }

    // Priority 3: Primary samples.
    key_type j = this->hash(i);
    if(this->samples[j].first == i)
    {
      sample.second = this->samples[j].second - offset;
      break;
    }

    // Priority 4: Candidate samples.
    if(this->use_candidates)
    {
      j = this->secondaryHash(i);
      if(this->candidate_samples[j].first == i)
      {
        sample.second = this->candidate_samples[j].second - offset;
        this->addSample(this->candidate_samples[j]);
        break;
      }
    }

    i = this->index.psi(i); offset++;
  }

  if(offset > 0) { this->trySample(sample, offset); }

  // Update the average number of steps taken.
  if(this->promote_probability > 0.0)
  {
    this->sum = this->sum + offset - this->window[this->pos];
    this->window[this->pos] = offset; this->pos = (this->pos + 1) % this->window_size;
  }

  return (steps ? offset : sample.second);
}

usint*
AdaptiveSamples::locate(pair_type range, bool steps)
{
  if(length(range) == 0) { return 0; }

  usint* data = new usint[length(range)];
  for(usint i = 0, j = range.first; j <= range.second; i++, j++)
  {
    data[i] = this->locate(j, steps);
  }

  return data;
}

void
AdaptiveSamples::display(pair_type range, uchar* data)
{
  if(!(this->supportsDisplay()) || isEmpty(range) || range.second >= this->index.getSize() || data == 0) { return; }

  pair_type res = this->regular_samples->inverseSA(range.first);
  usint i = res.first, pos = res.second;
  for(; i < range.first; i++)
  {
    pos = this->index.psi(pos);
    this->index.convertToSAIndex(pos);
  }
  for(; i <= range.second; i++)
  {
    data[i - range.first] = this->index.getCharacter(pos);
    pos = this->index.psi(pos);
    this->index.convertToSAIndex(pos);
  }
}

//--------------------------------------------------------------------------

void
AdaptiveSamples::addSample(sample_type sample)
{
  key_type i = this->hash(sample.first);
  if(this->samples[i].first >= this->index.getSize()) { this->items++; }
  else if(this->use_candidates) { this->addCandidate(this->samples[i]); }
  this->samples[i] = sample;
}

void
AdaptiveSamples::addCandidate(sample_type sample)
{
  key_type i = this->secondaryHash(sample.first);
  if(this->candidate_samples[i].first >= this->index.getSize()) { this->c_items++; }
  this->candidate_samples[i] = sample;
}

void
AdaptiveSamples::trySample(sample_type sample, usint offset)
{
  if(this->promote_probability > 0.0)
  {
    double temp = (this->promote_probability * this->window_size * offset) / this->sum;
    if(rand() / (RAND_MAX + 1.0) >= temp) { return; }
  }

  if(this->use_candidates) { this->addCandidate(sample); }
  else                     { this->addSample(sample); }
}

AdaptiveSamples::key_type
AdaptiveSamples::hash(key_type key)
{
  key = ~key + (key << 15);
  key = key ^ (key >> 12);
  key = key + (key << 2);
  key = key ^ (key >> 4);
  key = key * 2057;
  key = key ^ (key >> 16);
  return key % this->size;
}

AdaptiveSamples::key_type
AdaptiveSamples::secondaryHash(key_type key)
{
  key = (key + 0x7ed55d16) + (key << 12);
  key = (key ^ 0xc761c23c) ^ (key >> 19);
  key = (key + 0x165667b1) + (key << 5);
  key = (key + 0xd3a2646c) ^ (key << 9);
  key = (key + 0xfd7046c5) + (key << 3);
  key = (key ^ 0xb55a4f09) ^ (key >> 16);
  return key % this->size;
}

//--------------------------------------------------------------------------

} // namespace CSA
