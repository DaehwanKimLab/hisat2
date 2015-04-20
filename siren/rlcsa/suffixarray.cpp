#include <fstream>
#include <iostream>

#include "misc/utils.h"
#include "suffixarray.h"


namespace CSA
{

//--------------------------------------------------------------------------

SuffixArray::SuffixArray(const std::string& base_name, bool print) :
  ok(false),
  data(0), sa(0), original_sa(0), ranks(0), data_size(0),
  sequences(0)
{
  std::ifstream data_file(base_name.c_str(), std::ios_base::binary);
  if(!data_file)
  {
    std::cerr << "Error: Cannot open data file " << base_name << std::endl;
    return;
  }
  this->data_size = fileSize(data_file);
  this->data = new uchar[this->data_size];
  data_file.read((char*)(this->data), data_size);
  data_file.close();

  std::string sa_name = base_name + SA_EXTENSION;
  std::ifstream sa_file(sa_name.c_str(), std::ios_base::binary);
  if(!sa_file)
  {
    std::cerr << "Error: Cannot open suffix array file " << sa_name << std::endl;
    return;
  }
  sa_file.read((char*)&(this->sequences), sizeof(uint));
  this->sa = new uint[this->data_size];
  sa_file.read((char*)(this->sa), this->data_size * sizeof(uint));
  sa_file.close();

  this->original_sa = this->sa; this->sa += this->sequences;
  this->ok = true;
}

SuffixArray::SuffixArray(uchar* _data, uint bytes, uint threads) :
  ok(false),
  data(_data), sa(0), original_sa(0), ranks(0), data_size(bytes),
  sequences(0)
{
  if(_data == 0 || bytes == 0)
  {
    std::cerr << "Error: No input data given for suffix array construction!" << std::endl;
    return;
  }

  for(uint i = 0; i < this->data_size; i++) { if(this->data[i] == '\0') { this->sequences++; } }
  if(this->data[this->data_size - 1] != '\0')
  {
    std::cerr << "Error: Input data must end with \\0!" << std::endl;
    return;
  }

  short_pair* pairs = simpleSuffixSort(this->data, this->data_size, this->sequences, threads);
  this->sa = new uint[this->data_size];
  for(uint i = 0; i < this->data_size; i++) { this->sa[i] = pairs[i].first; }
  delete[] pairs;

  this->original_sa = this->sa; this->sa += this->sequences;
  this->ok = true;
}

SuffixArray::SuffixArray(uchar* _data, usint* ra, uint bytes, uint threads) :
  ok(false),
  data(_data), sa(0), original_sa(0), ranks(ra), data_size(bytes),
  sequences(0)
{
  if(_data == 0 || ra == 0 || bytes == 0)
  {
    std::cerr << "Error: No input data given for suffix array construction!" << std::endl;
    return;
  }

  for(uint i = 0; i < this->data_size; i++) { if(this->data[i] == '\0') { this->sequences++; } }
  if(this->data[this->data_size - 1] != '\0')
  {
    std::cerr << "Error: Input data must end with \\0!" << std::endl;
    return;
  }

  pair_type* pairs = new pair_type[this->data_size];
  for(usint i = 0; i < this->data_size; i++) { pairs[i] = pair_type(this->ranks[i], i); }
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif
  parallelSort(pairs, pairs + this->data_size);
  this->sa = new uint[this->data_size];
  for(uint i = 0; i < this->data_size; i++) { this->sa[i] = pairs[i].second; }
  delete[] pairs;

  this->original_sa = this->sa; this->sa += this->sequences;
  this->ok = true;
}

SuffixArray::~SuffixArray()
{
  delete[] this->data; this->data = 0;
  delete[] this->original_sa; this->sa = 0; this->original_sa = 0;
  delete[] this->ranks; this->ranks = 0;
  this->ok = false;
}

void
SuffixArray::writeTo(const std::string& base_name, bool write_data) const
{
  if(!this->ok) { return; }

  if(write_data)
  {
    std::ofstream data_file(base_name.c_str(), std::ios_base::binary);
    if(!data_file)
    {
      std::cerr << "Error: Cannot open data file " << base_name << std::endl;
      return;
    }
    data_file.write((char*)(this->data), this->data_size);
    data_file.close();
  }

  std::string sa_name = base_name + SA_EXTENSION;
  std::ofstream sa_file(sa_name.c_str(), std::ios_base::binary);
  if(!sa_file)
  {
    std::cerr << "Error: Cannot open suffix array file " << sa_name << std::endl;
    return;
  }
  sa_file.write((char*)&(this->sequences), sizeof(uint));
  sa_file.write((char*)(this->original_sa), this->data_size * sizeof(uint));
  sa_file.close();
}

//--------------------------------------------------------------------------

pair_type
SuffixArray::count(const std::string& pattern) const
{
  pair_type range;
  for(uint i = 0; i < pattern.length(); i++)
  {
    if(pattern[i] == '\0')
    {
      std::cerr << "Error: Pattern must not contain \\0!" << std::endl;
      return EMPTY_PAIR;
    }
  }

  // Lower bound for the range.
  sint low = 0, high = this->getSize() - 1;
  uint last_high = high;
  while(low < high)
  {
    uint mid = low + (high - low) / 2;
    uint matched = this->match(pattern, mid);
    if(matched >= pattern.length()) { high = mid; }
    else if((uchar)(pattern[matched]) < this->data[this->sa[mid] + matched])
    {
      high = ((sint)mid) - 1; last_high = high;
    }
    else { low = mid + 1; }
  }
  if(this->match(pattern, low) == pattern.length()) { range.first = low; }
  else { return EMPTY_PAIR; }

  // Upper bound for the range.
  high = last_high;
  while(low < high)
  {
    uint mid = low + (high - low + 1) / 2;
    uint matched = this->match(pattern, mid);
    if(matched >= pattern.length()) { low = mid; }
    else                            { high = mid - 1; }
  }
  range.second = high;

  return range;
}

uint
SuffixArray::match(const std::string& pattern, uint index) const
{
  uint matched = 0;
  for(uint pos = this->sa[index]; matched < pattern.length() && (uchar)(pattern[matched]) == this->data[pos]; matched++, pos++);
  return matched;
}

uint*
SuffixArray::locate(pair_type range) const
{
  if(isEmpty(range) || range.second >= this->getSize()) { return 0; }

  uint* result = new uint[length(range)];
  for(usint i = range.first; i <= range.second; i++) { result[i - range.first] = this->sa[i]; }
  return result;
}

uint
SuffixArray::locate(uint index) const
{
  if(index >= this->getSize()) { return this->data_size; }
  return this->sa[index];
}

//--------------------------------------------------------------------------

uint*
SuffixArray::getLCPArray(bool getPLCP) const
{
  if(!(this->isOk())) { return 0; }

  uint* plcp = new uint[this->data_size];
  for(uint i = 0; i < this->data_size; i++) { plcp[i] = this->data_size; }

  // Fill in the irreducible PLCP values.
  for(uint i = 0; i < this->sequences; i++)
  {
    plcp[this->original_sa[i]] = 0;
  }
  for(uint i = this->sequences; i < this->data_size; i++)
  {
    uint j = this->original_sa[i - 1], k = this->original_sa[i];
    if(j == 0 || k == 0 || this->data[j - 1] != this->data[k - 1] || this->data[j - 1] == 0)
    {
      uint temp = 0;
      while(this->data[j + temp] == this->data[k + temp] && this->data[j + temp] != 0) { temp++; }
      plcp[k] = temp;
    }
  }

  // Fill in the rest of the values.
  for(uint i = 1; i < this->data_size; i++)
  {
    if(plcp[i] >= this->data_size) { plcp[i] = plcp[i - 1] - 1; }
  }
  if(getPLCP) { return plcp; }

  // Convert to LCP array.
  uint* lcp = new uint[this->getSize()];
  for(uint i = 0; i < this->getSize(); i++)
  {
    lcp[i] = plcp[this->sa[i]];
  }
  delete[] plcp;

  return lcp;
}

//--------------------------------------------------------------------------

usint
SuffixArray::reportSize(bool print) const
{
  usint bytes = sizeof(*this) + this->data_size * (1 + sizeof(uint));
  if(this->ranks != 0) { bytes += this->data_size + sizeof(usint); }
  if(print)
  {
    std::cout << "Data:            " << (this->data_size / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Suffix array:    " << (this->data_size * sizeof(uint) / (double)MEGABYTE) << " MB" << std::endl;
    if(this->ranks != 0)
    {
      std::cout << "Rank array:      " << (this->data_size * sizeof(usint) / (double)MEGABYTE) << " MB" << std::endl;
    }
    std::cout << "Total size:      " << (bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << std::endl;
  }
  return bytes;
}

//--------------------------------------------------------------------------

} // namespace CSA
