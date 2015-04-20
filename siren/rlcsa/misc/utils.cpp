#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>

#include <sys/resource.h>

#include "utils.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif


namespace CSA
{

//--------------------------------------------------------------------------

std::streamoff
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

std::streamoff
fileSize(std::ofstream& file)
{
  std::streamoff curr = file.tellp();

  file.seekp(0, std::ios::end);
  std::streamoff size = file.tellp();
  file.seekp(0, std::ios::beg);
  size -= file.tellp();

  file.seekp(curr, std::ios::beg);
  return size;
}

void
largeWrite(std::ofstream& file, char* data, std::streamoff size, std::streamoff block_size)
{
	block_size = std::max(block_size, (std::streamoff)1);
	block_size *= MEGABYTE;

	for(std::streamoff i = 0; i < size; i += block_size)
	{
		file.write(data + i, std::min(size - i, block_size));
	}
}

//--------------------------------------------------------------------------

usint
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  usint chars = readRows(input, rows, skip_empty_rows);
  input.close();
  return chars;
}

usint
readRows(std::ifstream& file, std::vector<std::string>& rows, bool skip_empty_rows)
{
  usint chars = 0;
  while(file)
  {
    std::string buf;
    std::getline(file, buf);
    if(skip_empty_rows && buf.length() == 0) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }
  return chars;
}

usint
readPizzaChili(const std::string& filename, std::vector<std::string>& patterns)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readPizzaChili(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  usint chars = readPizzaChili(input, patterns);
  input.close();
  return chars;
}

usint
readPizzaChili(std::ifstream& file, std::vector<std::string>& patterns)
{
  std::string header;
  std::getline(file, header);

  std::size_t start = header.find("number=");
  if(start == std::string::npos) { return 0; }
  int temp = std::atoi(header.substr(start + 7).c_str());
  if(temp <= 0) { return 0; }
  usint n = temp;

  start = header.find("length=");
  if(start == std::string::npos) { return 0; }
  temp = std::atoi(header.substr(start + 7).c_str());
  if(temp <= 0) { return 0; }
  usint l = temp;

  char buffer[l];
  for(usint i = 0; i < n; i++)
  {
    file.read(buffer, l);
    patterns.push_back(std::string(buffer, l));
  }

  return n * l;
}

//--------------------------------------------------------------------------

double
readTimer()
{
  #ifdef MULTITHREAD_SUPPORT
  return omp_get_wtime();
  #else
  return clock() / (double)CLOCKS_PER_SEC;
  #endif
}

usint
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss;
}

//--------------------------------------------------------------------------

typedef std::pair<uint, uint> ss_range;
typedef std::pair<uint, usint> skew_pair;

struct SSSkewComparator
{
  bool operator() (const skew_pair& a, const skew_pair& b) const
  {
    return (a.second < b.second);
  }
} skew_comparator;

struct SSKeyComparator
{
  bool operator() (const short_pair& a, const short_pair& b) const
  {
    return (a.second < b.second);
  }
} key_comparator;

template<class T>
uint
setRanks(T* pairs, uint* keys, uint n, std::vector<ss_range>& unsorted, uint threads, uint chunk)
{
  uint total = 0;

  #ifdef MULTITHREAD_SUPPORT
  std::vector<ss_range> buffers[threads];
  uint subtotals[threads];
  for(uint i = 0; i < threads; i++) { subtotals[i] = 0; }
  #else
  std::vector<ss_range> next_unsorted;
  #endif

  #pragma omp parallel for schedule(dynamic, chunk)
  for(uint i = 0; i < unsorted.size(); i++)
  {
    #ifdef MULTITHREAD_SUPPORT
    uint thread = omp_get_thread_num();
    #endif
    T* prev = pairs + unsorted[i].first;
    keys[prev->first] = unsorted[i].first;
    T* limit = pairs + unsorted[i].second;
    while(prev < limit)
    {
      T* curr = prev + 1;
      if(curr->second != prev->second)
      {
        ++prev; keys[prev->first] = prev - pairs;
        continue;
      }
      keys[curr->first] = prev - pairs;

      for(++curr; curr <= limit && curr->second == prev->second; ++curr)
      {
        keys[curr->first] = prev - pairs;
      }

      #ifdef MULTITHREAD_SUPPORT
      buffers[thread].push_back(ss_range(prev - pairs, (curr - 1) - pairs));
      subtotals[thread] += curr - prev;
      #else
      next_unsorted.push_back(ss_range(prev - pairs, (curr - 1) - pairs));
      total += curr - prev;
      #endif
      prev = curr; if(prev <= limit) { keys[prev->first] = prev - pairs; }
    }
  }

  #ifdef MULTITHREAD_SUPPORT
  unsorted.clear();
  for(uint i = 0; i < threads; i++)
  {
    unsorted.insert(unsorted.end(), buffers[i].begin(), buffers[i].end());
    total += subtotals[i];
  }
  #else
  unsorted.assign(next_unsorted.begin(), next_unsorted.end());
  #endif

  return total;
}

uint
initialSort(short_pair* pairs, uint* keys, std::vector<ss_range>& unsorted, uint n, uint threads, uint h)
{
  // Sort according to first h characters.
  parallelSort(pairs, pairs + n, key_comparator);
  unsorted.push_back(ss_range(0, n - 1));
  uint total = setRanks(pairs, keys, n, unsorted, threads, 1);
//  std::cout << "Sorted with h = " << h << ", unsorted total = " << total << " (" << unsorted.size() << " ranges)" << std::endl;

  return total;
}

short_pair*
packPairs(skew_pair* old_pairs, uint n)
{
  short_pair* new_pairs = (short_pair*)old_pairs;
  for(uint i = 0; i < n; i++)
  {
    new_pairs[i].first = old_pairs[i].first;
    new_pairs[i].second = old_pairs[i].second;
  }
  return new_pairs;
}

short_pair*
prefixDoubling(short_pair* pairs, uint* keys, std::vector<ss_range>& unsorted, uint n, uint threads, uint total, uint h)
{
  // Double prefix length until sorted.
  while(total > 0)
  {
    uint chunk = std::max<size_t>((size_t)1, unsorted.size() / (threads * threads));
    #pragma omp parallel for schedule(dynamic, chunk)
    for(uint i = 0; i < unsorted.size(); i++)
    {
      // Set sort keys for the current range.
      short_pair* limit = pairs + unsorted[i].second;
      for(short_pair* curr = pairs + unsorted[i].first; curr <= limit; ++curr)
      {
        curr->second = keys[curr->first + h];
      }
      sequentialSort(pairs + unsorted[i].first, pairs + unsorted[i].second + 1, key_comparator);
    }
    total = setRanks(pairs, keys, n, unsorted, threads, chunk);
    h *= 2;
//    std::cout << "Sorted with h = " << h << ", unsorted total = " << total << " (" << unsorted.size() << " ranges)" << std::endl;
  }

  #pragma omp parallel for schedule(static)
  for(uint i = 0; i < n; i++) { pairs[i].second = keys[i]; }
  delete[] keys; keys = 0;
  return pairs;
}

short_pair*
prefixTripling(skew_pair* pairs, uint* keys, std::vector<ss_range>& unsorted, uint n, uint threads, uint total, uint h)
{
  const usint PACK_FACTOR = sizeof(usint) * CHAR_BIT / 2;

  // Triple prefix length until sorted.
  while(total > 0)
  {
    uint chunk = std::max<size_t>((size_t)1, unsorted.size() / (threads * threads));
    #pragma omp parallel for schedule(dynamic, chunk)
    for(uint i = 0; i < unsorted.size(); i++)
    {
      // Set sort keys for the current range.
      skew_pair* limit = pairs + unsorted[i].second;
      for(skew_pair* curr = pairs + unsorted[i].first; curr <= limit; ++curr)
      {
        curr->second = (usint)(keys[curr->first + h]) << PACK_FACTOR;
        if(n - h > curr->first + h) { curr->second += keys[curr->first + 2 * h]; }
      }
      sequentialSort(pairs + unsorted[i].first, pairs + unsorted[i].second + 1, skew_comparator);
    }
    total = setRanks(pairs, keys, n, unsorted, threads, chunk);
    h *= 3;
//    std::cout << "Sorted with h = " << h << ", unsorted total = " << total << " (" << unsorted.size() << " ranges)" << std::endl;
  }

  #pragma omp parallel for schedule(static)
  for(uint i = 0; i < n; i++) { pairs[i].second = keys[i]; }
  delete[] keys; keys = 0;
  return packPairs(pairs, n);
}

short_pair*
simpleSuffixSort(const usint* sequence, uint n, uint threads)
{
  if(sequence == 0 || n == 0) { return 0; }

  skew_pair* pairs = (skew_pair*)new short_pair[n * sizeof(skew_pair) / sizeof(short_pair) + 1];
  uint* keys = new uint[n];               // In text order.
  std::vector<ss_range> unsorted;
  threads = std::max(threads, (uint)1);
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif

  // Initialize pairs.
  #pragma omp parallel for schedule(static)
  for(uint i = 0; i < n; i++) { pairs[i].first = i; pairs[i].second = sequence[i]; }

  // Sort according to first character.
  parallelSort(pairs, pairs + n, skew_comparator);
  unsorted.push_back(ss_range(0, n - 1));
  uint total = setRanks(pairs, keys, n, unsorted, threads, 1);

  if(sizeof(usint) < 2 * sizeof(uint))
  {
    return prefixDoubling(packPairs(pairs, n), keys, unsorted, n, threads, total, 1);
  }
  else
  {
    return prefixTripling(pairs, keys, unsorted, n, threads, total, 1);
  }
}

short_pair*
simpleSuffixSort(short_pair* pairs, uint n, uint threads)
{
  if(pairs == 0 || n == 0) { return 0; }

  uint* keys = new uint[n];   // In text order.
  std::vector<ss_range> unsorted;
  threads = std::max(threads, (uint)1);
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif

  // Initialize pairs.
  #pragma omp parallel for schedule(static)
  for(uint i = 0; i < n; i++) { pairs[i].first = i; }

  uint total = initialSort(pairs, keys, unsorted, n, threads, 1);
  return prefixDoubling(pairs, keys, unsorted, n, threads, total, 1);
}

short_pair*
simpleSuffixSort(const uchar* sequence, uint n, uint sequences, uint threads)
{
  if(sequence == 0 || n == 0) { return 0; }

  short_pair* pairs = new short_pair[n];  // In sorted order.
  uint* keys = new uint[n];               // In text order.
  std::vector<ss_range> unsorted;
  threads = std::max(threads, (uint)1);
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif

  // Remap alphabet.
  uint alphabet[CHARS];
  for(uint c = 0; c < CHARS; c++) { alphabet[c] = 0; }
  for(uint i = 0; i < n; i++) { alphabet[sequence[i]]++; }
  uint alphabet_size = sequences;
  for(uint c = 1; c < CHARS; c++)
  {
    uint temp = alphabet_size;
    if(alphabet[c] > 0) { alphabet_size++; }
    alphabet[c] = temp;
  }

  // Determine pack factor.
  uint limit = std::numeric_limits<uint>::max() / alphabet_size;
  uint h = 1, pack_multiplier = 1;
  if(alphabet_size > 1)
  {
    while(pack_multiplier * alphabet_size <= limit) { h++; pack_multiplier *= alphabet_size; }
  }

  // Initialize pairs.
  uint zeros = 0, value = 0;
  for(uint i = 0; i < h; i++)
  {
    value *= alphabet_size;
    if(sequence[i] == 0) { value += zeros; zeros++; }
    else { value += alphabet[sequence[i]]; }
  }
  for(uint i = 0; i < n - h; i++)
  {
    pairs[i].first = i; pairs[i].second = value;
    value = (value % pack_multiplier) * alphabet_size;
    if(sequence[i + h] == 0) { value += zeros; zeros++; }
    else { value += alphabet[sequence[i + h]]; }
  }
  for(uint i = n - h; i < n; i++)
  {
    pairs[i].first = i; pairs[i].second = value;
    value = (value % pack_multiplier) * alphabet_size;
  }

  // Initialize pairs.
/*  uint zeros = 0;
  for(uint i = 0; i < n; i++)
  {
    pairs[i].first = i;
    pairs[i].second = (sequence[i] == '\0' ? zeros++ : sequence[i] + sequences);
  }
  uint h = 1;
  if(length(CHARS + sequences - 1) <= sizeof(uint) * CHAR_BIT / 2)
  {
    for(uint i = 0; i < n - 1; i++)
    {
      pairs[i].second = (CHARS + sequences) * pairs[i].second + pairs[i + 1].second;
    }
    pairs[n - 1].second *= CHARS + sequences;
    h = 2;
  }*/

  uint total = initialSort(pairs, keys, unsorted, n, threads, h);
  return prefixDoubling(pairs, keys, unsorted, n, threads, total, h);
}

//--------------------------------------------------------------------------

void
mergeRanges(std::vector<pair_type>* vec, bool parallel)
{
  if(vec == 0 || vec->size() <= 1) { return; }
  if(parallel) { parallelSort(vec->begin(), vec->end()); }
  else         { sequentialSort(vec->begin(), vec->end()); }

  std::vector<pair_type>::iterator prev = vec->begin();
  for(std::vector<pair_type>::iterator curr = prev + 1; curr != vec->end(); ++curr)
  {
    if(prev->second + 1 >= curr->first) { prev->second = std::max(curr->second, prev->second); }
    else { ++prev; *prev = *curr; }
  }
  vec->resize((prev - vec->begin()) + 1);
}

//--------------------------------------------------------------------------

} // namespace CSA
