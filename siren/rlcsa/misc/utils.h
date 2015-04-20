#ifndef _RLCSA_UTILS_H
#define _RLCSA_UTILS_H

#include <algorithm>
#include <fstream>
#include <vector>

#include "definitions.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif


namespace CSA
{



template<class T>
struct GenericTriple
{
  T first;
  T second;
  T third;

  GenericTriple() {}
  GenericTriple(T a, T b, T c) : first(a), second(b), third(c) {}
};

typedef GenericTriple<usint> Triple;


//--------------------------------------------------------------------------

std::streamoff fileSize(std::ifstream& file);
std::streamoff fileSize(std::ofstream& file);

// This is a workaround to a bug in ext4 file system. Writes larger than 2 gigabytes
// can fail, if done directly. Probably not relevant anymore.
// Block size is in megabytes.
void largeWrite(std::ofstream& file, char* data, std::streamoff size, std::streamoff block_size);

//--------------------------------------------------------------------------

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A, B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

// Returns the total length of the rows, excluding line ends.
usint readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);
usint readRows(std::ifstream& file, std::vector<std::string>& rows, bool skip_empty_rows);

// The same as above, but for the Pizza & Chili Corpus pattern format.
usint readPizzaChili(const std::string& filename, std::vector<std::string>& patterns);
usint readPizzaChili(std::ifstream& file, std::vector<std::string>& patterns);

//--------------------------------------------------------------------------

double readTimer();
usint  memoryUsage(); // Peak memory usage in kilobytes.

//--------------------------------------------------------------------------

// Input: sequence of length n. sequence[n-1] must be unique.
// Output: (SA[i], SA^-1[i])
typedef std::pair<uint, uint> short_pair;
short_pair* simpleSuffixSort(const usint* sequence, uint n, uint threads = 1);

// Input: sequence of pairs of length n.
//        character value is the second element of a pair.
//        character n-1 must be unique.
// Output: (SA[i], SA^-1[i])
short_pair* simpleSuffixSort(short_pair* pairs, uint n, uint threads = 1);

// Input: sequence of length n. sequence[n-1] must be '\0'.
// Output: (SA[i], SA^-1[i])
short_pair* simpleSuffixSort(const uchar* sequence, uint n, uint sequences, uint threads = 1);

//--------------------------------------------------------------------------

template <class Iterator>
void parallelSort(Iterator first, Iterator last)
{
  #ifdef MULTITHREAD_SUPPORT
    #ifdef _GLIBCXX_PARALLEL
      std::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
    #else
      std::sort(first, last);
    #endif
  #else
    std::sort(first, last);
  #endif
}

template <class Iterator, class Compare>
void
parallelSort(Iterator first, Iterator last, const Compare& comp)
{
  #ifdef MULTITHREAD_SUPPORT
    #ifdef _GLIBCXX_PARALLEL
      std::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
    #else
      std::sort(first, last, comp);
    #endif
  #else
    std::sort(first, last, comp);
  #endif
}

template <class Iterator>
void sequentialSort(Iterator first, Iterator last)
{
  #ifdef MULTITHREAD_SUPPORT
    #ifdef _GLIBCXX_PARALLEL
      std::sort(first, last, __gnu_parallel::sequential_tag());
    #else
      std::sort(first, last, mcstl::sequential_tag());
    #endif
  #else
    std::sort(first, last);
  #endif
}

template <class Iterator, class Compare>
void
sequentialSort(Iterator first, Iterator last, const Compare& comp)
{
  #ifdef MULTITHREAD_SUPPORT
    #ifdef _GLIBCXX_PARALLEL
      std::sort(first, last, comp, __gnu_parallel::sequential_tag());
    #else
      std::sort(first, last, comp, mcstl::sequential_tag());
    #endif
  #else
    std::sort(first, last, comp);
  #endif
}

template<class T>
void
removeDuplicates(std::vector<T>* vec, bool parallel = true)
{
  if(vec == 0) { return; }
  if(parallel) { parallelSort(vec->begin(), vec->end()); }
  else         { sequentialSort(vec->begin(), vec->end()); }
  vec->resize(std::unique(vec->begin(), vec->end()) - vec->begin());
}

template<class T>
void
removeDuplicates(std::vector<T>& vec, bool parallel = true)
{
  removeDuplicates(&vec, parallel);
}

void mergeRanges(std::vector<pair_type>* vec, bool parallel = true);

//--------------------------------------------------------------------------

} // namespace CSA


#endif // _RLCSA_UTILS_H
