#include <cstdlib>
#include <fstream>
#include <iostream>

#include "misc/utils.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
const int MAX_THREADS = 64;
#endif


using namespace CSA;


int main(int argc, char** argv)
{
  std::cout << "Parallel suffix sort test" << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage: ss_test base_name [threads]" << std::endl;
    return 1;
  }

  std::cout << "Base name: " << argv[1] << std::endl;
  uint threads = 1;
  if(argc >= 3) { threads = std::min(std::max(atoi(argv[2]), 1), MAX_THREADS); }
  std::cout << "Threads: " << threads << std::endl;
  std::cout << std::endl;

  std::ifstream input(argv[1], std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error: Cannot open input file " << argv[1] << std::endl;
    return 2;
  }
  uint n = fileSize(input);
  uchar* data = new uchar[n];
  input.read((char*)data, n);
  input.close();

  n++;
  usint* sequence = new usint[n];
  for(uint i = 0; i < n - 1; i++) { sequence[i] = data[i] + 1; }
  sequence[n - 1] = 0;
  delete[] data; data = 0;

  double start = readTimer();
  short_pair* result = simpleSuffixSort(sequence, n, threads);
  delete[] result; result = 0;
  double time = readTimer() - start;

  double megabytes = (n - 1) / (double)MEGABYTE;
  std::cout << megabytes << " megabytes in " << time << " seconds (" << (megabytes / time) << " MB/s)" << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return 0;
}
