#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "rlcsa.h"
#include "misc/utils.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
const int MAX_THREADS = 64;
#endif


using namespace CSA;


int main(int argc, char** argv)
{
  std::cout << "RLCSA display test" << std::endl;
  if(argc < 4)
  {
    std::cout << "Usage: extract_sequence base_name sequence_number output" << std::endl;
    return 1;
  }

  std::cout << "Base name: " << argv[1] << std::endl;
  usint sequence = atoi(argv[2]);
  std::cout << "Sequence number: " << sequence << std::endl;
  std::cout << "Output: " << argv[3] << std::endl;
  std::cout << std::endl;

  RLCSA rlcsa(argv[1]);
  if(!rlcsa.supportsDisplay())
  {
    std::cerr << "Error: Display is not supported!" << std::endl;
    return 2;
  }
  rlcsa.printInfo();
  rlcsa.reportSize(true);
  if(sequence >= rlcsa.getNumberOfSequences())
  {
    std::cerr << "Error: Invalid sequence number!" << std::endl;
    return 3;
  }

  std::ofstream output(argv[3], std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error: Cannot open output file!" << std::endl;
    return 4;
  }

  double start = readTimer();
  uchar* buffer = rlcsa.display(sequence);
  usint bytes = length(rlcsa.getSequenceRange(sequence));
  output.write((char*)buffer, bytes);
  delete[] buffer;
  output.close();

  double time = readTimer() - start;
  double megabytes = bytes / (double)MEGABYTE;
  std::cout << megabytes << " megabytes in " << time << " seconds (" << (megabytes / time) << " MB/s)" << std::endl;
  std::cout << std::endl;

  return 0;
}
