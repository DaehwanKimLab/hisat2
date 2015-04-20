#include <iostream>

#include "suffixarray.h"
#include "misc/utils.h"


using namespace CSA;


int
main(int argc, char** argv)
{
  std::cout << "Suffix array builder" << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage: build_sa base_name" << std::endl;
    return 1;
  }

  std::string base_name = argv[1];
  std::cout << "Input: " << base_name << std::endl;
  std::cout << std::endl;

  std::ifstream input(base_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error: Cannot open input file!" << std::endl;
    return 2;
  }
  usint bytes = fileSize(input);
  uchar* data = new uchar[bytes];
  input.read((char*)data, bytes);
  input.close();

  double start = readTimer();
  SuffixArray sa(data, bytes);
  if(!sa.isOk()) { return 3; }
  double time = readTimer() - start;
  double megabytes = bytes / (double)MEGABYTE;
  std::cout << "Indexed " << megabytes << " megabytes in " << time << " seconds (" << (megabytes / time) << " MB/s)." << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  sa.writeTo(base_name);
  return 0;
}
