#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "rlcsa_builder.h"


using namespace CSA;



const int MAX_THREADS = 64;


int
main(int argc, char** argv)
{
  std::cout << "RLCSA builder" << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage: build_rlcsa base_name [threads]" << std::endl;
    return 1;
  }

  std::string base_name = argv[1];
  int threads = 1;
  if(argc > 2)
  {
    threads = std::max(threads, atoi(argv[2]));
    threads = std::min(threads, MAX_THREADS);
  }
  std::cout << "Base name: " << base_name << std::endl;
  std::cout << "Threads: " << threads << std::endl;
  std::cout << std::endl;

  std::string parameters_name = base_name + PARAMETERS_EXTENSION;
  Parameters parameters;
  parameters.set(RLCSA_BLOCK_SIZE);
  parameters.set(SAMPLE_RATE);
  parameters.set(SUPPORT_LOCATE);
  parameters.set(SUPPORT_DISPLAY);
  parameters.set(WEIGHTED_SAMPLES);
  parameters.read(parameters_name);
  parameters.print();

  double start = readTimer();
  std::ifstream input(base_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error opening input file " << base_name << std::endl;
    return 2;
  }
  usint size = fileSize(input);
  if(size <= 0)
  {
    std::cerr << "Error: Input file is empty!" << std::endl;
    input.close();
    return 3;
  }
  uchar* data = new uchar[size];
  input.read((char*)data, size);
  input.close();

  RLCSA rlcsa(data, size, parameters.get(RLCSA_BLOCK_SIZE), parameters.get(SAMPLE_RATE), threads, true);
  data = 0; // RLCSA constructor deleted the data!
  if(!(rlcsa.isOk())) { return 4; }
  rlcsa.printInfo();
  rlcsa.reportSize(true);
  rlcsa.writeTo(base_name);
  double total = readTimer() - start;

  double megabytes = size / (double)MEGABYTE;
  std::cout << "Indexed " << megabytes << " megabytes in " << total << " seconds (" << (megabytes / total) << " MB/s)." << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return 0;
}
