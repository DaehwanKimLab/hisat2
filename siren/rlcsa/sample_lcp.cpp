#include <algorithm>
#include <fstream>
#include <iostream>

#include "rlcsa.h"
#include "misc/utils.h"


using namespace CSA;


int
main(int argc, char** argv)
{
  std::cout << "LCP sampler" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: sample_lcp base_name sample_rate" << std::endl;
    return 1;
  }

  std::string base_name = argv[1];
  std::string lcp_name = base_name + LCP_SAMPLES_EXTENSION;
  std::cout << "LCP samples: " << lcp_name << std::endl;
  std::ofstream lcp_file(lcp_name.c_str(), std::ios_base::binary);
  if(!lcp_file)
  {
    std::cerr << "Error creating LCP sample file!" << std::endl;
    return 2;
  }

  usint sample_rate = atoi(argv[2]);
  std::cout << "Sample rate: " << sample_rate << std::endl;
  std::cout << std::endl;
  RLCSA* rlcsa = new RLCSA(base_name);

  double start = readTimer();
  pair_type* sampled_values = 0;
  usint samples = rlcsa->sampleLCP(sample_rate, sampled_values, true);
  usint data_size = rlcsa->getSize();
  delete rlcsa; // Saves memory.
  LCPSamples lcp(sampled_values, data_size, samples, true, true);
  lcp.writeTo(lcp_file);

  double time = readTimer() - start;
  double megabytes = data_size / (double)MEGABYTE;
  double size = lcp.reportSize() / (double)MEGABYTE;
  std::cout << megabytes << " megabytes in " << time << " seconds (" << (megabytes / time) << " MB/s)" << std::endl;
  std::cout << "Sampled LCP size: " << size << " MB" << std::endl;
  std::cout << std::endl;

  lcp_file.close();

  return 0;
}
