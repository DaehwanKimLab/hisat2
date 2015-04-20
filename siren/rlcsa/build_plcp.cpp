#include <cstdlib>
#include <fstream>
#include <iostream>

#include "rlcsa.h"


using namespace CSA;


/*
  This program writes run-length encoded PLCP of the collection into a file.
*/


int
main(int argc, char** argv)
{
  std::cout << "PLCP builder" << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage: build_plcp base_name [block_size]" << std::endl;
    return 1;
  }

  std::string base_name = argv[1];
  std::string plcp_name = base_name + PLCP_EXTENSION;
  std::cout << "PLCP: " << plcp_name << std::endl;
  std::ofstream plcp_file(plcp_name.c_str(), std::ios_base::binary);
  if(!plcp_file)
  {
    std::cerr << "Error creating PLCP file!" << std::endl;
    return 2;
  }

  usint block_size = 32;
  if(argc > 2) { block_size = atoi(argv[2]); }
  std::cout << "Block size: " << block_size << std::endl;
  std::cout << std::endl;
  RLCSA rlcsa(base_name);

  clock_t start = clock();
  PLCPVector* plcp = rlcsa.buildPLCP(block_size);
  plcp->writeTo(plcp_file);
  clock_t stop = clock();

  double time = ((stop - start) / (double)CLOCKS_PER_SEC);
  double megabytes = rlcsa.getSize() / (double)MEGABYTE;
  double size = plcp->reportSize() / (double)MEGABYTE;
  PLCPVector::Iterator iter(*plcp);
  usint runs = iter.countRuns();

  std::cout << megabytes << " megabytes in " << time << " seconds (" << (megabytes / time) << " MB/s)" << std::endl;
  std::cout << "PLCP size: " << size << " MB" << std::endl;
  std::cout << "Runs: " << runs << std::endl;
  std::cout << std::endl;

  plcp_file.close();
  delete plcp;

  return 0;
}
