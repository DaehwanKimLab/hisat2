#include <cstdlib>
#include <fstream>
#include <iostream>

#include "rlcsa.h"


using namespace CSA;


/*
  This program writes run-length encoded PLCP of the collection into a file.
*/


usint
countRuns(usint prev, uchar* buffer, usint length)
{
  usint runs = 0;

  for(usint i = 0; i < length; i++)
  {
    if(buffer[i] != prev)
    {
      prev = buffer[i];
      if(buffer[i] != 0) { runs++; }
    }
  }

  return runs;
}


int
main(int argc, char** argv)
{
  std::cout << "RLCSA to BWT converter" << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage: read_bwt base_name [buffer_size]" << std::endl;
    return 1;
  }

  std::string base_name = argv[1];
  std::string bwt_name = base_name + ".bwt";
  std::cout << "BWT: " << bwt_name << std::endl;
  std::ofstream bwt_file(bwt_name.c_str(), std::ios_base::binary);
  if(!bwt_file)
  {
    std::cerr << "Error creating BWT file!" << std::endl;
    return 2;
  }
  std::cout << std::endl;

  RLCSA rlcsa(base_name);
  clock_t start = clock();
  usint buffer_size = 0;
  if(argc > 2) { buffer_size = atoi(argv[2]); }
  usint n = rlcsa.getSize() + rlcsa.getNumberOfSequences();

  usint runs = 0, prev = CHARS;
  if(buffer_size > 0)
  {
    for(usint i = 0; i < n; i += buffer_size)
    {
      pair_type range(i, std::min(i + buffer_size - 1, n - 1));
      uchar* bwt = rlcsa.readBWT(range);
      if(bwt != 0)
      {
        runs += countRuns(prev, bwt, length(range));
        prev = bwt[length(range) - 1];
        bwt_file.write((char*)bwt, length(range));
        delete[] bwt;
      }
    }
  }
  else
  {
    uchar* bwt = rlcsa.readBWT();
    if(bwt != 0)
    {
      runs = countRuns(prev, bwt, n);
      bwt_file.write((char*)bwt, n);
      delete[] bwt;
    }
  }

  clock_t stop = clock();
  double time = ((stop - start) / (double)CLOCKS_PER_SEC);
  double megabytes = n / (double)MEGABYTE;
  std::cout << megabytes << " megabytes in " << time << " seconds (" << (megabytes / time) << " MB/s)" << std::endl;
  std::cout << std::endl;

  // Testing direct reporting of the number of runs.
  // This is as expensive as reading the BWT.
  std::cout << "Number of runs: " << runs << std::endl;
  runs = rlcsa.countRuns();
  std::cout << "Number of runs (direct count): " << runs << std::endl;
  std::cout << std::endl;

  bwt_file.close();
  return 0;
}
