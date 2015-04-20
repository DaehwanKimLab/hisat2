#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "rlcsa.h"
#include "misc/utils.h"


using namespace CSA;


int main(int argc, char** argv)
{
  std::cout << "RLCSA locate test" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: locate_test basename [begin end] output [d|s]" << std::endl;
    std::cout << "  d  Use direct locate." << std::endl;
    std::cout << "  s  Output a 32-bit suffix array." << std::endl;
    return 1;
  }

  usint begin = 0, end = 0;
  int output_arg = 2;
  std::cout << "Base name: " << argv[1] << std::endl;
  if(argc >= 5)
  {
    output_arg = 4;
    begin = atoi(argv[2]), end = atoi(argv[3]);
    std::cout << "Begin: " << begin << std::endl;
    std::cout << "End: " << end << std::endl;
    if(begin > end)
    {
      std::cerr << "Error: Empty range!" << std::endl;
      return 2;
    }
  }

  bool direct = false, do_output = false;
  for(int i = output_arg + 1; i < argc; i++)
  {
    switch(argv[i][0])
    {
      case 'd':
        direct = true; break;
      case 's':
        do_output = true; break;
    }
  }
  if(!direct)
  {
    std::cout << "Using run-based optimizations." << std::endl;
  }
  else
  {
    std::cout << "Using direct locate." << std::endl;
  }
  if(do_output)
  {
    std::cout << "Writing plain SA to output file, ignoring begin and end." << std::endl;
  }
  std::cout << std::endl;

  RLCSA rlcsa(argv[1]);
  if(!rlcsa.supportsLocate())
  {
    std::cerr << "Error: Locate is not supported!" << std::endl;
    return 3;
  }

  if(output_arg == 4 && !do_output)
  {
    if(end >= rlcsa.getSize())
    {
      std::cerr << "Error: Invalid range!" << std::endl;
      return 4;
    }
  }
  else
  {
    begin = 0; end = rlcsa.getSize() - 1;
  }

  std::ofstream output(argv[output_arg], std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error: Cannot open output file!" << std::endl;
    return 5;
  }

  usint* buffer = new usint[MILLION];
  usint temp;
  std::clock_t start = std::clock();
  if(do_output) { temp = rlcsa.getSize(); output.write((char*)&temp, sizeof(temp)); }
  for(usint curr = begin; curr <= end; curr += MILLION)
  {
    pair_type range(curr, std::min(end, curr + MILLION - 1));
    if(direct)
    {
      for(usint i = 0; i < range.second + 1 - range.first; i++)
      {
        buffer[i] = rlcsa.locate(curr + i, false);
      }
    }
    else { rlcsa.locate(range, buffer); }
    for(usint i = 0; i < range.second + 1 - range.first; i++)
    {
      if(do_output) { temp = buffer[i]; output.write((char*)&temp, sizeof(temp)); }
      else          { output.write((char*)&(buffer[i]), sizeof(usint)); }
    }
  }
  std::clock_t stop = std::clock();
  delete[] buffer;

  double size = (end + 1 - begin);
  double time = (stop - start) / (double)CLOCKS_PER_SEC;
  std::cout << size << " locates in " << time << " seconds (" << (size / time) << " locates/s)" << std::endl;
  output.close();

  return 0;
}
