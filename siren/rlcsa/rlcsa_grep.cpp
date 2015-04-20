#include <iostream>

#include "rlcsa.h"
#include "misc/utils.h"

using namespace CSA;


enum mode_type { COUNT, DISPLAY, CONTEXT };


void printUsage()
{
  std::cout << "Usage: rlcsa_grep [-c] pattern base_name" << std::endl;
  std::cout << "  -c    print the number of matching lines" << std::endl;
  std::cout << "  -NUM  display NUM characters of leading and trailing context instead of" << std::endl;
  std::cout << "        the entire line" << std::endl;
}


int main(int argc, char** argv)
{
  int base_arg = 2, pattern_arg = 1;
  mode_type mode = DISPLAY;
  usint context = 0;

  if(argc < base_arg + 1)
  {
    printUsage();
    return 1;
  }

  if(argv[1][0] == '-')
  {
    base_arg++; pattern_arg++;
    if(std::string("-c").compare(argv[1]) == 0)
    {
      mode = COUNT;
    }
    else
    {
      mode = CONTEXT;
      context = atoi(&(argv[1][1]));
    }
    if(argc < base_arg + 1)
    {
      printUsage();
      return 2;
    }
  }

  RLCSA rlcsa(argv[base_arg], false);
  if(!rlcsa.isOk())
  {
    return 3;
  }

  usint len = std::string(argv[pattern_arg]).length();
  pair_type result_range = rlcsa.count(argv[pattern_arg]);
  usint occurrences = length(result_range);
  if(occurrences == 0)
  {
    if(mode == COUNT)
    {
      std::cout << 0 << std::endl;
    }
    return 0;
  }

  usint last_row = 0;
  usint* results = rlcsa.locate(result_range);
  if(mode != CONTEXT)
  {
    rlcsa.getSequenceForPosition(results, occurrences);
  }
  std::sort(results, results + occurrences);
  if(mode != CONTEXT)
  {
    for(usint i = 1; i < occurrences; i++)
    {
      if(results[i] != results[last_row])
      {
        last_row++; results[last_row] = results[i];
      }
    }
  }

  if(mode == COUNT)
  {
    std::cout << (last_row + 1) << std::endl;
  }
  else if(mode == DISPLAY)
  {
    for(usint i = 0; i <= last_row; i++)
    {
      uchar* row = rlcsa.display(results[i]);
      std::cout.write((char*)row, length(rlcsa.getSequenceRange(results[i])));
      std::cout << std::endl;
      delete[] row;
    }
  }
  else if(mode == CONTEXT)
  {
    usint result_length = 0;
    for(usint i = 0; i < occurrences; i++)
    {
      uchar* text = rlcsa.display(results[i], len, context, result_length);
      std::cout.write((char*)text, result_length);
      std::cout << std::endl;
      delete[] text;
    }
  }

  delete[] results;
  return 0;
}
