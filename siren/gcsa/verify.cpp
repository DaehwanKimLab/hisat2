#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include <misc/utils.h>

#include "gcsa.h"


using namespace CSA;


int main(int argc, char** argv)
{
  std::cout << "GCSA verifier" << std::endl;
  std::cout << std::endl;

  if(argc < 3)
  {
    std::cerr << "Usage: verify base_name input" << std::endl;
    return 1;
  }

  const GCSA gcsa(argv[1]);
  if(!gcsa.isOk()) { return 2; }
  gcsa.reportSize(true);

  std::ifstream input(argv[2], std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error: Cannot open input file " << argv[2] << std::endl;
    return 2;
  }

  DeltaVector::Iterator** array_iters = new DeltaVector::Iterator*[CHARS];
  for(usint i = 0; i < CHARS; i++)
  {
    array_iters[i] = gcsa.getIterator(i);
  }
  RLEVector::Iterator* edge_iter = gcsa.getEdgeIterator();

  std::cout << "Starting verification..." << std::endl;
  usint lines = 0, prev = edge_iter->select(0);
  while(!input.eof())
  {
    std::string line; getline(input, line);
    if(line.empty()) { continue; }
    lines++;
    bool ok = true;

    bool predecessors[CHARS];
    bool correct_predecessors[CHARS];
    for(usint i = 0; i < CHARS; i++)
    {
      predecessors[i] = false;
      correct_predecessors[i] = (array_iters[i] != 0 ? array_iters[i]->isSet(lines - 1) : false);
    }
    for(usint i = 0; i < line.size(); i++)
    {
      predecessors[(uchar)line[i]] = false;
    }
    for(usint i = 0; i < CHARS; i++)
    {
      if(predecessors[i] != correct_predecessors[i]) { ok = false; }
    }

    usint ones = 0;
    for(usint i = 0; i < line.size(); i++)
    {
      if(line[i] == '1') { ones++; }
    }
    usint next = edge_iter->selectNext();
    if(next + 1 - prev != ones) { ok = false; }

    if(!ok)
    {
      std::cout << "On line " << line << ", predecessors = ";
      for(usint i = 0; i < CHARS; i++)
      {
        if(predecessors[i]) { std::cout << (char)i; }
      }
      std::cout << " / ";
      for(usint i = 0; i < CHARS; i++)
      {
        if(correct_predecessors[i]) { std::cout << (char)i; }
      }
      std::cout << ", edges = " << (next - prev) << " / " << (ones - 1) << std::endl;
    }
    prev = next;
    if(lines % (10 * MILLION) == 0)
    {
      std::cout << "Processed " << lines << " lines..." << std::endl;
    }
  }
  std::cout << "Done! Total lines: " << lines << std::endl;

  for(usint i = 0; i < CHARS; i++)
  {
    delete array_iters[i];
  }
  delete[] array_iters; delete edge_iter;
  input.close();
  return 0;
}
