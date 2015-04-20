#include <cstdlib>
#include <fstream>
#include <iostream>

#include <misc/definitions.h>

const char letters[4] = { 'A', 'C', 'G', 'T' };


using namespace CSA;


int
main(int argc, char** argv)
{
  std::cout << "Data generator" << std::endl;
  std::cout << std::endl;
  if(argc < 6)
  {
    std::cout << "Usage: generator input output length sequences rate [seed]" << std::endl;
    std::cout << "    input          input file" << std::endl;
    std::cout << "    output         output file" << std::endl;
    std::cout << "    length         sequence length" << std::endl;
    std::cout << "    sequences      number of sequences" << std::endl;
    std::cout << "    rate           mutation rate" << std::endl;
    std::cout << "    seed           random seed (optional; default 0xDEADBEEF)" << std::endl;
    std::cout << std::endl;
    return 1;
  }

  std::cout << "Input: " << argv[1] << std::endl;
  std::ifstream input(argv[1], std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error opening input file!" << std::endl;
    return 2;
  }
  std::cout << "Output: " << argv[2] << std::endl;
  std::ofstream output(argv[2], std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error opening output file!" << std::endl;
    return 3;
  }
  usint bases = atol(argv[3]);
  std::cout << "Sequence length: " << bases << " bases" << std::endl;
  usint sequences = atoi(argv[4]);
  std::cout << "Number of sequences: " << sequences << std::endl;
  double rate = atof(argv[5]);
  std::cout << "Mutation rate: " << rate << std::endl;
  usint seed = 0xDEADBEEF;
  if(argc >= 7)
  {
    seed = atoi(argv[6]);
  }
  std::cout << "Random seed: " << seed << std::endl;
  std::cout << std::endl;

  srand(seed);
  usint mutations = 0;
  char* source = new char[bases];
  char* buffer = new char[sequences + 1]; buffer[sequences] = '\n';
  input.read(source, bases);
  usint size = input.gcount();

  bases = 0;
  for(usint i = 0; i < size; i++)
  {
    if(source[i] != letters[0] && source[i] != letters[1] &&
       source[i] != letters[2] && source[i] != letters[3]) { continue; }

    bases++;
    buffer[0] = source[i];
    for(usint j = 1; j < sequences; j++)
    {
      if(rand() / (RAND_MAX + 1.0) < rate)
      {
        usint tmp;
        switch(source[i])
        {
          case 'A':
            tmp = rand() % 3 + 1;
            break;
          case 'C':
            tmp = rand() % 3 + 2;
            break;
          case 'G':
            tmp = rand() % 3 + 3;
            break;
          case 'T':
            tmp = rand() % 3;
            break;
          default:
            tmp = rand() & 3;
            break;
        }
        buffer[j] = letters[tmp & 3];
        mutations++;
      }
      else
      {
        buffer[j] = source[i];
      }
    }
    output.write(buffer, sequences + 1);
  }

  std::cout << "Actual sequence length: " << bases << " bases" << std::endl;
  std::cout << "Number of mutations: " << mutations << std::endl;
  std::cout << std::endl;

  return 0;
}
