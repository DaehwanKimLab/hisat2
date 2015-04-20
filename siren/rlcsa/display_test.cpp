#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "rlcsa.h"
#include "sampler.h"
#include "misc/utils.h"


using namespace CSA;


usint*
generatePositions(weight_type* weights, usint size, usint number);

int main(int argc, char** argv)
{
  std::cout << "RLCSA display test" << std::endl;
  if(argc < 5)
  {
    std::cout << "Usage: display_test basename weights sequences length [random_seed]" << std::endl;
    return 1;
  }
  std::cout << std::endl;

  std::cout << "Base name: " << argv[1] << std::endl;
  std::cout << "Weights: " << argv[2] << std::endl; 

  usint sequences = std::max(atoi(argv[3]), 1);
  std::cout << "Sequences: " << sequences << std::endl;

  usint seq_length = std::max(atoi(argv[4]), 1);
  std::cout << "Sequence length: " << seq_length << std::endl;

  usint seed = 0xDEADBEEF;
  if(argc > 5)
  {
    seed = atoi(argv[5]);
  }
  std::cout << "Random seed: " << seed << std::endl; 
  std::cout << std::endl;

  RLCSA rlcsa(argv[1]);
  if(!rlcsa.isOk() || !rlcsa.supportsDisplay())
  {
    std::cerr << "Error: Display is not supported!" << std::endl;
    return 2;
  }
  rlcsa.printInfo();
  rlcsa.reportSize(true);
  usint size = rlcsa.getSize();


  // Read weights and determine positions.
  weight_type* weights = new weight_type[size];
  std::ifstream weight_file(argv[2], std::ios_base::binary);
  if(!weight_file)
  {
    std::cerr << "Error: Cannot open weight file!" << std::endl;
    delete[] weights;
    return 3;
  }
  weight_file.read((char*)weights, size * sizeof(weight_type));
  weight_file.close();
  srand(seed);
  usint* positions = generatePositions(weights, size, sequences);
  delete[] weights;
  std::cout << "Generated " << sequences << " positions." << std::endl;
  std::cout << std::endl;


  // Do & measure the actual work.
  uchar* buffer = new uchar[seq_length];
  double start = readTimer();
  for(usint i = 0; i < sequences; i++)
  {
    rlcsa.display(0, pair_type(positions[i], positions[i] + seq_length - 1), buffer);
  }
  double time = readTimer() - start;
  std::cout << "Total time: " << time << " seconds." << std::endl;
  std::cout << std::endl;


  delete[] buffer;
  delete[] positions;
  return 0;
}


usint*
generatePositions(weight_type* weights, usint size, usint number)
{
  DeltaEncoder weight_encoder(2 * sizeof(usint));
  DeltaEncoder pos_encoder(2 * sizeof(usint));
  usint sum = 0;
  for(usint i = 0; i < size; i++)
  {
    if(weights[i] > 0)
    {
      sum += weights[i];
      weight_encoder.addBit(sum);
      pos_encoder.addBit(i);
    }
  }
  weight_encoder.flush();
  pos_encoder.flush();
  DeltaVector weight_vector(weight_encoder, sum);
  DeltaVector pos_vector(pos_encoder, size);

  usint* result = new usint[number];
  DeltaVector::Iterator weight_iter(weight_vector);
  DeltaVector::Iterator pos_iter(pos_vector);
  for(usint i = 0; i < number; i++)
  {
    result[i] = pos_iter.select(weight_iter.rank(rand() % sum));
  }

  return result;
}
