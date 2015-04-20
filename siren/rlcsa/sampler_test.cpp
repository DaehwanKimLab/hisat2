#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "rlcsa.h"


using namespace CSA;


const parameter_type INITIAL_ADJUSTMENT = parameter_type("INITIAL_ADJUSTMENT", 0);


void
printUsage()
{
  std::cout << "Usage: sampler_test [parameters] base_name weights" << std::endl;
  std::cout << "  -b   Use backward edge weights (locate for LF, display for Psi)." << std::endl;
  std::cout << "  -g#  Select # (float) part of samples greedily, and the rest at regular intervals." << std::endl;
  std::cout << "  -m#  Multiply weights by # (float)." << std::endl;
  std::cout << "  -t#  Use # threads." << std::endl;
  std::cout << "  -w   Write the sampled positions into sample file." << std::endl;
}



int main(int argc, char** argv)
{
  std::cout << "Sampler test" << std::endl;
  std::cout << std::endl;

  // Handle command line parameters.
  double multiplier = 1.0, greediness = 0.0;
  bool greedy = false, forward = true, write = false;
  char* base_name = 0;
  char* weights_name = 0;
  usint threads = 1;
  for(int i = 1; i < argc; i++)
  {
    if(argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
        case 'b':
          forward = false; break;
        case 'g':
          greediness = atof(argv[i] + 2); greedy = true; break;
        case 'm':
          multiplier = atof(argv[i] + 2); break;
        case 't':
          threads = atoi(argv[i] + 2); break;
        case 'w':
          write = true; break;
        default:
          std::cout << "Invalid option: " << argv[i] << std::endl << std::endl;
          printUsage();
          return 1;
      }
    }
    else if(base_name == 0)
    {
      base_name = argv[i];
    }
    else if(weights_name == 0)
    {
      weights_name = argv[i];
    }
  }
  if(weights_name == 0) { printUsage(); return 2; }

  std::cout << "Options:";
  greediness = std::min(1.0, std::max(0.0, greediness));
  threads = std::max(threads, (usint)1); threads = std::min(threads, (usint)64);
  if(!forward) { std::cout << " backward"; }
  if(greedy) { std::cout << " greediness=" << greediness; }
  if(multiplier < 0) { multiplier = 1.0; }
  if(multiplier != 1.0) { std::cout << " multiplier=" << multiplier; }
  std::cout << " threads=" << threads;
  if(write) { std::cout << " write"; }
  std::cout << std::endl;
  std::cout << "Base name: " << base_name << std::endl;
  std::cout << "Weights: " << weights_name << std::endl;
  std::cout << std::endl;

  // Determine index parameters.
  std::string parameters_name = std::string(base_name) + PARAMETERS_EXTENSION;
  Parameters parameters;
  parameters.set(RLCSA_BLOCK_SIZE);
  parameters.set(SAMPLE_RATE);
  parameters.read(parameters_name);
  parameters.set(SUPPORT_LOCATE.first, 1);
  parameters.set(SUPPORT_DISPLAY.first, 1);
  parameters.set(WEIGHTED_SAMPLES.first, 1);
  // parameters.print();


  std::ifstream weights_file(weights_name, std::ios_base::binary);
  if(!weights_file)
  {
    std::cerr << "Error opening weights file!" << std::endl;
    return 3;
  }

  usint n = fileSize(weights_file) / sizeof(weight_type);
  weight_type* weights = new weight_type[n];
  weights_file.read((char*)weights, n * sizeof(weight_type));
  weights_file.close();
  for(usint i = 0; i < n; i++)
  {
    if(weights[i] >= 0)
    {
      weights[i] = ceil(weights[i] * multiplier);
    }
    else
    {
      std::cerr << "Negative value found: weights[" << i << "] == " << weights[i] << std::endl;
      delete[] weights;
      return 4;
    }
  }
  std::cout << "Input size:   " << n << std::endl;
  std::cout << "Sample rate:  " << parameters.get(SAMPLE_RATE) << std::endl;
  std::cout << std::endl;


  // Do the actual sampling.
  double start = readTimer();
  Sampler* sampler = 0;
  if(!greedy)
  {
    WeightedSampler* wsampler = new WeightedSampler(weights, n, forward);
    if(!(wsampler->buildSamples(parameters.get(SAMPLE_RATE), parameters.get(INITIAL_ADJUSTMENT), threads)))
    {
      delete wsampler; wsampler = 0;
      return 5;
    }
    sampler = wsampler; wsampler = 0;
  }
  else
  {
    SemiGreedySampler* gsampler = new SemiGreedySampler(weights, n);
    if(!(gsampler->buildSamples(parameters.get(SAMPLE_RATE), greediness)))
    {
      delete gsampler; gsampler = 0;
      return 5;
    }
    sampler = gsampler; gsampler = 0;
  }
  double sampled = readTimer();
  usint memory = memoryUsage();


  // Write the samples or build RLCSA.
  double stop;
  if(write)
  {
    stop = readTimer();
    sampler->writeTo(base_name);
  }
  else
  {
    std::ifstream input_file(base_name, std::ios_base::binary);
    if(!input_file)
    {
      std::cerr << "Error opening input file " << base_name << std::endl;
      return 5;
    }
    uchar* data = new uchar[n];
    input_file.read((char*)data, n);
    input_file.close();
    RLCSA rlcsa(data, n, parameters.get(RLCSA_BLOCK_SIZE), parameters.get(SAMPLE_RATE), threads, sampler, true);
    stop = readTimer();
    if(rlcsa.isOk())
    {
      rlcsa.printInfo();
      rlcsa.reportSize(true);
      rlcsa.writeTo(base_name);
    }
  }
  usint total_memory = memoryUsage();


  double sampling_time = sampled - start;
  double time = stop - start;
  double megabytes = n / (double)MEGABYTE;
  std::cout << "Total: " << time << " seconds (" << (megabytes / time) << " MB/s), " << total_memory << " kB" << std::endl;
  std::cout << "Sampling: " << sampling_time << " seconds (" << (megabytes / sampling_time) << " MB/s), " << memory << " kB" << std::endl;
  std::cout << std::endl;

  delete sampler; sampler = 0;
  return 0;
}
