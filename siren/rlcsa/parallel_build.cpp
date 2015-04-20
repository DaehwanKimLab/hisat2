#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif

#include "rlcsa_builder.h"
#include "misc/utils.h"


using namespace CSA;


double getRLCSA(RLCSABuilder& builder, const std::string& base_name);
double indexParts(std::vector<std::string>& filenames, usint threads, Parameters& parameters, double& build_time);

const int MAX_THREADS = 64;


int
main(int argc, char** argv)
{
  std::cout << "Parallel RLCSA builder" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: parallel_build [-f|-n] listname output [threads]" << std::endl;
    std::cout << "  -f   use fast algorithm with larger memory usage" << std::endl;
    std::cout << "  -n   do not merge the indexes" << std::endl;
    return 1;
  }

  int list_parameter = 1, output_parameter = 2, threads_parameter = 3;
  bool do_merge = true;
  bool fast_algorithm = false;
  if(argv[1][0] == '-')
  {
    if(argv[1][1] == 'f')
    {
      list_parameter++; output_parameter++; threads_parameter++;
      fast_algorithm = true;
      std::cout << "Using fast algorithm." << std::endl;
    }
    else if(argv[1][1] == 'n')
    {
      list_parameter++; output_parameter++; threads_parameter++;
      do_merge = false;
      std::cout << "Option '-n' specified. Partial indexes will not be merged." << std::endl;
    }
  }

  std::ifstream filelist(argv[list_parameter], std::ios_base::binary);
  if(!filelist)
  {
    std::cerr << "Error opening file list!" << std::endl;
    return 2;
  }
  std::vector<std::string> files;
  readRows(filelist, files, true);
  filelist.close();
  std::cout << "Input files: " << files.size() << std::endl;

  std::string base_name = argv[output_parameter];
  std::cout << "Output: " << base_name << std::endl;

  usint threads = 1;
  if(argc > threads_parameter)
  {
    threads = std::min(MAX_THREADS, std::max(atoi(argv[threads_parameter]), 1));
  }
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
  double megabytes = 0.0, build_time = 0.0;
  RLCSABuilder builder(parameters.get(RLCSA_BLOCK_SIZE), parameters.get(SAMPLE_RATE), 0, threads);

  if(fast_algorithm)
  {
    for(std::vector<std::string>::iterator iter = files.begin(); iter != files.end(); ++iter)
    {
      double mark = readTimer();
      std::cout << "Increment: " << *iter; std::cout.flush();
      builder.insertCollection(*iter);
      std::cout << " (" << (readTimer() - mark) << " seconds)" << std::endl;
    }
    std::cout << std::endl;
    megabytes = getRLCSA(builder, base_name);
  }
  else
  {
    megabytes = indexParts(files, threads, parameters, build_time);
    if(do_merge)
    {
      std::cout << "Phase 2: Merging the indexes" << std::endl;
      for(std::vector<std::string>::iterator iter = files.begin(); iter != files.end(); ++iter)
      {
        double mark = readTimer();
        std::cout << "Increment: " << *iter; std::cout.flush();
        builder.insertFromFile(*iter);
        std::cout << " (" << (readTimer() - mark) << " seconds)" << std::endl;
      }
      std::cout << std::endl;
      megabytes = getRLCSA(builder, base_name);
    }
  }

  double stop = readTimer();
  std::cout << megabytes << " megabytes indexed in " << (stop - start) << " seconds (" << (megabytes / (stop - start)) << " MB/s)." << std::endl;
  if(do_merge)
  {
    std::cout << "Build time:    " << build_time + builder.getBuildTime() << " seconds" << std::endl;
    std::cout << "Search time:   " << builder.getSearchTime() << " seconds" << std::endl;
    std::cout << "Sort time:     " << builder.getSortTime() << " seconds" << std::endl;
    std::cout << "Merge time:    " << builder.getMergeTime() << " seconds" << std::endl;
  }
  std::cout << "Memory usage:  " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return 0;
}


double
getRLCSA(RLCSABuilder& builder, const std::string& base_name)
{
  double megabytes = 0.0;

  RLCSA* index = builder.getRLCSA();
  if(index != 0 && index->isOk())
  {
    index->printInfo();
    index->reportSize(true);
    index->writeTo(base_name);
    megabytes = index->getSize() / (double)MEGABYTE;
  }

  delete index;
  return megabytes;
}


double
indexParts(std::vector<std::string>& filenames, usint threads, Parameters& parameters, double& build_time)
{
  double start = readTimer();
  std::cout << "Phase 1: Building indexes for input files" << std::endl;
  usint block_size = parameters.get(RLCSA_BLOCK_SIZE);
  usint sample_rate = parameters.get(SAMPLE_RATE);
  usint total_size = 0;

  for(usint i = 0; i < filenames.size(); i++)
  {
    double init = readTimer();
    uchar* data = 0; usint size = 0;
    std::cout << "Input: " << filenames[i]; std::cout.flush();
    std::ifstream input(filenames[i].c_str(), std::ios_base::binary);
    if(!input)
    {
      std::cerr << "Error opening input file " << filenames[i] << std::endl;
    }
    else
    {
      size = fileSize(input);
      data = new uchar[size];
      input.read((char*)data, size);
      input.close();
    }

    if(size > 0)
    {
      double mark = readTimer();
      RLCSA* index = new RLCSA(data, size, block_size, sample_rate, threads, true);
      double done = readTimer();
      std::cout << " (" << (done - init) << " seconds)" << std::endl;
      build_time += done - mark;
      if(index != 0 && index->isOk()) { index->writeTo(filenames[i]); total_size += size; }
      delete index;
    }
    else
    {
      std::cerr << "Warning: Empty input file " << filenames[i] << std::endl;
    }
  }

  double total_time = readTimer() - start;
  double megabytes = total_size / (double)MEGABYTE;
  std::cout << "Indexed " << megabytes << " megabytes in " << total_time << " seconds (" << (megabytes / total_time) << " MB/s)." << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return megabytes;
}
