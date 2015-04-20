#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <misc/definitions.h>
#include <misc/utils.h>


using namespace CSA;


std::string gapChars = " =-_";

int
main(int argc, char** argv)
{
  std::cout << "Alignment transposer" << std::endl;
  std::cout << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: transpose alignment_file output_file [endmarker] [multiple files]" << std::endl;
    return 1;
  }

  std::cout << "Alignment file: " << argv[1] << std::endl;
  std::ifstream alignment_file(argv[1], std::ios_base::binary);
  if(!alignment_file)
  {
    std::cerr << "Error opening alignment file!" << std::endl;
    return 2;
  }

  std::string output_name = argv[2];
  std::cout << "Output file: " << output_name << std::endl;
  std::ofstream* output_file = 0;

  usint endmarker = '\n';
  if(argc >= 4) { endmarker = atoi(argv[3]); }
  std::cout << "End marker: " << endmarker << std::endl;
  bool multiple_files = (argc >= 5);
  if(multiple_files) { std::cout << "Using separate output files for the sequences." << std::endl; }

  // Determine the number of sequences.
  std::string line;
  std::getline(alignment_file, line);
  usint sequences = line.size();
  usint filesize = fileSize(alignment_file);
  std::cout << "Sequences: " << sequences << std::endl;
  std::cout << std::endl;

  usint counts[CHARS];
  for(usint c = 0; c < CHARS; c++) { counts[c] = 0; }

  if(!multiple_files) { output_file = new std::ofstream(output_name.c_str(), std::ios_base::binary); }
  for(usint i = 0; i < sequences; i++)
  {
    if(multiple_files)
    {
      std::stringstream filename;
      filename << output_name << "." << i;
      output_file = new std::ofstream(filename.str().c_str(), std::ios_base::binary);
    }

    alignment_file.clear();
    alignment_file.seekg(0, std::ios_base::beg);
    usint len = 0;
    if(!multiple_files && i == 0 && sequences > filesize / 2)
    {
      getline(alignment_file, line);
      for(; i < sequences; i++)
      {
        if(gapChars.find(line[i]) == std::string::npos && line[i] != 0)
        {
          output_file->put(line[i]);
          counts[(int)line[i]]++;
          output_file->put(endmarker);
        }
      }
      break;
    }
    while(!alignment_file.eof())
    {
      getline(alignment_file, line);
      if(line.empty()) { continue; }
      if(gapChars.find(line[i]) == std::string::npos && line[i] != 0)
      {
        output_file->put(line[i]); len++;
        counts[(int)line[i]]++;
      }
    }
    output_file->put(endmarker);
    if(multiple_files) { delete output_file; }

    std::cout << "Sequence " << i << ", length " << len << std::endl;
  }
  if(!multiple_files) { delete output_file; }
  std::cout << std::endl;

  std::cout << "Character counts: " << std::endl;
  for(usint c = 0; c < CHARS; c++)
  {
    if(counts[c] == 0) { continue; }
    std::cout << "  " << ((char)c) << " (" << c << "): " << counts[c] << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
