#include <cstdlib>
#include <fstream>
#include <iostream>

#include "../misc/definitions.h"


using namespace CSA;


std::string gapChars = " =-_";

void
concatenateStrain(std::string &text, std::ifstream &file, usint strain)
{
  file.clear();
  file.seekg(0, std::ios::beg);

  std::string line;
  while(!file.eof())
  {
    getline(file, line);
    if(line.empty() || line[0] == '>') { continue; }
    char t = toupper(line[strain * 2]);
    if(gapChars.find(t) == std::string::npos) { text.push_back(t); }
  }
}


int
main(int argc, char** argv)
{
  std::cout << "Extracting concatenated text from alignment" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: extract_text alignment_file output_file" << std::endl;
    return 1;
  }

  std::cout << "Alignment file: " << argv[1] << std::endl;
  std::ifstream alignment_file(argv[1], std::ios_base::binary);
  if(!alignment_file)
  {
    std::cerr << "Error opening alignment file!" << std::endl;
    return 2;
  }

  std::cout << "Output file: " << argv[2] << std::endl;
  std::ofstream output_file(argv[2], std::ios_base::binary);
  if(!output_file)
  {
    std::cerr << "Error creating output file!" << std::endl;
    return 3;
  }
  std::cout << std::endl;

  std::string line;
  std::getline(alignment_file, line); // Skipping header line.
  std::getline(alignment_file, line);
  usint strains = line.size() / 2 - 1;
  std::cout << "Number of strains: " << strains << std::endl;

  line.clear();
  concatenateStrain(line, alignment_file, 0);
  usint n = line.size();
  std::cout << "Reference sequence length: " << n << std::endl;

  for(usint i = 1; i <= strains; i++)
  {
    concatenateStrain(line, alignment_file, i);
  }
  usint N = line.size() + 1;
  std::cout << "Total length: " << N << std::endl;

  output_file.write(line.c_str(), N - 1);
  char end = '\0';
  output_file.write(&end, 1);

  output_file.close();
  alignment_file.close();
  return 0;
}
