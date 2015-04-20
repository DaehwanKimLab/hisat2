#include <fstream>
#include <iostream>
#include <vector>

#include "../misc/definitions.h"
#include "../misc/utils.h"


using namespace CSA;


int
main(int argc, char** argv)
{
  std::cout << "Converting lines into a patterns file in Pizza & Chili format" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: convert_patterns input output" << std::endl;
    return 1;
  }

  std::cout << "Input file: " << argv[1] << std::endl;
  std::ifstream input_file(argv[1], std::ios_base::binary);
  if(!input_file)
  {
    std::cerr << "Error opening input file!" << std::endl;
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

  std::vector<std::string> rows;
  readRows(input_file, rows, true);

  usint patterns = rows.size();
  usint pattern_length = (*(rows.begin())).length();
  std::cout << "Number of patterns: " << patterns << std::endl;
  std::cout << "Pattern length:     " << pattern_length << std::endl;
  std::cout << std::endl;

  output_file << "# number=" << patterns << " length=" << pattern_length << " file=foo forbidden=" << std::endl;
  for(std::vector<std::string>::iterator iter = rows.begin(); iter != rows.end(); iter++)
  {
    output_file.write((*iter).c_str(), pattern_length);
  }

  output_file.close();
  input_file.close();
  return 0;
}
