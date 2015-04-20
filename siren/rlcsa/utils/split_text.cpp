#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../misc/definitions.h"
#include "../misc/utils.h"


using namespace CSA;


int
main(int argc, char** argv)
{
  std::cout << "Text splitter" << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: split_text file parts" << std::endl;
    return 1;
  }

  std::string base_name = argv[1];
  std::cout << "Input file: " << base_name << std::endl;
  std::ifstream input_file(base_name.c_str(), std::ios_base::binary);
  if(!input_file)
  {
    std::cerr << "Error opening input file!" << std::endl;
    return 2;
  }

  usint parts = atoi(argv[2]);
  usint size = fileSize(input_file);
  usint part_size = (size + parts - 1) / parts;
  std::cout << "Parts: " << parts << std::endl;
  std::cout << "Part size: " << part_size << " bytes" << std::endl;
  if(parts == 0 || part_size == 0)
  {
    std::cerr << "Invalid number of parts!" << std::endl;
    return 3;
  }
  std::cout << std::endl;

  char* buffer = new char[part_size + 1];
  usint len = 0;
  for(usint n = parts; n > 0; n /= 10) { len++; }

  for(usint i = 1; i <= parts; i++)
  {
    usint llen = len;
    for(usint num = i; num > 0; num /= 10) { llen--; }

    std::ostringstream out;
    out << base_name << "." << std::string(llen, '0') << i;
    std::ofstream part_file(out.str().c_str(), std::ios_base::binary);
    input_file.read(buffer, part_size);
    usint n = input_file.gcount();
    buffer[n] = 0;
    part_file.write(buffer, n + 1);
    part_file.close();
  }

  input_file.close();
  return 0;
}
