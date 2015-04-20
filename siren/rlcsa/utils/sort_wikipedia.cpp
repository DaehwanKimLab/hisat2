#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "../misc/utils.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif


using namespace CSA;


struct Document
{
  std::string key;
  char* start;
  char* stop;
};

struct DocumentComparator
{
  bool operator() (const Document& a, const Document& b)
  {
    return (a.key.compare(b.key) < 0);
  }
} document_comparator;


const std::string REVISION = "<revision>";
const std::string REVISION_END = "</revision>";
const std::string TIMESTAMP = "<timestamp>";
const std::string TIMESTAMP_END = "</timestamp>";


char*
find(char** from, const std::string& pattern, bool after)
{
  if(*from == 0) { return 0; }
  char* match = strstr(*from, pattern.c_str());
  if(match == 0) { *from = 0; return 0; }

  if(after) // Next line after the pattern.
  {
    char* tmp = strchr(match, '\n');
    if(tmp == 0) { *from = 0; return 0;}
    *from = match = tmp + 1;
  }
  else  // Lines until the one with the pattern.
  {
    *match = '\0';
    char* tmp = strrchr(*from, '\n');
    *match = pattern[0];
    if(tmp == 0) { *from = 0; return 0; }
    match = tmp + 1;
    tmp = strchr(tmp + 1, '\n');
    if(tmp == 0) { *from = 0; } else { *from = tmp + 1; }
  }

  return match;
}

std::string
between(char* from, char* to, const std::string& start, const std::string& stop)
{
  char tmp = *to;
  *to = '\0';

  from = strstr(from, start.c_str());
  if(from == 0) { *to = tmp; return ""; }
  from += start.length();
  char* limit = strstr(from, stop.c_str());
  if(limit == 0) { *to = tmp; return ""; }

  *to = tmp;
  return std::string(from, limit - from);
}


int
main(int argc, char** argv)
{
  double start = readTimer();

  std::cout << "Wikipedia sorter" << std::endl;
  std::cout << std::endl;

  if(argc < 3)
  {
    std::cout << "Usage: sort_wikipedia input part_size" << std::endl;
    return 1;
  }
  std::cout << "Input: " << argv[1] << std::endl;
  usint part_size = atoi(argv[2]);
  std::cout << "Part size: " << part_size << " megabytes" << std::endl;
  part_size *= MEGABYTE;
  if(part_size == 0)
  {
    std::cerr << "Error: Invalid part size!" << std::endl;
    return 2;
  }

  std::ifstream input(argv[1], std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error: Cannot open input file!" << std::endl;
    return 3;
  }
  usint size = fileSize(input);
  std::cout << "Size: " << (size / (double)MEGABYTE) << " megabytes" << std::endl;
  std::cout << std::endl;

  char* data = new char[size + 1];
  input.read(data, size); data[size] = '\0';
  input.close();

  std::vector<Document> documents;
  char* current = data;
  while(current != 0)
  {
    Document buffer;
    buffer.start = find(&current, REVISION, true);
    buffer.stop = find(&current, REVISION_END, false);
    if(buffer.start == 0 || buffer.stop == 0) { break; }
    buffer.key = between(buffer.start, buffer.stop, TIMESTAMP, TIMESTAMP_END);
    if(buffer.key.length() > 0) { documents.push_back(buffer); }
  }
  std::cout << "Sequences: " << documents.size() << std::endl;
  std::cout << std::endl;

  std::sort(documents.begin(), documents.end(), document_comparator);


  usint output_files = 0;
  std::ofstream* output = 0;
  usint output_size = 0;
  for(usint i = 0; i < documents.size(); i++)
  {
    usint doc_size = documents[i].stop - documents[i].start;
    if(doc_size + 1 + output_size > part_size && output_size > 0)
    {
      delete output; output = 0;
    }
    if(output == 0)
    {
      output_files++; output_size = 0;
      std::stringstream str(std::stringstream::out);
      str << "part." << output_files; str.flush();
      std::cout << "Writing " << str.str() << std::endl;
      output = new std::ofstream(str.str().c_str(), std::ios_base::binary);
    }
    output->write(documents[i].start, doc_size);
    output->put('\0');
    output_size += doc_size + 1;
  }
  delete output; output = 0;
  std::cout << std::endl;

  std::cout << "Files: " << output_files << std::endl;
  std::cout << "Time: " << (readTimer() - start) << " seconds" << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;
}
