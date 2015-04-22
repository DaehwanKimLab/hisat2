#include <cstdlib>
#include <fstream>
#include <iostream>

#include <misc/utils.h>


using namespace CSA;


std::string gapChars = " =-_";
std::string bases = "ACGNT";

char*
readAlignment(std::ifstream& input, usint& lines, usint sequences, bool quality_included)
{
  input.clear();
  input.seekg(0, std::ios_base::beg);

  usint chars[CHARS];
  for(usint c = 0; c < CHARS; c++) { chars[c] = 0; }
  char* data = new char[fileSize(input)];

  std::string line;
  if(quality_included) { getline(input, line); }
  usint multi = (quality_included ? 2 : 1);
  while(!input.eof())
  {
    getline(input, line);
    if(line.empty()) { continue; }
    std::string output(sequences + 1, '\n');

    for(usint i = 0; i < sequences; i++)
    {
      char t = toupper(line[i * multi]);
      if(gapChars.find(t) != std::string::npos) { t = '-'; }
      chars[(int)t]++;
      data[lines * (sequences + 1) + i] = t;
    }
    data[lines * (sequences + 1) + sequences] = '\n';
    lines++;
  }

  std::cout << "Character counts:" << std::endl;
  for(usint c = 0; c < CHARS; c++)
  {
    if(chars[c] > 0) { std::cout << "  " << (char)c << " (" << c << "): " << chars[c] << std::endl; }
  }
  std::cout << std::endl;

  std::cout << "Lines (initial): " << lines << std::endl;
  return data;
}

// Replace IUPAC codes denoting multiple possible bases with 'N'.
// Replace garbage with gaps.
usint
filterRow(char* data, usint sequences)
{
  usint replaced = 0;
  for(usint i = 0; i < sequences; i++)
  {
    switch(data[i])
    {
      case 'A':
      case 'C':
      case 'G':
      case 'N':
      case 'T':
      case '-':
        break;  // These are ok.

      case 'R':
      case 'Y':
      case 'K':
      case 'M':
      case 'S':
      case 'W':
      case 'B':
      case 'D':
      case 'H':
      case 'V':
        data[i] = 'N'; replaced++; break;

      default:
        data[i] = '-'; replaced++; break;
    }
  }

  return replaced;
}


struct AlignmentPosition
{
  char  current, previous;
  usint prev_line;
  usint id;

  void findPrevious(char* data, usint sequences)
  {
    this->current = this->previous;
    while(this->prev_line > 0)
    {
      this->prev_line--;
      if(data[this->prev_line * (sequences + 1) + this->id] != '-')
      {
        this->previous = data[this->prev_line * (sequences + 1) + this->id];
        this->prev_line++; break;
      }
    }
    if(this->prev_line == 0) { this->previous = 'Z'; }
    else                     { this->prev_line--; }
  }

  void print()
  {
    std::cout << "(" << this->current << ", " << this->previous << ", " << this->prev_line << ", " << this->id << ")" << std::endl;
  }
};

struct PositionComparator
{
  bool operator() (const AlignmentPosition* a, const AlignmentPosition* b) const
  {
    if(a->current < b->current)   { return true; }
    if(a->current > b->current)   { return false; }

    if(a->previous < b->previous) { return true; }
    if(a->previous > b->previous) { return false; }

    return (a->prev_line < b->prev_line);
  }
} pos_comparator;

usint
handleSequences(AlignmentPosition** active, pair_type range, char* data, usint sequences)
{
  if(active[range.first]->previous == 'Z') { return 0; }

  usint moved = 0;
  usint max_prev = active[range.second]->prev_line;
  for(usint i = range.first; i < range.second; i++)
  {
    if(active[i]->prev_line != max_prev)
    {
      std::swap(data[active[i]->prev_line * (sequences + 1) + active[i]->id],
                data[max_prev * (sequences + 1) + active[i]->id]);
      active[i]->prev_line = max_prev;
      moved++;
    }
  }

  return moved;
}

// Each gap must be preceded a character that does not occur at the last line of the gap.
// Additionally, two gaps ending at the same line but beginning at different lines must be
// preceded by different characters.
usint
handleGaps(char* data, std::ofstream& output, usint sequences, usint lines)
{
  AlignmentPosition seqs[sequences];
  for(usint i = 0; i < sequences; i++)
  {
    seqs[i].previous = '\0';
    seqs[i].prev_line = lines;
    seqs[i].id = i;
  }
  bool* ok = new bool[lines + 1];

  usint moved = 0, deleted = 0;
  for(usint l = 0; l <= lines; l++)
  {
    usint i = lines - l;  // Actual line.

    // Find the previous non-gap character for each sequence.
    AlignmentPosition* active[sequences];
    usint active_seqs = 0;
    for(usint j = 0; j < sequences; j++)
    {
      if(seqs[j].prev_line != i) { continue; }  // Inside a gap.
      seqs[j].findPrevious(data, sequences);
      active[active_seqs] = &(seqs[j]); active_seqs++;
    }

    // Sort active sequences by current, previous, prev_line.
    // For sequences with same (current, previous), set prev_line to
    // maximum among those sequences, and update the alignment.
    if(active_seqs > 0)
    {
      sequentialSort(active, active + active_seqs, pos_comparator);
      usint prev = 0;
      for(usint j = 1; j < active_seqs; j++)
      {
        if(active[j]->current != active[prev]->current ||
           active[j]->previous != active[prev]->previous)
        {
          moved += handleSequences(active, pair_type(prev, j - 1), data, sequences);
          prev = j;
        }
      }
      moved += handleSequences(active, pair_type(prev, active_seqs - 1), data, sequences);
      ok[i] = true;
    }
    else { deleted++; ok[i] = false; }
  }

  usint truelines = 0;
  for(usint i = 0; i < lines; i++)
  {
    if(ok[i]) { output.write(data + i * (sequences + 1), sequences + 1); truelines++; }
  }
  output.close();

  std::cout << "Lines (final): " << truelines << std::endl;
  std::cout << "Moved " << moved << " characters." << std::endl;
  std::cout << "Deleted " << deleted << " lines." << std::endl;
  std::cout << std::endl;

  delete[] ok;
  return truelines;
}


int
main(int argc, char** argv)
{
  std::cout << "Alignment cleaner" << std::endl;
  std::cout << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: clean_alignment alignment_file output_file" << std::endl;
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

  double start = readTimer();

  // Determine the format and the number of sequences.
  std::string line;
  bool quality_included = false;
  std::getline(alignment_file, line);
  usint sequences = line.size();
  if(line[0] == '>')
  {
    quality_included = true; 
    std::getline(alignment_file, line);
    sequences = line.size() / 2;
  }
  std::cout << "Sequences: " << sequences << std::endl;
  std::cout << std::endl;

  std::vector<std::string> alignment;
  usint lines = 0;
  char* data = readAlignment(alignment_file, lines, sequences, quality_included);
  alignment_file.close();

  usint replaced = 0;
  for(usint i = 0; i < lines; i++)
  {
    replaced += filterRow(data + i * (sequences + 1), sequences);
  }
  std::cout << "Replaced " << replaced << " characters." << std::endl;

  handleGaps(data, output_file, sequences, lines);
  delete[] data;
  output_file.close();

  double time = readTimer() - start;
  std::cout << "Used " << time << " seconds." << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;
  
  return 0;
}
