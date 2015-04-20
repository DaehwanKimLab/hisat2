#include <fstream>
#include <iostream>

#include "alphabet.h"



namespace CSA
{

//--------------------------------------------------------------------------

Alphabet::Alphabet(const usint* counts) :
  ok(false)
{
  this->initialize(counts);
}

Alphabet::Alphabet(std::ifstream& file) :
  ok(false)
{
  usint counts[CHARS];
  file.read((char*)counts, CHARS * sizeof(usint));
  this->initialize(counts);
}

Alphabet::Alphabet(FILE* file) :
  ok(false)
{
  if(file == 0) { return; }
  usint counts[CHARS];
  if(!std::fread(counts, CHARS * sizeof(usint), 1, file)) { return; }
  this->initialize(counts);
}

void
Alphabet::initialize(const usint* counts)
{
  if(counts == 0) { return; }

  this->size = 0; this->chars = 0;
  for(usint c = 0; c < CHARS; c++)
  {
    this->index_ranges[c] = ((counts[c] > 0 || this->size > 0) ?
                             pair_type(this->size, this->size + counts[c] - 1) :
                             EMPTY_PAIR);
    if(counts[c] > 0)
    {
      this->text_chars[this->chars] = c;
      this->chars++;
    }
    size += counts[c];
  }

  this->index_rate = std::max((this->size + CHARS - 1) / CHARS, (usint)1);
  usint current = 0;

  for(usint c = 0, i = 0; c < this->chars; c++)
  {
    pair_type range = this->index_ranges[this->text_chars[c]];
    while(current <= range.second)
    {
      this->index_pointers[i] = c;
      current += this->index_rate;
      i++;
    }
  }

  this->ok = true;
}

void
Alphabet::writeTo(std::ofstream& file) const
{
  for(usint c = 0; c < CHARS; c++)
  {
    usint temp = this->countOf(c);
    file.write((char*)&temp, sizeof(temp));
  }
}

void
Alphabet::writeTo(FILE* file) const
{
  if(file == 0) { return; }
  for(usint c = 0; c < CHARS; c++)
  {
    usint temp = this->countOf(c);
    std::fwrite(&temp, sizeof(temp), 1, file);
  }
}

usint
Alphabet::reportSize() const
{
  return sizeof(*this);
}

//--------------------------------------------------------------------------

} // namespace CSA
