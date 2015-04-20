#ifndef ALPHABET_H
#define ALPHABET_H

#include <cstdio>
#include <fstream>

#include "misc/definitions.h"


namespace CSA
{


class Alphabet
{
  public:
    explicit Alphabet(const usint* counts);
    explicit Alphabet(std::ifstream& file);
    explicit Alphabet(FILE* file);
    ~Alphabet() {}

    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;
    usint reportSize() const;
    inline bool isOk() const { return this->ok; }

    inline bool hasChar(usint c) const { return !isEmpty(this->index_ranges[c]); }
    inline usint countOf(usint c) const { return length(this->index_ranges[c]); }
    // cumulative() is array C in CSA.
    inline usint cumulative(usint c) const { return this->index_ranges[c].first; }
    inline pair_type getRange(usint c) const { return this->index_ranges[c]; }

    inline usint getDataSize() const { return this->size; }
    inline usint getAlphabetSize() const { return this->chars; }

    inline usint getFirstChar() const { return this->text_chars[0]; }
    inline usint getTextChar(usint i) const { return this->text_chars[i]; }

    inline usint charAt(usint i) const
    {
      const usint* curr = &(this->text_chars[this->index_pointers[i / this->index_rate]]);
      while(i > this->index_ranges[*curr].second) { curr++; }
      return *curr;
    }

  private:
    // index_ranges[c] is range for suffixes starting with 'c'.
    pair_type index_ranges[CHARS];
    usint size;

    usint text_chars[CHARS];  // which characters are present in the text
    usint chars;

    usint index_pointers[CHARS]; // which of the above is at i * index_rate
    usint index_rate;

    bool ok;

    void initialize(const usint* counts);

    // These are not allowed.
    Alphabet();
    Alphabet(const Alphabet&);
    Alphabet& operator = (const Alphabet&);
};


} // namespace CSA


#endif // ALPHABET_H
