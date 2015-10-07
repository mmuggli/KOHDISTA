#ifndef ALPHABET_H
#define ALPHABET_H

#include <cstdio>
#include <fstream>

#include "misc/definitions.h"
#include <map>


namespace CSA
{


class Alphabet
{
  public:
    explicit Alphabet(const std::map<usint, usint>& counts);
    explicit Alphabet(std::ifstream& file);
    explicit Alphabet(FILE* file);
    ~Alphabet() {}

    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;
    usint reportSize() const;
    inline bool isOk() const { return this->ok; }

    inline bool hasChar(usint c) const { return this->index_ranges.count(c); }
    inline usint countOf(usint c) const { return this->index_ranges.count(c) ?  length(this->index_ranges.at(c)) : 0; }
    // cumulative() is array C in CSA.
    inline usint cumulative(usint c) const { return this->index_ranges.at(c).first; }
    inline pair_type getRange(usint c) const { return this->index_ranges.at(c); }

    inline usint getDataSize() const { return this->size; }
    inline usint getAlphabetSize() const { return this->chars; }


    inline usint getFirstChar() const { return this->text_chars.at(0); }
    inline usint getTextChar(usint i) const { return this->text_chars.at(i); }
    inline std::map<usint, pair_type>::const_iterator begin() const { return index_ranges.begin();}
    inline std::map<usint, pair_type>::const_iterator end() const { return index_ranges.end(); }
    inline usint charAt(usint i) const
    {
        const usint* curr = &(this->text_chars.at(this->index_pointers.at(i / this->index_rate)));
      while(i > this->index_ranges.at(*curr).second) { curr++; }
      return *curr;
    }

  private:
    // index_ranges[c] is range for suffixes starting with 'c'.
    std::map<usint, pair_type> index_ranges;
    usint size;

    std::map<usint, usint> text_chars;  // which characters are present in the text
    usint chars;

    std::map<usint, usint> index_pointers; // which of the above is at i * index_rate
    usint index_rate;

    bool ok;

    void initialize(const std::map<usint, usint>& counts);

    // These are not allowed.
    Alphabet();
    Alphabet(const Alphabet&);
    Alphabet& operator = (const Alphabet&);
};


} // namespace CSA


#endif // ALPHABET_H
