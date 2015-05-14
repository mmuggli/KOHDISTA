#include <fstream>
#include <iostream>

#include "alphabet.h"
#include <assert.h>  


namespace CSA
{

//--------------------------------------------------------------------------

    Alphabet::Alphabet(const std::map<usint, usint>& counts) :
  ok(false)
{
  this->initialize(counts);
}

Alphabet::Alphabet(std::ifstream& file) :
  ok(false)
{
//  usint counts[CHARS];
  std::map<usint,usint> counts;
//  assert (!"FIXME: next comment line doesn't use new counts properly");
  while (1) {
      usint key = 0;
      usint val = 0;
      file.read((char*)&key, sizeof(key));
      file.read((char*)&val, sizeof(val));
      if (key==0 && val == 0) break;
      counts[key] = val;
  }
  //file.read((char*)counts, CHARS * sizeof(usint));
  this->initialize(counts);
}

Alphabet::Alphabet(FILE* file) :
  ok(false)
{
  if(file == 0) { return; }
  //usint counts[CHARS];
  std::map<usint,usint> counts;
  //assert (!"FIXME: next comment line doesn't use new counts properly");
//  if(!std::fread(counts, CHARS * sizeof(usint), 1, file)) { return; }
  while (1) {
      usint key = 0;
      usint val = 0;
      std::fread((char*)&key, sizeof(key), 1, file);
      std::fread((char*)&val, sizeof(val), 1, file);
      if (key==0 && val == 0) break;
      counts[key] = val;
  }

  this->initialize(counts);
}

void
Alphabet::initialize(const std::map<usint, usint>& counts)
{
    if(counts.empty()) { return; }

  this->size = 0; this->chars = 0;
  //for(usint c = 0; c < CHARS; c++)
  int symcnt = 0;
  for(std::map<usint, usint>::const_iterator mapiter = counts.begin(); mapiter != counts.end(); ++mapiter)
  {
      symcnt++;
      if (symcnt % 1000 == 0) std::cout << "alpha::init-1'd " << symcnt << " symbols." << std::endl;

    this->index_ranges[mapiter->first] = ((mapiter->second > 0 || this->size > 0) ?
                             pair_type(this->size, this->size + mapiter->second - 1) :
                             EMPTY_PAIR);
    if(mapiter->second > 0)
    {
      this->text_chars[this->chars] = mapiter->first;
      this->chars++;
    }
    size += mapiter->second;
  }

  this->index_rate = std::max((this->size + CHARS - 1) / CHARS, (usint)1);
  std::cout << "Alphabet::index_rate: " << this->index_rate << std::endl;
  usint current = 0;
  symcnt = 0;
  for(usint c = 0, i = 0; c < this->chars; c++)
  {
      symcnt++;
      if (symcnt % 1000 == 0) std::cout << "alpha::init-2'd " << symcnt << " symbols." << std::endl;
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
//    for(usint c = 0; c < 256/*FIXME:CHARS*/; c++) {
    for( std::map<usint, pair_type>::const_iterator itr = this->begin(); itr != this->end(); ++itr){
        usint c = itr->first;
        usint temp = this->countOf(c);
        if (temp != 0) {
            file.write((char*)&c, sizeof(c));//FIXME: possibly not portable wrt endianness
            file.write((char*)&temp, sizeof(temp));//FIXME: possibly not portable wrt endianness
        }
    }
    const usint zero = 0;
    //write a centinel map[0]=0 entry
    file.write((char*)&zero, sizeof(zero));//FIXME: possibly not portable wrt endianness
    file.write((char*)&zero, sizeof(zero));//FIXME: possibly not portable wrt endianness
  
  // for(usint c = 0; c < CHARS; c++)
  // {
  //   usint temp = this->countOf(c);
    
  // }
}

void
Alphabet::writeTo(FILE* file) const
{
    if(file == 0) { return; }
    for(usint c = 0; c < CHARS; c++) {
        usint temp = this->countOf(c);
        if (temp != 0) {
            std::fwrite((char*)&c, sizeof(c), 1, file);//FIXME: possibly not portable wrt endianness
            std::fwrite((char*)&temp, sizeof(temp),1, file);//FIXME: possibly not portable wrt endianness
        }
    }
    const usint zero = 0;
        //write a centinel map[0]=0 entry
    std::fwrite((char*)&zero, sizeof(zero), 1, file);//FIXME: possibly not portable wrt endianness
    std::fwrite((char*)&zero, sizeof(zero),1, file);//FIXME: possibly not portable wrt endianness
  

        // for(usint c = 0; c < CHARS; c++)
 //  {
 //    usint temp = this->countOf(c);
 //    std::fwrite(&temp, sizeof(temp), 1, file);
 //  }
}

usint
Alphabet::reportSize() const
{
  return sizeof(*this);
}

//--------------------------------------------------------------------------

} // namespace CSA
