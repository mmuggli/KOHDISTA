#include <cstdlib>

#include "sdslvector.h"
#include "../misc/utils.h"


namespace CSA
{
  BitVector::Iterator*  SDSLVector::newIterator()
  {
    return new SDSLVector::Iterator(*this);
  }


SDSLVector::SDSLVector(std::ifstream& file) :
  BitVector(file)
{
}

SDSLVector::SDSLVector(FILE* file) :
  BitVector(file)
{
}

SDSLVector::SDSLVector(Encoder& encoder, usint universe_size) :
  BitVector(encoder, universe_size)
{
    backing_vector = encoder.backing_vector;
}

SDSLVector::~SDSLVector()
{
}

//--------------------------------------------------------------------------

usint
SDSLVector::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

SDSLVector::Iterator::Iterator(const SDSLVector& par) :
  BitVector::Iterator(par)
{
    current = 0;
}

SDSLVector::Iterator::~Iterator()
{
}

static inline    bool min(usint i, int j)
{
    if (i < j) {
        return i;
    }else{
        return j;
    }
}
    
usint
SDSLVector::Iterator::rank(usint value, bool at_least)
{
    //FIXME: handle at_least arg
    usint ret =0;
    for(usint i = 0; i < min(value, backing_vector.size()); ++i)
    {
        if (backing_vector[i]) {
            ++ret;
            current = i;
    }
    return ret;
        
//   const SDSLVector& par = (const SDSLVector&)(this->parent);

//   if(value >= par.size) { return par.items; }

// //  this->valueLoop(value);

//   usint idx = this->sample.first + this->cur + 1;
//   if(!at_least && this->val > value)
//   {
//     idx--;
//   }
//   if(at_least && this->val < value)
//   {
//     this->getSample(this->block + 1);
//     this->run = 0;
//     idx = this->sample.first + this->cur + 1;
//   }
//   return idx;
}

usint
SDSLVector::Iterator::select(usint index)
{

    usint seen =0;
    for(usint i = 0; i < backing_vector.size(); ++i)
    {
        seen += backing_vector[i];
        if (seen == index) {
            current = i;
            return i;
        }
    }
    return backing_vector.size();
    
  // const SDSLVector& par = (const SDSLVector&)(this->parent);

  // if(index >= par.items) { return par.size; }
  // this->getSample(this->sampleForIndex(index));
  // this->run = 0;

  // usint lim = index - this->sample.first;
  // while(this->cur < lim)
  // {
  //   this->val += this->buffer.readDeltaCode();
  //   usint temp = this->buffer.readDeltaCode();
  //   this->val += temp - 1;
  //   this->cur += temp;
  // }
  // if(this->cur > lim)
  // {
  //   this->run = this->cur - lim;
  //   this->cur -= this->run;
  //   this->val -= this->run;
  // }

  // return this->val;
}

usint
SDSLVector::Iterator::selectNext()
{
    
    for(usint i = current + 1; i < backing_vector.size(); ++i)
    {
        if (backing_vector[i]) {
            current = i;
            return i;
        }
    }

    return backing_vector.size();
  // if(this->cur >= this->block_items)
  // {
  //   this->getSample(this->block + 1);
  //   this->run = 0;
  //   return this->val;
  // }

  // this->cur++;
  // if(this->run > 0)
  // {
  //   this->val++;
  //   this->run--;
  // }
  // else
  // {
  //   this->val += this->buffer.readDeltaCode();
  //   this->run = this->buffer.readDeltaCode() - 1;
  // }

  // return this->val;
}

// pair_type
// SDSLVector::Iterator::valueBefore(usint value)
// {
//   const SDSLVector& par = (const SDSLVector&)(this->parent);

//   if(value >= par.size) { return pair_type(par.size, par.items); }

//   this->getSample(this->sampleForValue(value));
//   if(this->val > value) { return pair_type(par.size, par.items); }
//   this->run = 0;

//   while(this->cur < this->block_items && this->val < value)
//   {
//     usint temp = this->buffer.readDeltaCode(value - this->val);
//     if(temp == 0) { break; }
//     this->val += temp;

//     temp = this->buffer.readDeltaCode();
//     this->cur += temp;
//     this->val += temp - 1;
//   }
//   if(this->val > value)
//   {
//     this->run = this->val - value;
//     this->val = value;
//     this->cur -= this->run;
//   }

//   return pair_type(this->val, this->sample.first + this->cur);
// }

// pair_type
// SDSLVector::Iterator::valueAfter(usint value)
// {
//   const SDSLVector& par = (const SDSLVector&)(this->parent);

//   if(value >= par.size) { return pair_type(par.size, par.items); }

//   this->valueLoop(value);

//   if(this->val < value)
//   {
//     this->getSample(this->block + 1);
//     this->run = 0;
//   }

//   return pair_type(this->val, this->sample.first + this->cur);
// }

// pair_type
// SDSLVector::Iterator::nextValue()
// {
//   if(this->cur >= this->block_items)
//   {
//     this->getSample(this->block + 1);
//     this->run = 0;
//     return pair_type(this->val, this->sample.first);
//   }

//   this->cur++;
//   if(this->run > 0)
//   {
//     this->val++;
//     this->run--;
//   }
//   else
//   {
//     this->val += this->buffer.readDeltaCode();
//     this->run = this->buffer.readDeltaCode() - 1;
//   }

//   return pair_type(this->val, this->sample.first + this->cur);
// }

// pair_type
// SDSLVector::Iterator::selectRun(usint index, usint max_length)
// {
//   usint value = this->select(index);

//   usint len = std::min(max_length, this->run);
//   this->run -= len; this->cur += len; this->val += len;

//   return pair_type(value, len);
// }

// pair_type
// SDSLVector::Iterator::selectNextRun(usint max_length)
// {
//   usint value = this->selectNext();

//   usint len = std::min(max_length, this->run);
//   this->run -= len; this->cur += len; this->val += len;

//   return pair_type(value, len);
// }



bool
SDSLVector::Iterator::isSet(usint value)
{

//  this->valueLoop(value);

    return backing_vector[value];
}

// usint
// SDSLVector::Iterator::countRuns()
// {
//   const SDSLVector& par = (const SDSLVector&)(this->parent);

//   if(par.items == 0) { return 0; }

//   usint runs = 1;
//   pair_type res = this->selectRun(0, par.items);
//   usint last = res.first + res.second;

//   while(last < par.size)
//   {
//     res = this->selectNextRun(par.items);
//     if(res.first < par.size && res.first > last + 1) { runs++; }
//     last = res.first + res.second;
//   }

//   return runs;
// }

// void
// SDSLVector::Iterator::valueLoop(usint value)
// {
//   this->getSample(this->sampleForValue(value));
//   this->run = 0;

//   if(this->val >= value) { return; }
//   while(this->cur < this->block_items)
//   {
//     this->val += this->buffer.readDeltaCode();
//     this->cur++;
//     this->run = this->buffer.readDeltaCode() - 1;
//     if(this->val >= value) { break; }

//     this->cur += this->run;
//     this->val += this->run;
//     if(this->val >= value)
//     {
//       this->run = this->val - value;
//       this->val = value;
//       this->cur -= this->run;
//       break;
//     }
//     this->run = 0;
//   }
// }

//--------------------------------------------------------------------------

// SDSLEncoder::SDSLEncoder(usint block_bytes, usint superblock_size) :
//   VectorEncoder(block_bytes, superblock_size),
//   run(EMPTY_PAIR)
// {
// }

// SDSLEncoder::~SDSLEncoder()
// {
// }

// void
// SDSLEncoder::setBit(usint value)
// {
//   this->setRun(value, 1);
// }

// void
// SDSLEncoder::setRun(usint start, usint len)
// {
//   if(this->items == 0)
//   {
//     this->setFirstBit(start);
//     if(len > 1)
//     {
//       this->SDSLncode(1, len - 1);
//     }
//     return;
//   }
//   if(start < this->size || len == 0) { return; }

//   // Write as much into the buffer as possible.
//   usint diff = start + 1 - this->size;
//   usint free_bits = this->buffer->bitsLeft();
//   usint code_bits = this->buffer->deltaCodeLength(diff);
//   if(free_bits > code_bits) // At least a part of the run fits into the block.
//   {
//     free_bits -= code_bits;
//     usint run_bits = this->buffer->deltaCodeLength(len);
//     if(run_bits <= free_bits)
//     {
//       this->SDSLncode(diff, len);
//       return;
//     }

//     // Encode as much as possible and let the rest spill.
//     usint llen = 1;
//     while(llen + 2 * length(llen + 1) - 1 <= free_bits) { llen++; }
//     llen = ((usint)1 << llen) - 1;

//     this->SDSLncode(diff, llen);
//     len -= llen;

//     // A new sample will be added.
//     this->size++;
//     this->items++;
//   }
//   else
//   {
//     this->size = start + 1;
//     this->items++;
//   }

//   // Didn't fit into the block. A new sample & block required.
//   this->addNewBlock();
//   if(len > 1)
//   {
//     this->SDSLncode(1, len - 1);
//   }
// }

// void
// SDSLEncoder::addBit(usint value)
// {
//   this->addRun(value, 1);
// }

// void
// SDSLEncoder::addRun(usint start, usint len)
// {
//   if(this->run.second == 0)
//   {
//     this->run = pair_type(start, len);
//   }
//   else if(start == this->run.first + this->run.second)
//   {
//     this->run.second += len;
//   }
//   else
//   {
//     this->setRun(this->run.first, this->run.second);
//     this->run = pair_type(start, len);
//   }
// }

// void
// SDSLEncoder::flush()
// {
//   this->setRun(this->run.first, this->run.second);
//   this->run.second = 0;
// }
    
SDSLEncoder::SDSLEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size)
{
}

SDSLEncoder::~SDSLEncoder()
{
}

void
SDSLEncoder::setBit(usint value)
{
    backing_vector[value] = 1;
}

void
SDSLEncoder::setRun(usint start, usint len)
{
  for(usint i = start; i < start + len; i++) { this->setBit(i); }
}

void
SDSLEncoder::addBit(usint value)
{
  this->setBit(value);
}

void
SDSLEncoder::addRun(usint start, usint len)
{
  this->setRun(start, len);
}

void
SDSLEncoder::flush()
{
}

} // namespace CSA


} // namespace CSA
