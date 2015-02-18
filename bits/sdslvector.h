#ifndef SDSLVECTOR_H
#define SDSLVECTOR_H

#include <fstream>

#include "bitvector.h"


namespace CSA
{

class SDSLVector : public BitVector
{
  public:
    typedef VectorEncoder Encoder;
    BitVector::Iterator* newIterator();
    explicit SDSLVector(std::ifstream& file);
    explicit SDSLVector(FILE* file);
    SDSLVector(Encoder& encoder, usint universe_size);


    ~SDSLVector();

//--------------------------------------------------------------------------

    usint reportSize() const;

//--------------------------------------------------------------------------

    class Iterator : public BitVector::Iterator
    {
      public:
        explicit Iterator(const SDSLVector& par);
        ~Iterator();

        usint rank(usint value, bool at_least = false);

        usint select(usint index);
        usint selectNext();

        // pair_type valueBefore(usint value);
        // pair_type valueAfter(usint value);
        // pair_type nextValue();

        // pair_type selectRun(usint index, usint max_length);
        // pair_type selectNextRun(usint max_length);

        bool isSet(usint value);

        // usint countRuns();

      protected:

//        void valueLoop(usint value);

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------

  protected:

    // These are not allowed.
    SDSLVector();
    SDSLVector(const SDSLVector&);
    SDSLVector& operator = (const SDSLVector&);
};


} // namespace CSA


#endif // SDSLVECTOR_H

