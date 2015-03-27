#ifndef CHARVECTOR_H
#define CHARVECTOR_H
#include <fstream>
#include "../misc/definitions.h"
#include "bitvector.h"
#include <map>
#include "deltavector.h"
#include <sdsl/suffix_arrays.hpp>

namespace CSA
{

class CharVector
{
public:
    void populate(usint c, DeltaVector::Encoder*, usint offset);
    void populate(usint c, DeltaVector *);
    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;
    usint reportSize() const;
//    usint rank(usint c, usint num, bool at_least = false) const ;
    usint maxlength() const;
    class Iterator {
    public:
        Iterator(BitVector::Iterator *itr);
        ~Iterator();
        usint select(usint);
        usint selectNext();
        usint rank(usint, bool at_least = false);
        bool isSet(usint);

    private:
        BitVector::Iterator *itr;
    };

    Iterator *newIterator(usint c) const;
private:
    std::map<usint, CSA::BitVector*> array;
    sdsl::csa_wt<sdsl::wt_int<>, 
                 64, 
                 64, 
                 sdsl::sa_order_sa_sampling<>, 
                 sdsl::int_vector<>, 
                 sdsl::int_alphabet<>
                 > fm_index;


};

} // namespace CSA


#endif // CHARVECTOR_H
