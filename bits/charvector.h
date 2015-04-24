#ifndef CHARVECTOR_H
#define CHARVECTOR_H
#include <fstream>
#include "../misc/definitions.h"
#include "rlevector.h"
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
    void syncFMIndex();
    void writeTo(std::ofstream& file) const;
    void load(std::ifstream& file);
    void writeTo(FILE* file) const;
    usint reportSize() const;
//    usint rank(usint c, usint num, bool at_least = false) const ;
    std::vector<usint> restricted_unique_range_values(usint l, usint r, usint min, usint max) const;//FIXME: don't copy vector in return  
    usint maxlength() const;
    void constructF(CSA::RLEEncoder &inedges, unsigned int incomingedge_offset);
    void setwt(sdsl::int_vector<> &v);
    inline std::map<usint, CSA::BitVector*>::const_iterator begin() const { return array.begin();}
    inline std::map<usint, CSA::BitVector*>::const_iterator end() const { return array.end();}
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
    sdsl::wt_int<> *wt;
    std::map<usint, CSA::BitVector*> array;
    // sdsl::csa_wt<sdsl::wt_int<>, 
    //              64, 
    //              64, 
    //              sdsl::sa_order_sa_sampling<>, 
    //              sdsl::int_vector<>, 
    //              sdsl::int_alphabet<>
    //              > fm_index;

    CSA::BitVector *incoming; // F?
    CSA::BitVector::Iterator *incoming_itr;
};

} // namespace CSA


#endif // CHARVECTOR_H
