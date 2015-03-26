#include "charvector.h"
#include "deltavector.h"

namespace CSA
{

    void CharVector::populate(usint c, DeltaVector::Encoder* encoder, usint offset)
    {

        array[c] = new CSA::DeltaVector(*encoder, offset);
    }
    void CharVector::writeTo(std::ofstream& file) const
    {
        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        {
            if(array.count(i)) { array.at(i)->writeTo(file); }
        }

    }

    void CharVector::writeTo(FILE* file) const
    {
        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        {
            if(array.count(i)) { array.at(i)->writeTo(file); }
        }
    }
    usint CharVector::reportSize() const
    {
        usint array_size = 0;
        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        {
            if(array.count(i)) { array_size += array.at(i)->reportSize(); }
        }
        return array_size;


    }
//     usint CharVector::rank(usint c, usint num, bool at_least) const {
//     // BitVector::Iterator* array_iter = this->array.at(c)->newIterator();
//     // usint r = array_iter->rank(range.first, at_least);
//     // delete array_iter;
//     // return r;
//     return 1;
// }

    usint CharVector::Iterator::select(usint i){return itr->select(i);}
    usint  CharVector::Iterator::selectNext(){return itr->selectNext();}
    usint  CharVector::Iterator::rank(usint i, bool at_least){return itr->rank(i, at_least);}
    bool  CharVector::Iterator::isSet(usint i){return itr->isSet(i);}


    CharVector::Iterator *CharVector::newIterator(usint c) const{return new Iterator(array.at(c)->newIterator());}

    CharVector::Iterator::Iterator(BitVector::Iterator *itr) : itr(itr)
    {
    }
    CharVector::Iterator::~Iterator() 
    {
        delete itr;
    }
        
} //namespace CSA
