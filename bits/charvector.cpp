#include "charvector.h"

namespace CSA
{
void CharVector::writeTo(std::ofstream& file) const
{
}

void CharVector::writeTo(FILE* file) const{}
usint CharVector::reportSize() const{return 0;}
//     usint CharVector::rank(usint c, usint num, bool at_least) const {
//     // BitVector::Iterator* array_iter = this->array.at(c)->newIterator();
//     // usint r = array_iter->rank(range.first, at_least);
//     // delete array_iter;
//     // return r;
//     return 1;
// }

usint CharVector::Iterator::select(usint){return 0;}
usint  CharVector::Iterator::selectNext(){return 0;}
    usint  CharVector::Iterator::rank(usint, bool ){return 0;}
bool  CharVector::Iterator::isSet(usint){return 0;}


    CharVector::Iterator *CharVector::newIterator(usint c) const{return NULL;}
} //namespace CSA
