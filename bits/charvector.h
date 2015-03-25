#ifndef CHARVECTOR_H
#define CHARVECTOR_H
#include <fstream>
#include "../misc/definitions.h"

namespace CSA
{

class CharVector
{
public:
    
    void writeTo(std::ofstream& file) const;
    void writeTo(FILE* file) const;
    usint reportSize() const;
//    usint rank(usint c, usint num, bool at_least = false) const ;

    class Iterator {
    public:
        usint select(usint);
        usint selectNext();
        usint rank(usint, bool at_least = false);
        bool isSet(usint);
};

    Iterator *newIterator(usint c) const;

};

} // namespace CSA


#endif // CHARVECTOR_H
