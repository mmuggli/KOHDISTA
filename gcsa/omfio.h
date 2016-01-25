#ifndef OMFIO_H
#define OMFIO_H

#include <iostream>
#include <vector>

class Omfio {
public:
    class Rmap {
    public:
        Rmap(std::string, std::string, std::string, std::vector<long unsigned int> fragments);
        void dump();

        std::string id;
        std::string enzyme_name;
        std::string enzyme_acronym;
        std::vector<long unsigned int> fragments;

    };


    Omfio(std::string fname);
    void dump();




    std::vector<Rmap> rmaps;
};

#endif
