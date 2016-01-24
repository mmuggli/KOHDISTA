#ifndef OMFIO_H
#define OMFIO_H

#include <iostream>
#include <vector>
class Omfio {
public:
    Omfio(std::string fname);
    void dump();

private:
    class Rmap {
    public:
        Rmap(std::string, std::string, std::string, std::vector<float> fragments);
        void dump();
    private:        
        std::string id;
        std::string enzyme_name;
        std::string enzyme_acronym;
        std::vector<float> fragments;

    };

    std::vector<Rmap> rmaps;
};

#endif
