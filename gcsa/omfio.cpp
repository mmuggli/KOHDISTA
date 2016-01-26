#include "omfio.h"
#include <fstream>
Omfio::Omfio(std::string fname)
{
    std::cout << "Loading Rmaps from " << fname << "..." << std::endl;
    std::ifstream rmap_file(fname.c_str());
    std::string id, enzyme_name, enzyme_acronym;

    while (rmap_file >> id) {
        std::vector<long unsigned int> fragments;
        rmap_file >> enzyme_name >> enzyme_acronym;
        float fragment;
        while (rmap_file >> fragment) {
            fragments.push_back((long unsigned int)(fragment*1000.0));
        }
        rmap_file.clear();
        Rmap rmap(id, enzyme_name, enzyme_acronym, fragments);
        rmaps.push_back(rmap);
    }
    std::cout << "Loaded " << rmaps.size() << " Rmaps." << std::endl;

}
void Omfio::dump()
{
    for (std::vector<Omfio::Rmap>::iterator it = rmaps.begin(); it != rmaps.end(); ++it) {
        it->dump();
    }
}

Omfio::Rmap::Rmap(std::string a_id, std::string a_enzn, std::string a_enza, std::vector<long unsigned int> a_fragments) : id(a_id), enzyme_name(a_enzn), enzyme_acronym(a_enza), fragments(a_fragments)
{
}

void Omfio::Rmap::dump()
{
    std::cout << id << std::endl;
    std::cout << "\t" << enzyme_name << "\t" << enzyme_acronym;
    for (std::vector<long unsigned int>::iterator it = fragments.begin(); it != fragments.end(); ++it) {
        std::cout << "\t" << ((float)*it)/1000.0;
    }
    std::cout << std::endl << std::endl;
}

