#include <cstdlib>
#include <fstream>
#include <iostream>

#include "gcsa.h"
#include "graph.h"

#include <misc/utils.h>


//using namespace CSA;
typedef CSA::usint usint;

//typedef CSA::pair_type pair_type;

int
main(int argc, char** argv)
{
  std::cout << "GCSA builder" << std::endl;
  std::cout << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage: build_index [-b] base_name" << std::endl;
    	std::cout << "  -b  create backbone information" << std::endl;
    return 1;
  }

	bool backbone = false;
	usint name_arg = 1;
	if(argc >= 3 && argv[1][0] == '-' && argv[1][1] == 'b')
	{
		backbone = true; name_arg = 2;
	}

    std::string base_name(argv[name_arg]);
    std::cout << "Input:  " << base_name << std::endl;
    std::cout << std::endl;

    double start = CSA::readTimer();

    GCSA::Graph* graph = new GCSA::Graph(base_name);
    double loaded_time = CSA::readTimer();
    std::cout << "Graph loaded in " << loaded_time - start << " seconds." << std::endl;
    if(!graph->ok) { return 2; }
    double backboned_time = loaded_time;
    if(backbone)
    {
        std::cout << "Generating backbone... "; std::cout.flush();
        graph->createBackbone();
        backboned_time = CSA::readTimer();
		std::cout << "done in " << backboned_time - loaded_time << " seconds." << std::endl << std::endl;
	}
    graph->printInfo();

    std::cout << "Constructing pathgraph..."  ;
    GCSA::PathGraph* pg = new GCSA::PathGraph(*graph);
    double pg_time = CSA::readTimer();
    std::cout << "done in " << pg_time - backboned_time  << " seconds." << std::endl << std::endl;
    delete graph; graph = 0;
    pg->printInfo();
    while(pg->status != GCSA::PathGraph::sorted)
    {
        if(pg->status != GCSA::PathGraph::ok)
        {
            std::cerr << "Error: Invalid PathGraph!" << std::endl;
            delete pg; pg = 0;
            return 3;
        }
        GCSA::PathGraph* next = new GCSA::PathGraph(*pg);
        delete pg; pg = next;
        pg->printInfo();
    }
    std::cout << std::endl;

    graph = new GCSA::Graph(base_name);
    if(backbone) { graph->createBackbone(); }
    GCSA::GCSA gcsa(*pg, *graph, true);
    gcsa.writeTo(base_name);
    //delete graph; graph = 0;
    //delete pg; pg = 0;


    double time = CSA::readTimer() - start;
    std::cout << "Used " << time << " seconds." << std::endl;
    std::cout << "Memory: " << CSA::memoryUsage() << " kB" << std::endl;
    std::cout << std::endl;

    return 0;
}
