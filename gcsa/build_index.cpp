#include <cstdlib>
#include <fstream>
#include <iostream>

#include "gcsa.h"
#include "graph.h"

#include <misc/utils.h>


typedef CSA::usint usint;


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
    CSA::start_time = start;

    GCSA::Graph* graph = new GCSA::Graph(base_name);
    double loaded_time = CSA::readTimer();
    std::cout << "Graph loaded in " << loaded_time - start << " seconds." << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    if(!graph->ok) { return 2; }
    double backboned_time = loaded_time;
    if(backbone)
    {
        std::cout << "Generating backbone... "; std::cout.flush();
        graph->createBackbone();
        backboned_time = CSA::readTimer();
		std::cout << "done in " << backboned_time - loaded_time << " seconds." << std::endl << std::endl;
        std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
	}
    graph->printInfo();

    std::cout << "Constructing pathgraph..."  ;
    GCSA::PathGraph* pg = new GCSA::PathGraph(*graph);
    double pg_time = CSA::readTimer();
    std::cout << "done in " << pg_time - backboned_time  << " seconds." << std::endl << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    delete graph; graph = 0;
    pg->printInfo();
    while(pg->status != GCSA::PathGraph::sorted)
    {
        std::cout << "Starting next generation of pathgraph" << std::endl;
        std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
        if(pg->status != GCSA::PathGraph::ok)
        {
            std::cerr << "Error: Invalid PathGraph!" << std::endl;
            delete pg; pg = 0;
            return 3;
        }
        GCSA::PathGraph* next = new GCSA::PathGraph(*pg);
        delete pg; pg = next;
        pg->printInfo();
        std::cout << "done with construction generation" << std::endl;
        std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    }
    std::cout << std::endl;
    std::cout << "constructing new graph" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    graph = new GCSA::Graph(base_name);
    if(backbone) { 
        std::cout << "constructing backbone" << std::endl;
        std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
        graph->createBackbone(); 
        std::cout << "done constructing backbone" << std::endl;
        std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    }
    std::cout << "constructing gcsa" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    GCSA::GCSA gcsa(*pg, *graph, true);
    std::cout << "writing gcsa" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    gcsa.writeTo(base_name);
    //delete graph; graph = 0;
    //delete pg; pg = 0;


    double time = CSA::readTimer() - start;
    std::cout << "Used " << time << " seconds." << std::endl;
    std::cout << "Memory: " << CSA::memoryUsage() << " kB" << std::endl;
    std::cout << std::endl;

    return 0;
}
