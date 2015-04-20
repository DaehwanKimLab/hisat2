#include <cstdlib>
#include <fstream>
#include <iostream>

#include "graph.h"

#include <misc/utils.h>


using namespace CSA;


int
main(int argc, char** argv)
{
  std::cout << "Determinize and merge automata" << std::endl;
  std::cout << std::endl;
  if(argc < 3)
  {
    std::cout << "Usage: determinize [-b|-r] automaton1 [automaton2 [..]] base_name" << std::endl;
    std::cout << "  -b  Use the backbone information stored in the automaton." << std::endl;
    std::cout << "  -r  Discard the backbone information stored in the automaton." << std::endl;
    return 1;
  }

  pair_type range(1, argc - 2);
  bool use_backbone = false, remove_backbone = false;
  if(argv[1][0] == '-')
  {
    switch(argv[1][1])
    {
      case 'b':
        range.first++;
        use_backbone = true;
        std::cout << "Using backbone information." << std::endl;
        break;
      case 'r':
        range.first++;
        remove_backbone = true;
        std::cout << "Removing backbone information." << std::endl;
        break;
      default:
        std::cout << "Invalid option: " << argv[1] << std::endl << std::endl;
        return 2;
    }
  }

  std::cout << "Number of automata: " << length(range) << std::endl;
  std::string base_name = argv[argc - 1];
  std::cout << "Base name: " << base_name << std::endl;
  std::cout << std::endl;

  double start = readTimer();
  Graph* result = 0;
  for(usint arg = range.first; arg <= range.second; arg++)
  {
    std::string automaton_name(argv[arg]);
    std::cout << "Automaton: " << automaton_name << std::endl;
    Graph* graph = new Graph(automaton_name, true);
    if(!(graph->ok)) { delete graph; std::cout << std::endl; continue; }
    graph->printInfo();
    if(remove_backbone) { graph->removeBackbone(); }
    if(graph->isDeterministic())
    {
      std::cout << "The automaton is already deterministic." << std::endl;
    }
    else
    {
      std::cout << "Determinizing..."; std::cout.flush();
      if(use_backbone)    { graph->createBackbone(); }
      graph->determinize();
      if(use_backbone)    { graph->encodeBackbone(); }
      std::cout << " done!" << std::endl;
      graph->printInfo();
    }
    if(result == 0) { result = graph; }
    else            { result->addGraph(*graph); delete graph; }
    std::cout << std::endl;
  }
  if(result != 0)
  {
    std::cout << "Final graph:" << std::endl;
    result->printInfo();
    std::cout << std::endl;
    result->write(base_name);
    delete result;
  }
  double stop = readTimer();

  std::cout << "Used " << (stop - start) << " seconds." << std::endl;
  std::cout << "Memory: " << memoryUsage() << " kB" << std::endl;
  std::cout << std::endl;

  return 0;
}
