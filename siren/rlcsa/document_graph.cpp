#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stack>

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif

#include "rlcsa.h"
#include "docarray.h"

using namespace CSA;


/*
  This is a temporary program for experimenting with subset compression in document listing.
*/


int testDocarray(const std::string& base_name);
int buildIndex(const std::string& base_name, const std::string& prefix);


int
main(int argc, char** argv)
{
  std::cout << "Document graph extractor" << std::endl;
  if(argc < 2)
  {
    std::cout << "Usage:" << std::endl;
    std::cout << "document_graph basename block_size storing_factor" << std::endl;
    std::cout << "  to build the document graph and sparse suffix tree" << std::endl;
    std::cout << "document_graph basename prefix" << std::endl;
    std::cout << "  to compress the grammar and build the index" << std::endl;
    std::cout << "document_graph basename" << std::endl;
    std::cout << "  to test the document array" << std::endl;
    return 1;
  }
  std::cout << std::endl;

  std::string base_name = argv[1];
  std::cout << "Base name: " << base_name << std::endl;

  // Load the index.
  if(argc == 2)
  {
    std::cout << std::endl;
    return testDocarray(base_name);
  }

  // Load the index and compress
  if(argc == 3)
  {
    std::string prefix = argv[2];
    std::cout << "Prefix: " << prefix << std::endl;
    std::cout << std::endl;
    return buildIndex(base_name, prefix);
  }

  usint block_size = atoi(argv[2]);
  std::cout << "Block size: " << block_size << std::endl;
  usint storing_factor = atoi(argv[3]);
  std::cout << "Storing factor: " << storing_factor << std::endl;
  std::cout << std::endl;

  RLCSA rlcsa(argv[1]);
  if(!rlcsa.isOk() || !rlcsa.supportsLocate())
  {
    std::cerr << "Error: Locate is not supported!" << std::endl;
    return 2;
  }
  rlcsa.printInfo();
  rlcsa.reportSize(true);
  if(block_size == 0 || storing_factor == 0) { return 0; }


  double start = readTimer();

  // Get the suffix array from RLCSA.
  usint size = rlcsa.getSize();
  uint docs = rlcsa.getNumberOfSequences();
  usint* documents = new usint[size];
  #pragma omp parallel for schedule(static)
  for(usint i = 0; i < size; i += MEGABYTE)
  {
    rlcsa.locate(pair_type(i, std::min(i + MEGABYTE, size - 1)), documents + i);
  }
  double locate = readTimer();
  std::cout << "Locate: " << (locate - start) << " seconds" << std::endl;

  // Get the LCP array.
  uint* lcp = new uint[size + 1]; lcp[size] = 0;
/*  #ifdef MULTITHREAD_SUPPORT
  usint threads = omp_get_max_threads();
  PLCPVector* plcpvec = rlcsa.buildPLCP(16, threads);
  std::cout << "Total: " << plcpvec->getNumberOfItems() << " / " << plcpvec->getSize() << std::endl;
  PLCPVector::Iterator* iters[threads];
  for(usint i = 0; i < threads; i++) { iters[i] = new PLCPVector::Iterator(*plcpvec); }
  #pragma omp parallel for schedule(static, MEGABYTE)
  for(usint i = 0; i < size; i++)
  {
    lcp[i] = iters[omp_get_thread_num()]->select(documents[i]) - 2 * documents[i];
  }
  for(usint i = 0; i < threads; i++) { delete iters[i]; iters[i] = 0; }
  #else*/
  PLCPVector* plcpvec = rlcsa.buildPLCP(16);
  PLCPVector::Iterator iter(*plcpvec);
  for(usint i = 0; i < size; i++)
  {
    lcp[i] = iter.select(documents[i]) - 2 * documents[i];
  }
//  #endif
  delete plcpvec;
  double lcparray = readTimer();
  std::cout << "LCP: " << (lcparray - locate) << " seconds" << std::endl;

  // Convert SA to document array.
  #pragma omp parallel for schedule(static)
  for(usint i = 0; i < size; i += MEGABYTE)
  {
    rlcsa.getSequenceForPosition(documents + i, std::min(MEGABYTE, size - i));
  }
  // Change document array to be 1-based.
  #pragma omp parallel for schedule(static, MEGABYTE)
  for(usint i = 0; i < size; i++) { documents[i]++; }
  double find = readTimer();
  std::cout << "Documents: " << (find - lcparray) << " seconds" << std::endl;

  // Build a truncated suffix tree.
  // A node with at most block_size suffixes should not have any children.
  std::stack<STNode*> nodestack; nodestack.push(new STNode(0, pair_type(0, 0)));
  STNode* prev = 0; STNode* root = 0;
  for(usint i = 1; i <= size; i++)
  {
    usint left = i - 1;
    while(lcp[i] < nodestack.top()->string_depth)
    {
      nodestack.top()->range.second = i - 1;
      prev = nodestack.top(); nodestack.pop();
      if(length(prev->range) <= block_size) { prev->deleteChildren(); }
      root = prev;  // Last processed node.
      left = prev->range.first;
      if(lcp[i] <= nodestack.top()->string_depth) { nodestack.top()->addChild(prev); prev = 0; }
    }
    if(lcp[i] > nodestack.top()->string_depth)
    {
      STNode* curr = new STNode(lcp[i], pair_type(left, left));
      if(prev != 0) { curr->addChild(prev); prev = 0; }
      nodestack.push(curr);
    }
  }
  while(root->parent != 0) { root = root->parent; }  // Find the actual root.
  root->range.second = size - 1; root->containsAllDocuments();
  delete[] lcp; lcp = 0;
  uint st_nodes = 0, leaf_nodes = 0;
  root->addLeaves();
  root->verifyTree();
  root->determineSize(st_nodes, leaf_nodes);
  double st = readTimer();
  std::cout << "Suffix tree: " << (st - find) << " seconds" << std::endl;
  std::cout << "Total " << size << " positions, " << docs << " documents" << std::endl;
  std::cout << "Truncated suffix tree has " << st_nodes << " nodes" << std::endl;
  std::cout << std::endl;


  // Start processing.
  uint tree_leaves = 0, tree_nodes = leaf_nodes, graph_nodes = 0, edges = 0;
  std::string graph_name = base_name + ".graph";
  std::ofstream output(graph_name.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error: Cannot open output file " << graph_name << std::endl;
    delete[] documents;
    return 3;
  }
  std::string singleton_name = base_name + ".singletons";
  std::ofstream singletons(singleton_name.c_str(), std::ios_base::binary);
  if(!singletons)
  {
    std::cerr << "Error: Cannot open output file " << singleton_name << std::endl;
    delete[] documents;
    output.close();
    return 3;
  }
  output.write((char*)&tree_nodes, sizeof(tree_nodes)); // Placeholder for nodes/subsets.
  output.write((char*)&edges, sizeof(edges)); // Placeholder for edges.

  // For each suffix tree node, produce the corresponding subset if necessary.
  double sets = readTimer();
  while(!nodestack.empty()) { nodestack.pop(); }
  for(prev = root; prev != 0; prev = prev->child) { nodestack.push(prev); }
  while(nodestack.top() != root)
  {
    STNode* curr = nodestack.top(); nodestack.pop();
    for(STNode* temp = curr->sibling; temp != 0; temp = temp->child) { nodestack.push(temp); }
    if(curr->contains_all)
    {
      curr->id = tree_nodes; tree_nodes++;
      continue;
    }

    curr->computeStoredDocuments(documents);
    if(curr->docs->size() >= docs) { curr->containsAllDocuments(); }

    if(curr->child == 0 || curr->contains_all || curr->stored_documents > storing_factor * curr->docs->size())
    {
      curr->storeThisSet();
      int curr_node = -1 * docs;
      if(curr->child == 0) { curr->id = tree_leaves; tree_leaves++; curr_node -= tree_leaves; }
      else                 { curr->id = tree_nodes; tree_nodes++; curr_node -= tree_nodes; }
      if(!(curr->contains_all))
      {
        // Store singletons in a separate file. The presence of singletons greatly slows down the
        // search for bicliques.
        if(curr->docs->size() == 1)
        {
          singletons << (-curr_node) << ": " << curr->docs->at(0) << std::endl;
        }
        else
        {
          output.write((char*)&curr_node, sizeof(curr_node));
          for(std::vector<uint>::iterator iter = curr->docs->begin(); iter != curr->docs->end(); ++iter)
          {
            uint temp = *iter;
            output.write((char*)&temp, sizeof(temp));
          }
          graph_nodes++;
          edges += curr->docs->size();
        }
      }
    }
    else
    {
      curr->remove();
      delete curr;
    }
  }
  root->id = tree_nodes; tree_nodes++;
  double block_stop = readTimer();
  std::cout << "Blocks: " << (block_stop - sets) << " seconds" << std::endl;
  std::cout << std::endl;

  // Write the number of nodes & edges.
  output.seekp(0, std::ios::beg);
  output.write((char*)&graph_nodes, sizeof(graph_nodes));
  output.write((char*)&edges, sizeof(edges));
  std::cout << "Graph has " << graph_nodes << " nodes, " << edges << " edges" << std::endl;
  std::cout << std::endl;
  output.close();
  singletons.close();
  delete[] documents;

  // Write the tree structure.
  DocArray docarray(root, rlcsa);
  if(docarray.isOk())
  {
    std::cout << "Final tree has " << docarray.getNumberOfNodes() << " nodes, " << tree_leaves << " leaves" << std::endl;
    std::cout << std::endl;
    docarray.reportSize(true);
    docarray.writeTo(base_name);
  }

  delete root;
  return 0;
}


int
testDocarray(const std::string& base_name)
{
  RLCSA rlcsa(base_name);
  if(!(rlcsa.isOk())) { return 4; }
  DocArray docarray(rlcsa, base_name, true);
  if(!(docarray.isOk())) { return 5; }

  rlcsa.printInfo();
  rlcsa.reportSize(true);
  docarray.reportSize(true);

  return 0;
}


int
buildIndex(const std::string& base_name, const std::string& prefix)
{
  RLCSA rlcsa(base_name);
  if(!(rlcsa.isOk())) { return 4; }
  DocArray docarray(rlcsa, base_name, false);
  if(!(docarray.isOk())) { return 5; }

  docarray.readRules(prefix, true);
  if(!(docarray.hasGrammar())) { return 6; }
  docarray.writeTo(base_name);

  rlcsa.printInfo();
  rlcsa.reportSize(true);
  docarray.reportSize(true);

  return 0;
}
