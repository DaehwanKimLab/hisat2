#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include <rlcsa.h>

#include "bwasearch.h"
#include "parameter_handler.h"
#include "pattern_classifier.h"


using namespace CSA;


int main(int argc, char** argv)
{
  std::cout << "RLCSA test" << std::endl;
  std::cout << std::endl;

  ParameterHandler handler(argc, argv, true, "Usage: rlcsa_test [options] base_name [patterns]");
  if(!(handler.ok)) { handler.printUsage(); return 1; }
  handler.printOptions();

  const RLCSA rlcsa(handler.index_name);
  if(!rlcsa.isOk()) { return 2; }
  rlcsa.reportSize(true);
  if(handler.patterns_name == 0) { return 0; }
  BWASearch<RLCSA> bwasearch(rlcsa);

  std::ifstream patterns(handler.patterns_name, std::ios_base::binary);
  if(!patterns)
  {
    std::cerr << "Error opening pattern file!" << std::endl;
    return 3;
  }
  PatternClassifier classifier(handler.write ? handler.patterns_name : "");

  std::vector<std::string> rows;
  readRows(patterns, rows, true);

  usint total = 0, n = rows.size();
  usint found = 0, forward = 0, reverse = 0;
  usint forward_occurrences = 0, reverse_occurrences = 0;

  usint* edit_distances = new usint[handler.k + 1];
  for(usint i = 0; i <= handler.k; i++) { edit_distances[i] = 0; }
  usint* match_counts = new usint[n];
  for(usint i = 0; i < n; i++) { match_counts[i] = 0; }

  double start = readTimer();
  for(usint i = 0; i < n; i++)
  {
    total += rows[i].length();
    bool match = false;

    // Always start with exact matching.
    pair_type result = bwasearch.find(rows[i], false, handler.skip);
    usint temp = bwasearch.handleOccurrences(result, handler.locate, handler.max_matches);
    if(temp > 0)
    {
      forward_occurrences += temp; match_counts[i] += temp;
      found++; forward++; match = true;
      edit_distances[0]++;
      if(handler.write) { classifier.forward(rows[i]); }
    }
    if(handler.reverse_complement)
    {
      pair_type reverse_result = bwasearch.find(rows[i], true, handler.skip);
      temp = bwasearch.handleOccurrences(reverse_result, handler.locate, handler.max_matches);
      if(temp > 0)
      {
        if(!match)         { found++; edit_distances[0]++; }
        if(handler.k == 0) { reverse_occurrences += temp; }
        else               { forward_occurrences += temp; }
        match_counts[i] += temp;
        if(handler.max_matches > 0 && match_counts[i] > handler.max_matches)
        {
          match_counts[i] = handler.max_matches + 1;
        }
        reverse++; match = true;
        if(handler.write)  { classifier.reverse(rows[i]); }
      }
    }
    if(handler.write && !match) { classifier.notfound(rows[i]); }

    // Do approximate matching only if there are no exact matches.
    if(handler.k > 0 && !match)
    {
      std::vector<MatchInfo*>* results = bwasearch.find(rows[i], handler.k, handler.indels, handler.penalties);
      temp = bwasearch.handleOccurrences(results, handler.locate, handler.max_matches);
      if(temp > 0)
      {
        forward_occurrences += temp; match_counts[i] += temp;
        found++;
        edit_distances[(*(results->begin()))->errors]++;
      }
      bwasearch.deleteResults(results);
    }
  }
  double time = readTimer() - start;
  double megabases = total / (double)MILLION;

  std::cout << "Patterns:     " << n << " (" << (n / time) << " / sec)" << std::endl;

  std::cout << "Total size:   " << megabases << " megabases";
  if(!(handler.locate))
  {
    std::cout << " (" << (megabases / time) << " megabases / sec)";
  }
  std::cout << std::endl;

  std::cout << "Found:        " << found;
  if(handler.k == 0 && handler.reverse_complement)
  {
    std::cout << " (" << forward << " forward, " << reverse << " reverse complement)";
  }
  std::cout << std::endl;

  if(handler.locate)
  {
    std::cout << "Occurrences:  " << (forward_occurrences + reverse_occurrences) << " (";
    if(handler.k == 0 && handler.reverse_complement)
    {
      std::cout << forward_occurrences << " forward, " << reverse_occurrences << " reverse complement, ";
    }
    std::cout << ((forward_occurrences + reverse_occurrences) / time) << " / sec)" << std::endl;
  }

  std::cout << "Time:         " << time << " seconds" << std::endl;
  std::cout << std::endl;

  if(handler.verbose)
  {
    std::cout << "MATCH STATISTICS" << std::endl << std::endl;
    std::cout << "Edit distances for matching patterns" << std::endl;
    for(usint i = 0; i <= handler.k; i++) { std::cout << "  " << i << "  " << edit_distances[i] << std::endl; }
    std::cout << std::endl;

    std::sort(match_counts, match_counts + n);
    std::cout << "Number of matches for matching patterns" << std::endl;
    usint curr = 0, count = 0;
    for(usint i = 0; i < n; i++)
    {
      if(match_counts[i] != curr)
      {
        if(curr > 0) { std::cout << "  " << curr << "  " << count << std::endl; }
        curr = match_counts[i]; count = 1;
      }
      else
      {
        count++;
      }
    }
    if(curr > 0) { std::cout << "  " << curr << "  " << count << std::endl; }
    std::cout << std::endl;
  }
  delete edit_distances; edit_distances = 0;
  delete match_counts; match_counts = 0;

  return 0;
}
