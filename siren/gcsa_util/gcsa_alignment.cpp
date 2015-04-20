#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include <misc/utils.h>

#include "gcsa.h"
#include "bwasearch.h"
#include "parameter_handler.h"
#include "pattern_classifier.h"

/**
 * Parts taken and modified from the GCSA package
Copyright (c) 2007-2011 Jouni Siren

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/

using namespace CSA;

usint locateBackbone(usint index, const GCSA& gcsa, const Backbone& backb)
{
    usint pos = gcsa.locate(index);
    while (pos > backb.getSize())
    {
        // index = gcsa.Psi(index);
        std::vector<usint>* tmp = gcsa.getSuccessors(index);
        index = (*tmp)[0]; // First outgoing edge
        delete tmp;
        pos = gcsa.locate(index);
    }
    return pos; // Returns the first backbone location found on the left-side of the occurrence.
}

usint outputOccurrences(usint no, const GCSA& gcsa, const Backbone& backb, const std::string& pattern, std::vector<MatchInfo*>* results, usint max_matches)
{
    if (results->empty()) { return 0; }
    usint temp = 0;

    for(std::vector<MatchInfo*>::const_iterator iter = results->begin(); iter != results->end(); ++iter)
    {
        pair_type range = (*iter)->range;
        gcsa.convertToSARange(range);
        for (usint i = range.first; i <= range.second; ++i)
        {
            std::cout << no+1 << '\t'
                      << pattern << '\t' 
                      << locateBackbone(i, gcsa, backb) << '\t' << std::endl;
        }
        temp += length(range);
    }
    if(max_matches > 0 && temp > max_matches) { return max_matches + 1; }
    return temp;
}

int main(int argc, char** argv)
{
  std::cerr << "GCSA alignment test" << std::endl;
  std::cerr << std::endl;

  ParameterHandler handler(argc, argv, false, "Usage: gcsa_alignment [options] base_name [patterns]");
  if(!(handler.ok)) { handler.printUsage(); std::cerr << "error: check stdout output for more details" << std::endl; return 1; }

  const GCSA gcsa(handler.index_name);
  const Backbone backb(handler.index_name, gcsa);
  if(!gcsa.isOk()) { return 2; }
  if(!backb.isOk()) { return 3; }
  if(handler.patterns_name == 0) { return 0; }
  BWASearch<GCSA> bwasearch(gcsa);

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
    std::vector<MatchInfo*>* results = bwasearch.find(rows[i], handler.k, handler.indels, handler.penalties);
    //std::cerr << "pattern = " << rows[i] << std::endl;
    usint temp = outputOccurrences(i, gcsa, backb, rows[i], results, handler.max_matches);
    if(temp > 0)
    {
        forward_occurrences += temp; match_counts[i] += temp;
        found++;
        edit_distances[(*(results->begin()))->errors]++;
    }
    bwasearch.deleteResults(results);
  }
  double time = readTimer() - start;
  double megabases = total / (double)MILLION;

  std::cerr << "Patterns:     " << n << " (" << (n / time) << " / sec)" << std::endl;

  std::cerr << "Total size:   " << megabases << " megabases";
  if(!(handler.locate))
  {
    std::cerr << " (" << (megabases / time) << " megabases / sec)";
  }
  std::cerr << std::endl;

  std::cerr << "Found:        " << found;
  if(handler.k == 0 && handler.reverse_complement)
  {
    std::cerr << " (" << forward << " forward, " << reverse << " reverse complement)";
  }
  std::cerr << std::endl;

  if(handler.locate)
  {
    std::cerr << "Occurrences:  " << (forward_occurrences + reverse_occurrences) << " (";
    if(handler.k == 0 && handler.reverse_complement)
    {
      std::cerr << forward_occurrences << " forward, " << reverse_occurrences << " reverse complement, ";
    }
    std::cerr << ((forward_occurrences + reverse_occurrences) / time) << " / sec)" << std::endl;
  }

  std::cerr << "Time:         " << time << " seconds" << std::endl;
  std::cerr << std::endl;

  if(handler.verbose)
  {
    std::cerr << "MATCH STATISTICS" << std::endl << std::endl;
    std::cerr << "Edit distances for matching patterns" << std::endl;
    for(usint i = 0; i <= handler.k; i++) { std::cerr << "  " << i << "  " << edit_distances[i] << std::endl; }
    std::cerr << std::endl;

    std::sort(match_counts, match_counts + n);
    std::cerr << "Number of matches for matching patterns" << std::endl;
    usint curr = 0, count = 0;
    for(usint i = 0; i < n; i++)
    {
      if(match_counts[i] != curr)
      {
        if(curr > 0) { std::cerr << "  " << curr << "  " << count << std::endl; }
        curr = match_counts[i]; count = 1;
      }
      else
      {
        count++;
      }
    }
    if(curr > 0) { std::cerr << "  " << curr << "  " << count << std::endl; }
    std::cerr << std::endl;
  }
  delete edit_distances; edit_distances = 0;
  delete match_counts; match_counts = 0;

  return 0;
}
