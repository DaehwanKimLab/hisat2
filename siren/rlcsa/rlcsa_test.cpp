#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "rlcsa.h"
#include "adaptive_samples.h"
#include "suffixarray.h"
#include "docarray.h"

#ifdef MULTITHREAD_SUPPORT
#include <omp.h>
#endif


using namespace CSA;


const int MAX_THREADS = 64;


struct Pattern
{
  Pattern(const std::string& _pattern, bool ignore_tab, usint ignore) :
    pattern(_pattern), weight(1), range(EMPTY_PAIR), steps(0), docc(0)
  {
    if(ignore_tab)
    {
      usint pos = this->pattern.find('\t');
      this->weight = atoi(this->pattern.substr(0, pos).c_str());
      this->pattern = this->pattern.substr(pos + 1);
    }
    else if(ignore > 0)
    {
      this->weight = atoi(this->pattern.substr(0, ignore).c_str());
      this->pattern = this->pattern.substr(ignore);
    }
  }

  inline bool found() const { return !isEmpty(this->range); }
  inline usint occ() const { return length(this->range); }
  inline double ratio() const { return this->occ() / (double)(this->docc); }

  std::string pattern;
  weight_type weight;
  pair_type range;
  usint steps;
  usint docc;
};

struct PatternOccDoccComparator
{
  bool operator() (const Pattern& a, const Pattern& b)
  {
    if(a.found() && b.found()) { return (a.ratio() > b.ratio()); }
    return a.found();
  }
} pod_comparator;



void printUsage();
void printPatterns(std::vector<Pattern>& patterns, const std::string& file_name, bool only_found = false);
void createPatterns(std::vector<std::string>& rows, std::vector<Pattern>& patterns, bool ignore_tab, usint ignore);

usint totalFound(const std::vector<Pattern>& patterns);
usint totalSize(const std::vector<Pattern>& patterns);
usint totalOcc(const std::vector<Pattern>& patterns);
usint totalSteps(const std::vector<Pattern>& patterns);
usint totalDocc(const std::vector<Pattern>& patterns);


int main(int argc, char** argv)
{
  std::cout << "RLCSA test" << std::endl;
  std::cout << std::endl;

  bool adaptive = false, direct = false, locate = false, pizza = false, count_steps = false;
  bool use_sa = false;
  bool listing = false, rle = false;
  usint ignore = 0, generate = 0;
  bool ignore_tab = false;
  bool write = false, write_patterns = false, sort_patterns = false;
  char* base_name = 0;
  char* patterns_name = 0;
  #ifdef MULTITHREAD_SUPPORT
  sint threads = 1;
  #endif
  for(int i = 1; i < argc; i++)
  {
    if(argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
        case 'a':
          adaptive = true; break;
        case 'd':
          direct = true; break;
        case 'g':
          generate = atoi(argv[i] + 2); break;
        case 'i':
          if(argv[i][2] == 't') { ignore = 1; ignore_tab = true; }
          else                  { ignore = atoi(argv[i] + 2); }
          break;
        case 'l':
          locate = true; break;
        case 'L':
          listing = true; break;
        case 'o':
          sort_patterns = true; break;
        case 'p':
          pizza = true; break;
        case 'r':
          rle = true; break;
        case 'S':
          use_sa = true; break;
        case 's':
          count_steps = true; break;
        case 'W':
          write_patterns = true; break;
        case 'w':
          write = true; break;
        default:
          std::cout << "Invalid option: " << argv[i] << std::endl << std::endl;
          printUsage();
          return 1;
      }
    }
    else if(base_name == 0)
    {
      base_name = argv[i];
    }
    else if(patterns_name == 0)
    {
      patterns_name = argv[i];
    }
    #ifdef MULTITHREAD_SUPPORT
    else
    {
      threads = std::min(MAX_THREADS, std::max(atoi(argv[i]), 1));
    }
    #endif
  }
  if(base_name == 0) { printUsage(); return 2; }

  std::cout << "Options:";
  if(generate > 0) { std::cout << " generate=" << generate; }
  if(ignore_tab) { std::cout << " ignore=tab"; }
  else if(ignore > 0) { std::cout << " ignore=" << ignore; }
  if(use_sa) { adaptive = direct = count_steps = listing = false; }
  if(listing)
  {
    locate = adaptive = count_steps = write = false;
    std::cout << " listing";
    if(direct || rle)
    {
      std::cout << "(";
      if(direct)      { std::cout << " direct"; }
      if(rle)         { std::cout << " rle"; }
      if(count_steps) { std::cout << " steps"; write = false; }
      std::cout << " )";
    }
  }
  else if(locate)
  {
    rle = sort_patterns = false;
    std::cout << " locate";
    if(adaptive || direct || count_steps)
    {
      std::cout << "(";
      if(adaptive)    { std::cout << " adaptive"; direct = true; }
      if(direct)      { std::cout << " direct"; }
      if(count_steps) { std::cout << " steps"; write = false; }
      std::cout << " )";
    }
  }
  else { adaptive = direct = count_steps = write = rle = sort_patterns = false; }
  if(sort_patterns) { std::cout << " sort_patterns"; }
  if(pizza) { std::cout << " pizza"; }
  if(use_sa) { std::cout << " sa"; }
  if(write_patterns) { std::cout << " write_patterns"; }
  if(write) { std::cout << " write"; }
  std::cout << std::endl;
  std::cout << "Base name: " << base_name << std::endl;
  if(patterns_name != 0) { std::cout << "Patterns: " << patterns_name << std::endl; }
  #ifdef MULTITHREAD_SUPPORT
  std::cout << "Threads: " << threads << std::endl; 
  #endif
  std::cout << std::endl;


  const RLCSA* rlcsa = (use_sa ? 0 : new RLCSA(base_name, false));
  const SuffixArray* sa = (use_sa ? new SuffixArray(base_name, false) : 0);
  usint size = 0, text_size = 0;
  if(use_sa)
  {
    if(!sa->isOk())
    {
      delete rlcsa; rlcsa = 0;
      delete sa; sa = 0;
      return 3;
    }
    sa->reportSize(true);
    size = sa->getSize();
    text_size = size;
  }
  else
  {
    if(!rlcsa->isOk())
    {
      delete rlcsa; rlcsa = 0;
      delete sa; sa = 0;
      return 3;
    }
    rlcsa->printInfo();
    rlcsa->reportSize(true);
    size = rlcsa->getSize();
    text_size = rlcsa->getTextSize();
  }
  if(patterns_name == 0)
  {
    delete rlcsa; rlcsa = 0;
    delete sa; sa = 0;
    return 0;
  }

  AdaptiveSamples* adaptive_samples = 0;
  if(adaptive)
  {
    adaptive_samples = new AdaptiveSamples(*rlcsa, base_name);
    adaptive_samples->report();
    if(!adaptive_samples->isOk())
    {
      delete rlcsa; rlcsa = 0;
      delete sa; sa = 0;
      delete adaptive_samples; adaptive_samples = 0;
      return 4;
    }
  }

  DocArray* docarray = 0;
  if(listing)
  {
    if(direct)
    {
      docarray = new DocArray(*rlcsa);
    }
    else
    {
      docarray = new DocArray(*rlcsa, base_name);
      docarray->reportSize(true);
      if(!docarray->isOk() || !docarray->hasGrammar())
      {
        delete rlcsa; rlcsa = 0;
        delete sa; sa = 0;
        delete adaptive_samples; adaptive_samples = 0;
        delete docarray; docarray = 0;
        return 5;
      }
    }
  }

  std::ifstream pattern_file(patterns_name, std::ios_base::binary);
  if(!pattern_file)
  {
    std::cerr << "Error opening pattern file!" << std::endl;
    delete rlcsa; rlcsa = 0;
    delete sa; sa = 0;
    delete adaptive_samples; adaptive_samples = 0;
    return 6;
  }
  std::vector<std::string> rows;
  if(!pizza) { readRows(pattern_file, rows, true); }
  else       { readPizzaChili(pattern_file, rows); }
  pattern_file.close();
  if(rows.size() == 0)
  {
    std::cerr << "No patterns found!" << std::endl;
    delete rlcsa; rlcsa = 0;
    delete sa; sa = 0;
    delete adaptive_samples; adaptive_samples = 0;
    return 7;
  }
  std::vector<Pattern> patterns;
  createPatterns(rows, patterns, ignore_tab, ignore);
  rows.clear();


  // Generate random patterns according to weights.
  if(generate > 0)
  {
    DeltaEncoder encoder(2 * sizeof(usint));
    usint sum = 0;
    for(std::vector<Pattern>::iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
    {
      encoder.addBit(sum);
      sum += iter->weight;
    }
    encoder.flush();
    DeltaVector vector(encoder, sum);

    std::vector<Pattern> new_patterns; new_patterns.reserve(generate);
    DeltaVector::Iterator viter(vector);
    srand(0xDEADBEEF);  // FIXME random seed
    for(usint i = 0; i < generate; i++)
    {
      usint tmp = viter.rank(rand() % sum);
      new_patterns.push_back(patterns[tmp - 1]);
    }
    patterns.clear();
    patterns.assign(new_patterns.begin(), new_patterns.end());
    new_patterns.clear();
    std::cout << "Generated " << patterns.size() << " patterns." << std::endl;
    std::cout << std::endl;
  }

  // Prepare weights for writing the number of occurrences.
  weight_type* totals = 0;
  if(write)
  {
    totals = new weight_type[text_size];
    for(usint i = 0; i < text_size; i++) { totals[i] = 0; }
  }


  // Actual work.
  #ifdef MULTITHREAD_SUPPORT
  omp_set_num_threads(threads);
  #endif
  double start = readTimer();

  #pragma omp parallel for schedule(dynamic, 1)
  for(usint i = 0; i < patterns.size(); i++)
  {
    patterns[i].range = (use_sa ? sa->count(patterns[i].pattern) : rlcsa->count(patterns[i].pattern));

    if(locate && patterns[i].found())
    {
      usint* matches = 0;
      uint* sa_matches = 0;
      if(adaptive)    { matches = adaptive_samples->locate(patterns[i].range, count_steps); }
      else if(use_sa) { sa_matches = sa->locate(patterns[i].range); }
      else            { matches = rlcsa->locate(patterns[i].range, direct, count_steps); }
      if(count_steps)
      {
        patterns[i].steps = 0;
        for(usint j = 0; j < patterns[i].occ(); j++)
        {
          patterns[i].steps += matches[j];
        }
      }

      if(write)
      {
        #pragma omp critical
        {
          if(use_sa) { for(usint j = 0; j < patterns[i].occ(); j++) { totals[sa_matches[j]] += patterns[i].weight; } }
          else       { for(usint j = 0; j < patterns[i].occ(); j++) { totals[matches[j]] += patterns[i].weight; } }
        }
      }
      delete[] sa_matches; sa_matches = 0;
      delete[] matches; matches = 0;
    }
    if(listing && patterns[i].found())
    {
      if(rle)
      {
        std::vector<pair_type>* docs =
          (direct ? docarray->directListingRLE(patterns[i].range)
                  : docarray->listDocumentsRLE(patterns[i].range));
        if(docs != 0)
        {
          patterns[i].docc = 0;
          for(std::vector<pair_type>::iterator iter = docs->begin(); iter != docs->end(); ++iter)
          {
            patterns[i].docc += length(*iter);
          }
        }
        delete docs;
      }
      else
      {
        std::vector<usint>* docs =
          (direct ? docarray->directListing(patterns[i].range)
                  : docarray->listDocuments(patterns[i].range));
        if(docs != 0)
        {
          patterns[i].docc = docs->size();
          delete docs; docs = 0;
        }
      }
    }
  }

  double seconds = readTimer() - start;
  double megabytes = totalSize(patterns) / (double)MEGABYTE;


  // Report.
  std::cout << "Patterns:      " << patterns.size() << " (" << (patterns.size() / seconds) << " / sec)" << std::endl;
  std::cout << "Total size:    " << megabytes << " MB";
  if(!locate)
  {
    std::cout << " (" << (megabytes / seconds) << " MB/s)";
  }
  std::cout << std::endl;

  usint total_occ = totalOcc(patterns), total_steps = totalSteps(patterns);
  std::cout << "Found:         " << totalFound(patterns) << std::endl;
  std::cout << "Occurrences:   " << total_occ;
  if(locate)
  {
    std::cout << " (" << (total_occ / seconds) << " / sec)";
    if(count_steps)
    {
      std::cout << std::endl
                << "Locate steps:  " << total_steps << " (" << (total_steps / (double)total_occ) << " / occurrence)";
    }
  }
  std::cout << std::endl;
  if(listing)
  {
    usint total_docc = totalDocc(patterns);
    std::cout << "Documents:     " << total_docc << " (" << (total_docc / seconds) << " / sec)" << std::endl;
  }
  std::cout << "Time:          " << seconds << " seconds" << std::endl;
  std::cout << std::endl;


  // Write.
  if(write)
  {
    std::string output_name(patterns_name); output_name += ".found";
    std::ofstream output(output_name.c_str(), std::ios_base::binary);
    if(!output)
    {
      std::cerr << "Cannot open output file (" << output_name << ")!" << std::endl;
    }
    else
    {
      output.write((char*)totals, text_size * sizeof(weight_type));
      output.close();
    }

    weight_type max_weight = 0, total_weight = 0;
    for(usint i = 0; i < size; i++) 
    {
      max_weight = std::max(max_weight, totals[i]);
      total_weight += totals[i];
    }
    std::cout << "Weights:       " << max_weight << " (max), " 
              << (((double)total_weight) / size) << " (average)" << std::endl;
    std::cout << std::endl;
  }

  if(write_patterns)
  {
    printPatterns(patterns, std::string(patterns_name) + ".selected");
    printPatterns(patterns, std::string(patterns_name) + ".found", true);
  }
  if(sort_patterns)
  {
    sequentialSort(patterns.begin(), patterns.end(), pod_comparator);
    printPatterns(patterns, std::string(patterns_name) + ".sorted", true);
  }

  // Cleanup.

  if(adaptive) { adaptive_samples->report(); }
  delete rlcsa; rlcsa = 0;
  delete sa; sa = 0;
  delete adaptive_samples; adaptive_samples = 0;
  delete docarray; docarray = 0;
  delete[] totals; totals = 0;
  return 0;
}


void
printUsage()
{
  std::cout << "Usage: rlcsa_test [options] base_name [patterns [threads]]" << std::endl;
  std::cout << "  -a   Use adaptive samples." << std::endl;
  std::cout << "  -d   Use direct locate / document listing." << std::endl;
  std::cout << "  -g#  Use the weights to generate # actual patterns." << std::endl;
  std::cout << "  -i#  Ignore first # characters of each pattern." << std::endl;
  std::cout << "  -l   Locate the occurrences." << std::endl;
  std::cout << "  -L   List the documents containing the pattern." << std::endl;
  std::cout << "  -o   Write the patterns sorted by occ/docc into patterns.sorted." << std::endl;
  std::cout << "  -p   Pattern file is in Pizza & Chili format." << std::endl;
  std::cout << "  -r   Run-length encode the results (requires -L)." << std::endl;
  std::cout << "  -S   Use a plain suffix array (negates adaptive, direct, steps)." << std::endl;
  std::cout << "  -s   Count the number of steps required for locate() (negates write)." << std::endl;
  std::cout << "  -W   Write selected and found patterns into patterns.selected and patterns.found." << std::endl;
  std::cout << "  -w   Write the weighted number of occurrences of each suffix into patterns.found." << std::endl;
  std::cout << std::endl;
  std::cout << "Default pattern weight for -g and -w is 1." << std::endl;
  std::cout << "If -i is specified, the ignored characters are used as weights." << std::endl;
  std::cout << "Use -it to ignore the pattern until the first \\t." << std::endl;
  std::cout << std::endl;
}

void
printPatterns(std::vector<Pattern>& patterns, const std::string& file_name, bool only_found)
{
  if(patterns.empty()) { return; }

  std::ofstream output(file_name.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error: Cannot open file " << file_name << " for writing patterns!" << std::endl;
    return;
  }

  for(std::vector<Pattern>::iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
  {
    if(!only_found || iter->found()) { output << iter->pattern << std::endl; }
  }
  output.close();
}

void
createPatterns(std::vector<std::string>& rows, std::vector<Pattern>& patterns, bool ignore_tab, usint ignore)
{
  for(std::vector<std::string>::iterator iter = rows.begin(); iter != rows.end(); ++iter)
  {
    patterns.push_back(Pattern(*iter, ignore_tab, ignore));
  }
}

usint
totalFound(const std::vector<Pattern>& patterns)
{
  double total = 0;
  for(std::vector<Pattern>::const_iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
  {
    if(iter->found()) { total++; }
  }
  return total;
}

usint
totalSize(const std::vector<Pattern>& patterns)
{
  double total = 0;
  for(std::vector<Pattern>::const_iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
  {
    total += iter->pattern.length();
  }
  return total;
}

usint
totalOcc(const std::vector<Pattern>& patterns)
{
  double total = 0;
  for(std::vector<Pattern>::const_iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
  {
    total += iter->occ();
  }
  return total;
}

usint
totalSteps(const std::vector<Pattern>& patterns)
{
  double total = 0;
  for(std::vector<Pattern>::const_iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
  {
    total += iter->steps;
  }
  return total;
}

usint
totalDocc(const std::vector<Pattern>& patterns)
{
  double total = 0;
  for(std::vector<Pattern>::const_iterator iter = patterns.begin(); iter != patterns.end(); ++iter)
  {
    total += iter->docc;
  }
  return total;
}
