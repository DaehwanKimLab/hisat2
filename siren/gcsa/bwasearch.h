#ifndef BWASEARCH_H
#define BWASEARCH_H


#include <algorithm>
#include <vector>

#include <misc/utils.h>


namespace CSA
{

/*
  This uses the RANGES part of external module interface.
*/

//--------------------------------------------------------------------------

struct MatchInfo
{
  const static uint TYPE_MASK = 0xF;
  const static uint MATCH     = 0x1;
  const static uint MISMATCH  = 0x2;
  const static uint INSERTION = 0x4; // Insertion in pattern.
  const static uint DELETION  = 0x8; // Deletion in pattern.

  const static uint INHERITED_MASK     = 0x30;  // Flags inherited from parent.
  const static uint REVERSE_COMPLEMENT = 0x10;
  const static uint IS_WEIGHTED        = 0x20;

  const static uint HAS_BEEN_VISITED   = 0x40;
  const static uint HAS_BEEN_EXPANDED  = 0x80;

  const static uint LOWEST_PENALTY   = 3;   // Lowest penalty for an error.
  const static uint MISMATCH_PENALTY = 3;
  const static uint GAP_OPEN         = 11;
  const static uint GAP_EXTEND       = 4;

  uint     score;
  uint     errors;
  uint     flags;
  uint     matched;

  pair_type range;  // This is a BWT range.
  uint     c;

  uint       children;
  MatchInfo* parent;

  // We have matched 0 characters.
  MatchInfo(pair_type rng, uint _flags) :
    score(0), errors(0), flags(MATCH | _flags),
    matched(0),
    range(rng), c(0),
    children(0), parent(0)
  {
  }

  // We have matched _matched characters.
  MatchInfo(uint _matched, pair_type rng, uint _flags) :
    score(0), errors(0), flags(MATCH | _flags),
    matched(_matched),
    range(rng), c(0),
    children(0), parent(0)
  {
  }

  MatchInfo(MatchInfo* par, pair_type rng, uint _c, uint type) :
    flags(type | (INHERITED_MASK & par->flags)),
    range(rng), c(_c),
    children(0), parent(par)
  {
    this->parent->children++;

    switch(type)
    {
      case MATCH:
        this->score = this->parent->score;
        this->errors = this->parent->errors;
        this->matched = this->parent->matched + 1;
        break;
      case MISMATCH:
        if(this->isWeighted()) { this->score = this->parent->score + MISMATCH_PENALTY; }
        else                   { this->score = this->parent->score + 1; }
        this->errors = this->parent->errors + 1;
        this->matched = this->parent->matched + 1;
        break;
      case INSERTION:
        if(this->isWeighted())
        {
          this->score = this->parent->score + ((this->parent->isInsertion()) ? GAP_EXTEND : GAP_OPEN);
        }
        else { this->score = this->parent->score + 1; }
        this->errors = this->parent->errors + 1;
        this->matched = this->parent->matched + 1;
        break;
      case DELETION:
        if(this->isWeighted())
        {
          this->score = this->parent->score + ((this->parent->isDeletion()) ? GAP_EXTEND : GAP_OPEN);
        }
        else { this->score = this->parent->score + 1; }
        this->errors = this->parent->errors + 1;
        this->matched = this->parent->matched;
        break;
    }
  }

  ~MatchInfo()
  {
    if(this->parent != 0) { this->parent->children--; }
  }

  MatchInfo* match(pair_type rng, usint next)
  {
    return new MatchInfo(this, rng, next, MATCH);
  }

  MatchInfo* mismatch(pair_type rng, usint next)
  {
    return new MatchInfo(this, rng, next, MISMATCH);
  }

  MatchInfo* insertion(usint next)
  {
    return new MatchInfo(this, this->range, next, INSERTION);
  }

  MatchInfo* deletion(pair_type rng, usint next)
  {
    return new MatchInfo(this, rng, next, DELETION);
  }

  bool isMatch()     { return this->flags & MATCH; }
  bool isMismatch()  { return this->flags & MISMATCH; }
  bool isInsertion() { return this->flags & INSERTION; }
  bool isDeletion()  { return this->flags & DELETION; }

  bool hasFlags(uint _flags) { return this->flags & _flags; }

  bool isReverseComplement() { return this->flags & REVERSE_COMPLEMENT; }
  bool isWeighted()          { return this->flags & IS_WEIGHTED; }
  bool hasBeenVisited()      { return this->flags & HAS_BEEN_VISITED; }
  bool hasBeenExpanded()     { return this->flags & HAS_BEEN_EXPANDED; }

  void visit()               { this->flags |= HAS_BEEN_VISITED; }
  void expand()              { this->flags |= HAS_BEEN_EXPANDED; }

  uint getEstimate()
  {
    if(!(this->hasBeenVisited())) { return this->score; }
    return this->score + (this->isWeighted() ? LOWEST_PENALTY : 1);
  }

  private:
    // These are not allowed.
    MatchInfo();
    MatchInfo(const MatchInfo&);
    MatchInfo& operator = (const MatchInfo&);
};


// STL has a maximum heap, so operator() should return true if a is worse than b.
struct MatchInfoComparator
{
  bool operator() (MatchInfo* a, MatchInfo* b) const
  {
    uint a_key = a->getEstimate();
    uint b_key = b->getEstimate();
    if(a_key != b_key) { return (a_key > b_key); }

    return (a->matched < b->matched);
//    return (a->score > b->score);
  }
};

//--------------------------------------------------------------------------

template<class Index>
class BWASearch
{
  public:

    const static usint INITIAL_STEP = 8;  // Should be >= log(n) / 4.

    BWASearch(const Index& _index) :
      index(_index), ALPHABET("ACGTN"), COMPLEMENT("TGCAN")
    {
    }

    pair_type find(const std::string& pattern, bool reverse_complement, usint skip = 0) const
    {
      std::string pat = pattern.substr(skip);
      if(pat.length() == 0) { return this->index.getSARange(); }

      uchar c = (reverse_complement ? this->complement(pat[0]) : pat[pat.length() - 1]);
      pair_type range = this->index.getCharRange(c);
      if(isEmpty(range)) { return range; }

      MatchInfo info(1, range, (reverse_complement ? MatchInfo::REVERSE_COMPLEMENT : 0));
      this->backwardSearch(pat, info);
      this->index.convertToSARange(info.range);

      return info.range;
    }

    /*
      Approximate search with at most k mismatches/errors for the pattern and its
      reverse complement. Returns only the best matches.
      The returned vector must be deleted with deleteResults().
      We only match characters to "ACGTN" in this search.

      For the best results, try exact matching before approximate matching.

      The algorithm is similar to the one in:

      H. Li, R. Durbin: Fast and accurate short read alignment with
      Burrows-Wheeler transform. Bioinformatics, 2009.
    */
    std::vector<MatchInfo*>* find(const std::string& pattern, usint k, bool allow_indels, bool weighted) const
    {
      std::vector<MatchInfo*>* results = new std::vector<MatchInfo*>;
      std::vector<MatchInfo*>  heap;
      MatchInfoComparator      comp;

      pair_type range = this->index.getBWTRange();
      usint* fw_bounds = this->errorBounds(pattern, k, false);
      uint flags = (weighted ? MatchInfo::IS_WEIGHTED : 0);
      if(fw_bounds[0] <= k)
      {
        heap.push_back(new MatchInfo(range, flags));
        push_heap(heap.begin(), heap.end(), comp);
      }
      usint* rc_bounds = this->errorBounds(pattern, k, true);
      if(rc_bounds[0] <= k)
      {
        heap.push_back(new MatchInfo(range, flags | MatchInfo::REVERSE_COMPLEMENT));
        push_heap(heap.begin(), heap.end(), comp);
      }

      uint best_match = (uint)-1;
      MatchInfo* info = 0;
      char c = 0;
      usint* bounds = 0;
      const std::string* alphabet = 0;
      while(true)
      {
        if(info == 0)
        {
          if(heap.empty()) { break; }
          pop_heap(heap.begin(), heap.end(), comp);
          info = heap.back(); heap.pop_back();
          if(info->getEstimate() > best_match) { info->expand(); this->deleteBranch(info); break; }
        }
        if(info->matched >= pattern.length())
        {
          results->push_back(info);
          best_match = std::min(best_match, info->score);
          info = 0; continue;
        }

        if(info->isReverseComplement())
        {
          c = this->complement(pattern[info->matched]);
          bounds = rc_bounds;
          alphabet = &(this->COMPLEMENT);
        }
        else
        {
          c = pattern[pattern.length() - info->matched - 1];
          bounds = fw_bounds;
          alphabet = &(this->ALPHABET);
        }

        // First visit to this MatchInfo. Continue by expanding the matching path.
        if(!(info->hasBeenVisited()))
        {
          info->visit();
          heap.push_back(info);
          push_heap(heap.begin(), heap.end(), comp);
          pair_type rng = this->index.LF(info->range, c);
          info = (isEmpty(rng) ? 0 : info->match(rng, c));
          continue;
        }

        // Second visit. Continue by expanding the edit operations.
        if(bounds[info->matched + 1] < k - info->errors)
        {
          for(usint i = 0; i < ALPHABET_SIZE; i++)
          {
            if(alphabet->at(i) == c) { continue; }
            pair_type rng = this->index.LF(info->range, alphabet->at(i));
            if(isEmpty(rng)) { continue; }
            heap.push_back(info->mismatch(rng, alphabet->at(i)));
            push_heap(heap.begin(), heap.end(), comp);

            if(allow_indels && bounds[info->matched] < k - info->errors &&
               info->hasFlags(MatchInfo::MATCH | MatchInfo::DELETION))
            {
              heap.push_back(info->deletion(rng, alphabet->at(i)));
              push_heap(heap.begin(), heap.end(), comp);
            }
          }
          if(allow_indels && info->hasFlags(MatchInfo::MATCH | MatchInfo::INSERTION))
          {
            heap.push_back(info->insertion(c));
            push_heap(heap.begin(), heap.end(), comp);
          }
        }
        info->expand(); this->deleteBranch(info); info = 0;
      }

      // Delete remaining suboptimal branches of the match tree.
      // Note that forcing deletion might delete parent before child.
      for(std::vector<MatchInfo*>::iterator iter = heap.begin(); iter != heap.end(); ++iter)
      {
        (*iter)->expand(); this->deleteBranch(*iter);
      }

      delete[] fw_bounds; fw_bounds = 0;
      delete[] rc_bounds; rc_bounds = 0;
      return results;
    }

    usint handleOccurrences(pair_type range, bool locate, usint max_matches)
    {
      if(isEmpty(range)) { return 0; }
      usint temp = length(range);

      if(locate)
      {
        std::vector<usint>* occurrences = this->index.locateRange(range);
        if(occurrences != 0)
        {
          temp = occurrences->size();
          delete occurrences;
        }
      }

      if(max_matches > 0 && temp > max_matches) { return max_matches + 1; }
      return temp;
    }

    usint handleOccurrences(std::vector<pair_type>* ranges, bool locate, usint max_matches)
    {
      if(ranges->empty()) { return 0; }
      usint temp = 0;

      if(locate)
      {
        std::vector<usint>* occurrences = this->index.locateRanges(*ranges);
        if(occurrences != 0)
        {
          temp = occurrences->size();
          delete occurrences;
        }
      }
      else
      {
        for(std::vector<pair_type>::iterator iter = ranges->begin(); iter != ranges->end(); ++iter)
        {
          temp += length(*iter);
        }
      }

      if(max_matches > 0 && temp > max_matches) { return max_matches + 1; }
      return temp;
    }

    usint handleOccurrences(std::vector<MatchInfo*>* results, bool locate, usint max_matches)
    {
      if(results->empty()) { return 0; }

      if(locate)
      {
        std::vector<pair_type> ranges;
        ranges.reserve(results->size());
        for(std::vector<MatchInfo*>::iterator iter = results->begin(); iter != results->end(); ++iter)
        {
          ranges.push_back((*iter)->range);
        }
        this->index.convertToSARange(ranges);
        return this->handleOccurrences(&ranges, locate, max_matches);
      }
      else
      {
        usint temp = 0;
        for(std::vector<MatchInfo*>::iterator iter = results->begin(); iter != results->end(); ++iter)
        {
          temp += length((*iter)->range);
        }
        if(max_matches > 0 && temp > max_matches) { return max_matches + 1; }
        return temp;
      }
    }

    // Deletes the match tree in the vector as well as the vector itself.
    void deleteResults(std::vector<MatchInfo*>* results)
    {
      if(results == 0) { return; }

      for(std::vector<MatchInfo*>::iterator iter = results->begin(); iter != results->end(); ++iter)
      {
        this->deleteBranch(*iter, true);
      }
      delete results;
    }

  private:
    const Index& index;

    const std::string ALPHABET;
    const std::string COMPLEMENT;
    static const usint ALPHABET_SIZE = 5;

    char complement(char c) const
    {
      size_t temp = this->COMPLEMENT.find(c);
      if(temp != std::string::npos) { c = this->ALPHABET[temp]; }
      return c;
    }

    void deleteBranch(MatchInfo* info, bool force = false) const
    {
      while(info != 0 && info->children == 0 && (force || info->hasBeenExpanded()))
      {
        MatchInfo* temp = info;
        info = temp->parent;
        delete temp;
      }
    }

    void backwardSearch(const std::string& pattern, MatchInfo& info) const
    {
      sint pos = (info.isReverseComplement() ? info.matched : pattern.length() - info.matched - 1);
      sint adj = (info.isReverseComplement() ? 1 : -1);

      for(; info.matched < pattern.length(); info.matched++, pos += adj)
      {
        uchar c = (info.isReverseComplement() ? this->complement(pattern[pos]) : pattern[pos]);
        info.range = this->index.LF(info.range, c);
        if(isEmpty(info.range)) { return; }
      }
    }

    usint* errorBounds(const std::string& pattern, usint k, bool complement) const
    {
      // bounds[i]: lower bound for errors after matching 'i' characters.
      usint* bounds = new usint[pattern.length() + 1];

      usint errors = 0;
      usint start = 0, limit = pattern.length() - 1;
      while(start <= limit)
      {
        // Search for the lowest 'low' such that there is no match for 'pattern[start, low]'.
        usint low = start, high = limit + 1;
        usint diff = INITIAL_STEP;

        if(errors < k)
        {
          while(low < high && low <= limit)
          {
            usint mid = std::min(low + 2 * diff, low + (high - low) / 2);
            diff = mid - low;
            MatchInfo info(this->index.getBWTRange(), (complement ? MatchInfo::REVERSE_COMPLEMENT : 0));
            usint temp = (complement ? pattern.length() - mid - 1 : start);
            this->backwardSearch(pattern.substr(temp, mid + 1 - start), info);
            if(isEmpty(info.range))
            {
              high = mid;
            }
            else
            {
              low = mid + 1;
            }
          }

          for(usint i = start; i < low; i++) { bounds[i] = errors; }
          errors++; start = low + 1;
          if(low <= limit) { bounds[low] = errors; }
        }
        else  // We cannot afford any more errors.
        {
          MatchInfo info(this->index.getBWTRange(), (complement ? MatchInfo::REVERSE_COMPLEMENT : 0));
          usint temp = (complement ? 0 : start);
          this->backwardSearch(pattern.substr(temp, pattern.length() - start), info);
          if(isEmpty(info.range)) { errors++; }
          for(usint i = start; i <= limit; i++) { bounds[i] = errors; }
          break;
        }
      }

      std::reverse(bounds, bounds + pattern.length());
      bounds[pattern.length()] = 0;
      return bounds;
    }

    // These are not allowed.
    BWASearch();
    BWASearch(const BWASearch&);
    BWASearch& operator = (const BWASearch&);
};

//--------------------------------------------------------------------------

};  // namespace CSA


#endif  // BWASEARCH_H
