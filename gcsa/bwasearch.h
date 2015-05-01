#ifndef BWASEARCH_H
#define BWASEARCH_H


#include <algorithm>
#include <vector>
#include <assert.h>
#include <misc/utils.h>


namespace GCSA
{
	typedef CSA::pair_type pair_type;
	typedef CSA::uchar uchar;
	typedef CSA::sint sint;
	typedef CSA::usint usint;
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

    pair_type find(const std::vector<usint>& pattern, bool reverse_complement, usint skip = 0) const
        {
            std::cout << "find(pattern, reverse_complement=" << reverse_complement << ", skip=" << skip << ")" << std::endl;
            std::vector<usint> pat, revpat;
            for (usint i = skip; i < pattern.size(); ++i) {
                pat.push_back(pattern[i]);
            }//= pattern.substr(skip);
            revpat = pat; //FIXME: do something more effient than copy the whole vector, change code to iterate in rev order
            std::reverse(revpat.begin(), revpat.end());  
      
            if(pat.size() == 0) { return this->index.getSARange(); }





            unsigned int myc = (reverse_complement ? this->complement(pat[0]) : pat[pat.size() - 1]);

            usint delta = 450;
            pair_type myinitrange = this->index.getSARange();
            std::vector<long unsigned int> hits = this->index.restricted_unique_range_values(myinitrange.first, myinitrange.second, 
                                                                                             myc - delta, myc + delta);

            // actual algo
            for(std::vector<long unsigned int>::iterator itr = hits.begin(); itr != hits.end(); ++itr) {
                



                pair_type myrange = /*this->index.getSARange()*/ this->index.getCharRange(*itr);
                //pair_type myrange = this->index.getSARange();// this->index.getCharRange(myc);
                myrange.second += 1;
                std::cout << "bootstrap range for initial query symbol candidate " <<*itr<< " is [" << myrange.first << ".." << myrange.second << "]" << std::endl;
                //pair_type retrange = 
                this->mybackwardSearch(pat, pat.size() - 1 , myrange);
                //if (!CSA::isEmpty(retrange)) return retrange; //fixme
                // //FIXME: reverse pattern here and rerun
                //std::cout << "DEBUG: my reverse search for pattern:" << std::endl;
                //this->mybackwardSearch(revpat, revpat.size() - 1 , myrange);
                
            }
            // std::cout << "DEBUG: normal search for pattern:" << std::endl;
            // unsigned int c = (reverse_complement ? this->complement(pat[0]) : pat[pat.size() - 1]);
            // pair_type range = /*this->index.getSARange()*/ this->index.getCharRange(c);
            // if(CSA::isEmpty(range)) { return range; }

            // MatchInfo info(1, range, (reverse_complement ? MatchInfo::REVERSE_COMPLEMENT : 0));
            // std::cout << "Normal search initial interval is [" << info.range.first <<".."<<info.range.second<<"]" << std::endl;
            // this->backwardSearch(pat, info);
            // std::cout << "Normal search final interval is [" << info.range.first <<".."<<info.range.second<<"]" << std::endl;

            // this->index.convertToSARange(info.range);  
            //return info.range;
            return pair_type(1,0);
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
    std::vector<MatchInfo*>* find(const std::vector<usint>& pattern, usint k, bool allow_indels, bool weighted) const
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
        if(info->matched >= pattern.size())
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
          c = pattern[pattern.size() - info->matched - 1];
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
          info = (CSA::isEmpty(rng) ? 0 : info->match(rng, c));
          continue;
        }

        // Second visit. Continue by expanding the edit operations.
        if(bounds[info->matched + 1] < k - info->errors)
        {
          for(usint i = 0; i < ALPHABET_SIZE; i++)
          {
            if(alphabet->at(i) == c) { continue; }
            pair_type rng = this->index.LF(info->range, alphabet->at(i));
            if(CSA::isEmpty(rng)) { continue; }
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
      if(CSA::isEmpty(range)) { return 0; }
      usint temp = CSA::length(range);

      if(locate)
      {
        std::vector<usint>* occurrences = this->index.locateRange(range);
        if(occurrences != 0)
        {
            
            std::cout << "Found " << occurrences->size() << " match(s) located at: " ;
            for (std::vector<usint>::iterator mi = occurrences->begin(); mi != occurrences->end(); ++mi) {
                std::cout << *mi << ", ";
            }
            std::cout << std::endl;
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
          temp += CSA::length(*iter);
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
          temp += CSA::length((*iter)->range);
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

    //TODO: convert this to recursive call
    void mybackwardSearch(const std::vector<usint>& pattern,  unsigned int it, pair_type range) const
        {

            if (it == 0) { // match complete
                for(int i=0; i < pattern.size() - it; ++i) std::cout << "\t";
                //std::cout << "Found match in interval [" << range.first << ".." << range.second << "]" << std::endl;
                std::vector<usint>* occurrences = this->index.locateRange(range);
                if(occurrences != 0)
                {
            
                    std::cout << "Found " << occurrences->size() << " match(s) located at: " ;
                    for (std::vector<usint>::iterator mi = occurrences->begin(); mi != occurrences->end(); ++mi) {
                        std::cout << *mi << ", ";
                    }
                    std::cout << std::endl;
                    
                    delete occurrences;
                }

                //return range; //FIXME - returns only first match this way
            } else {
                for(int i=0; i < pattern.size() - it; ++i) std::cout << "\t";
                std::cout << "mybackwardSsearch(pattern[" << it -1 << "] /* "<< pattern[it-1] << " */, range=<" <<range.first << "," << range.second << ">)" <<  std::endl;
                // int lookahead = 1;
                // //trim lookahead to max remaining
                // if ((int)it - 1 - lookahead < 0) {
                //     lookahead = it - 1;
                // }
                // for (int actv_la = 0; actv_la <= lookahead; ++actv_la) {
                //     for(int i=0; i < pattern.size() - it; ++i) std::cout << "\t";
                //     std::cout << "active lookahead: " << actv_la << std::endl;

                //     // compute the sum of the next lookahead fragments
                //     unsigned int c = 0;
                //     for (int j = 0; j <= actv_la; ++j) {
                //         int index = it - 1 - j;
                //         assert(index >= 0);
                //         c += pattern[index];
                //     }
                int actv_la = 0;
                int c = pattern[it - 1];
                    //wt stuff
                    usint delta = 450;
                    std::vector<long unsigned int> hits = this->index.restricted_unique_range_values(range.first, range.second, 
                                                                                                     c - delta, c + delta);
                    for(int i=0; i < pattern.size() - it; ++i) std::cout << "\t";
                    // actual algo
                    for(std::vector<long unsigned int>::iterator itr = hits.begin(); itr != hits.end(); ++itr) {
                        pair_type new_range = this->index.LF(range, *itr); //FIXME: renenable WT later
                        //pair_type new_range = this->index.LF(range, c);
                        //for(int i=0; i < pattern.size() - it; ++i) std::cout << "\t";
//                        std::cout << "LF(<" <<range.first << "," << range.second << ">, " << *itr <<") = <" << new_range.first <<  "," << new_range.second << ">" << std::endl;
                        if(!CSA::isEmpty(new_range)) {
                            //pair_type retrange = 
                            this->mybackwardSearch(pattern, it - 1 - actv_la, new_range);
                            // if(!CSA::isEmpty(retrange)) {

                            //     //return retrange; //FIXME - returns only first match this way
                            // }
                        }
                    }
              //}
            }
            //return pair_type(1,0);


        }

    void backwardSearch(const std::vector<usint>& pattern, MatchInfo& info) const
    {
        std::cout << "matched is " << info.matched << std::endl;
      sint pos = (info.isReverseComplement() ? info.matched : pattern.size() - info.matched - 1);
      sint adj = (info.isReverseComplement() ? 1 : -1);

      for(; info.matched < pattern.size(); info.matched++, pos += adj)
      {
        unsigned int c = (info.isReverseComplement() ? this->complement(pattern[pos]) : pattern[pos]);

        //wavelet tree stuff
        // usint delta = 13;
        // std::vector<long unsigned int> hits = this->index.restricted_unique_range_values(info.range.first, info.range.second, 
        //                                                                                  c - delta, c + delta);
        // std::cout << "DEBUG - wavelet tree query in SA interval [" << info.range.first << ".."<< info.range.second 
        //           << "] has the following symbols within " << delta << " alphabet symbols of " << c << ": " ;
        // for(std::vector<long unsigned int>::iterator itr = hits.begin(); itr != hits.end(); ++itr) {
        //     std::cout << *itr << " ";
        // }
        // std::cout << std::endl;

        //std::cout << "LF(<" <<info.range.first << "," << info.range.second << ">, " << c <<") = <" ;
        info.range = this->index.LF(info.range, c);
        //std::cout << info.range.first <<  "," << info.range.second << ">" << std::endl;
        if(CSA::isEmpty(info.range)) { return; }
      }
      std::cout << "Found standard backwardSearch match!" <<std::endl;
    }

    usint* errorBounds(const std::vector<usint>& pattern, usint k, bool complement) const
    {
      // bounds[i]: lower bound for errors after matching 'i' characters.
      usint* bounds = new usint[pattern.size() + 1];

      usint errors = 0;
      usint start = 0, limit = pattern.size() - 1;
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
            usint temp = (complement ? pattern.size() - mid - 1 : start);

            // pattern.substr(temp, mid + 1 - start)
            std::vector<usint> pat;
            for (usint i = temp; i < mid + 1 - start; ++i) {
                pat.push_back(pattern[i]);
            }

            this->backwardSearch(pat, info);
            if(CSA::isEmpty(info.range))
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

          // pattern.substr(temp, pattern.size() - start)
          std::vector<usint> pat;
          for (usint i = temp; i < pattern.size() - start; ++i) {
              pat.push_back(pattern[i]);
          }
          
          this->backwardSearch(pat, info);
          if(CSA::isEmpty(info.range)) { errors++; }
          for(usint i = start; i <= limit; i++) { bounds[i] = errors; }
          break;
        }
      }

      std::reverse(bounds, bounds + pattern.size());
      bounds[pattern.size()] = 0;
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
