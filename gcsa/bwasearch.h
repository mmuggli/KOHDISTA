#ifndef BWASEARCH_H
#define BWASEARCH_H


#include <algorithm>
#include <vector>
#include <assert.h>
#include <misc/utils.h>
#include <unordered_set>
#include <set>
#include <utility>
#include <boost/math/distributions/chi_squared.hpp>
#include <math.h>
namespace GCSA
{
	typedef CSA::pair_type pair_type;
	typedef CSA::uchar uchar;
	typedef CSA::sint sint;
	typedef CSA::usint usint;
const float OM_STDDEV = 2.45884783995 * 1000.0; // based on ~23k valuev paired cutsite alignments 
//const float OM_STDDEV = 2.28463258304 * 1000.0; // based on ~16k valuev 1:1 frag alignment
//const float OM_STDDEV = 0.150 * 1000.0; // based on ~16k valuev 1:1 frag alignment
    //const uint DELTA = OM_STDDEV *  3.0;

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
    typedef std::pair<unsigned int, pair_type> work_t;
    BWASearch(const Index& _index,  const CSA::DeltaVector &_rmap_starts, const std::vector<std::pair<unsigned int, std::string> > _frag2rmap) :
      index(_index), rmap_starts(_rmap_starts), frag2rmap(_frag2rmap), ALPHABET("ACGTN"), COMPLEMENT("TGCAN")
    {
    }
    unsigned int get_stddev(const unsigned int frag_bp) const {
        double sigma_kbp = .58;
        double frag_kbp = (double)frag_bp / 1000.0;
        double variance_kbp = powf(sigma_kbp, 2);
        double expect_var_kbp = 2* variance_kbp * frag_kbp; // FIXME: do we need to double this? or does valuev .58 already account for noise in both frags
        double expect_stddev_kbp = sqrt(expect_var_kbp);
        double expect_stddev_bp = expect_stddev_kbp * 1000.0;
//        return OM_STDDEV;
//        return 100;
//        std::cout << "(stddev " << expect_stddev_bp << ")" << std::endl;
        return expect_stddev_bp;
    }

    pair_type find(const std::vector<usint>& pattern, bool reverse_complement, usint skip = 0) const {
        if (VERBOSE) std::cout << "find(pattern, reverse_complement=" << reverse_complement << ", skip=" << skip << ")" << std::endl;

        std::cout << "backward searching matching orientation" << std::endl;
        std::vector<usint> pat;
        for (usint i = skip; i < pattern.size(); ++i) {
            pat.push_back(pattern[i]);
        }
        assert(pat.size() > 0);
        find_one_dir(pat);

        std::cout << "backward searching reverse orientation" << std::endl;
        std::vector<usint> revpat;
        for (long long int i = pattern.size() - skip - 1; i >= 0 ; --i) {
            revpat.push_back(pattern[i]);
        }
        assert(revpat.size() > 0);
        find_one_dir(revpat);        

        return pair_type(1,0);;
    }


    static inline unsigned int uint_max(unsigned int a, unsigned int b) {
        return a < b ? b : a;
    }

    pair_type find_one_dir(const std::vector<usint>& pat) const {

        unsigned int myc = (pat[pat.size() - 1] );


            pair_type myinitrange = this->index.getSARange();
            unsigned int delta = get_stddev(myc) * 3;
            std::vector<long unsigned int> hits = this->index.restricted_unique_range_values(myinitrange.first, myinitrange.second, 
                                                                                             myc <= delta ? 1 : myc - delta, // if subtracting results in less than 1, use 1
                                                                                             myc + delta);

            for(std::vector<long unsigned int>::iterator itr = hits.begin(); itr != hits.end(); ++itr) {
                pair_type myrange = /*this->index.getSARange()*/ this->index.getCharRange(*itr);
                myrange.second += 1;

                if (VERBOSE >= 2) std::cout << "bootstrap range for initial query symbol candidate " <<*itr<< " is [" << myrange.first << ".." << myrange.second << "]" << std::endl;

                int deviation = abs(*itr - myc);
                float chi_squared = std::pow((float)deviation / (float)get_stddev(uint_max(*itr, myc)), 2);

                std::vector<long unsigned int> hitvec;
                hitvec.push_back(*itr);

                std::vector<std::pair<long unsigned int, int> > qvec;
                std::pair<long unsigned int, int> q;
                q.first = myc;
                q.second = 0;
                qvec.push_back(q);

                std::vector<pair_type> ranges;



                this->mybackwardSearch(pat, // pattern to search for
                                       pat.size() - 1 , // index of next symbol to search for
                                       myrange, // SA interval
                                       chi_squared, // chi**2 sum
                                       1); // matched count
                
                
            }
            return pair_type(1,0);
        }



    const double MAX_CHISQUARED_CDF = .9;
    const int MAX_LOOKAHEAD = 2;
    //TODO: convert this to recursive call
    //TODO: use sigma*length as stddev

    void report_occurance(unsigned int val) const {
     
        bool backbone = true;
        if (val > this->index.getBackbone()->getSize()-1) {
            backbone = false;
            val -= this->index.getBackbone()->getSize()-1;
        }
        if (!backbone)  std::cout <<"(val-"<<this->index.getBackbone()->getSize()-1 << ")";
        std::cout << val << ", ";
        
        
        // lookup the name
        CSA::DeltaVector::Iterator rmap_iter(rmap_starts);
        unsigned int rmap_num = rmap_iter.rank(val) - 1;
        unsigned int offset = val - rmap_iter.select(rmap_num );
        std::cout << "<" << frag2rmap[rmap_num].second << "+" <<offset << ">" << std::endl;;
        
    }

    bool mybackwardSearch(const std::vector<usint>& pattern,  const unsigned int &it, const pair_type &range, const double &chi_squared_sum, const unsigned int &matched_count) const {
        // handle it=0 to prevent underrun in the other branch
        if (it == 0  || matched_count >= 15) { // stop the recurrsion
            boost::math::chi_squared cs(matched_count);
            double chisqcdf = boost::math::cdf(cs, chi_squared_sum);
            std::cout << "'chi_squared_cdf' : " << chisqcdf << std::endl;

            std::vector<usint>* occurrences = this->index.locateRange(range);
            std::cout << "Found " << occurrences->size() << " match(s) (BWT range " << range.first << ".." << range.second << ") located at: " ;
            assert(occurrences != 0);
            for (std::vector<usint>::iterator mi = occurrences->begin(); mi != occurrences->end(); ++mi) {
                report_occurance(*mi);
            }
            delete occurrences;
        } else {

            int lookahead = MAX_LOOKAHEAD;
            //trim lookahead to max remaining
            if ((int)it - 1 - lookahead < 0) {
                lookahead = it - 1;
            }

            for (int actv_la = 0; actv_la <= lookahead; ++actv_la) {

                // compute the sum of the next lookahead fragments
                unsigned int c = 0;
                for (int j = 0; j <= actv_la; ++j) {
                    int index = it - 1 - j;
                    assert(index >= 0);
                    c += pattern[index];
                }

                //wt stuff
                unsigned int delta = 3 * get_stddev(c) ;
                std::vector<long unsigned int> hits = this->index.restricted_unique_range_values(range.first, range.second, 
                                                                                                 c <= delta ? 1 : c - delta,  // if subtracting results in less than 1, use 1
                                                                                                 c + delta);

                for(std::vector<long unsigned int>::iterator itr = hits.begin(); itr != hits.end(); ++itr) {
                    pair_type new_range = this->index.LF(range, *itr); 

                    int deviation = abs(*itr - c);
                    float chi_squared = std::pow((float)deviation / (float)get_stddev(uint_max(c, *itr)), 2);
                    boost::math::chi_squared cs(matched_count + 1);
                    double chisqcdf = boost::math::cdf(cs, chi_squared_sum + chi_squared);
                    if (chisqcdf <=   MAX_CHISQUARED_CDF) {
                        unsigned int next_search_index = it - 1 - actv_la;
                        this->mybackwardSearch(pattern, next_search_index, new_range, chi_squared_sum + chi_squared, matched_count + 1);
                    }

                }
            }
        }
        return false;
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
    const CSA::DeltaVector& rmap_starts;
    const std::vector<std::pair<unsigned int, std::string> > frag2rmap;
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

    
    void backwardSearch(const std::vector<usint>& pattern, MatchInfo& info) const
    {
        std::cout << "matched is " << info.matched << std::endl;
      sint pos = (info.isReverseComplement() ? info.matched : pattern.size() - info.matched - 1);
      sint adj = (info.isReverseComplement() ? 1 : -1);

      for(; info.matched < pattern.size(); info.matched++, pos += adj)
      {
        unsigned int c = (info.isReverseComplement() ? this->complement(pattern[pos]) : pattern[pos]);

        //wavelet tree stuff

        // std::vector<long unsigned int> hits = this->index.restricted_unique_range_values(info.range.first, info.range.second, 
        //                                                                                  c - DELTA, c + DELTA);
        // std::cout << "DEBUG - wavelet tree query in SA interval [" << info.range.first << ".."<< info.range.second 
        //           << "] has the following symbols within " << DELTA << " alphabet symbols of " << c << ": " ;
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
