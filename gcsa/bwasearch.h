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
    typedef std::pair<unsigned int, pair_type> work_t;
    std::vector<usint> global_pattern;
    std::set<usint> global_occurrences;
    std::set<work_t > global_exhausted_nodes;
    
template<class Index>
class BWASearch
{
  public:

    const static usint INITIAL_STEP = 8;  // Should be >= log(n) / 4.

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



    void find(const std::vector<usint>& pattern) const {



        std::cout << "backward searching matching orientation" << std::endl;
        std::set<usint> occurrences;
        global_occurrences = occurrences;
        for (unsigned int skip = 0; skip < 3; ++skip) {
        


            std::vector<usint> pat;
            for (usint i = 0; i < pattern.size() - skip; ++i) {
                pat.push_back(pattern[i]);
            }
            assert(pat.size() > 0);
            find_one_dir(pat, occurrences);
        }
        if (global_occurrences.size()) {
            std::cout << "Found " << global_occurrences.size() << ":" << std::endl;
            for(std::set<usint>::iterator oi = global_occurrences.begin(); oi != global_occurrences.end(); ++oi) {
                report_occurrence(*oi);
            }
        }



        std::cout << "backward searching reverse orientation" << std::endl;
        std::set<usint> revoccurrences;
        global_occurrences = revoccurrences;
        for (unsigned int skip = 0; skip < 3; ++skip) {
        
        

            std::vector<usint> revpat;
            for (long long int i = pattern.size()  - 1; i >= skip ; --i) {
                revpat.push_back(pattern[i]);
            }
            assert(revpat.size() > 0);
            find_one_dir(revpat, revoccurrences);        
        }

        if (global_occurrences.size()) {
            std::cout << "Found " << global_occurrences.size() << ":" << std::endl;
            for(std::set<usint>::iterator oi = global_occurrences.begin(); oi != global_occurrences.end(); ++oi) {
                report_occurrence(*oi);
            }
        }

     
    }
    const float STDDEV_MULT = 3;

    static inline unsigned int uint_max(unsigned int a, unsigned int b) {
        return a < b ? b : a;
    }

    void local_restricted_unique_range_values(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long d, std::vector<long unsigned int> &ret) const {
        ret = this->index.restricted_unique_range_values(a,b,c,d);
    }


    
    pair_type find_one_dir(const std::vector<usint>& pattern, std::set<usint> &occurrences) const {

        int pat_cursor = pattern.size() - 1;

            int lookahead = MAX_LOOKAHEAD;
            if ((int)pat_cursor - 1 - lookahead < 0) {
                lookahead = pat_cursor - 1; //trim lookahead to max remaining, preventing pattern underrun

            }

            for (int actv_la = 0; actv_la <= lookahead; ++actv_la) {
                double search_start = CSA::readTimer();


                // compute the sum of the next lookahead fragments
                unsigned int c = 0;
                for (int j = 0; j <= actv_la; ++j) {
                    int index = pat_cursor  - j;
                    assert(index >= 0);
                    c += pattern[index];
                }

//                unsigned int c = (pat[pat.size() - 1] );
        



                pair_type myinitrange = this->index.getSARange();
                unsigned int delta = get_stddev(c) * STDDEV_MULT;
                double starttime = CSA::readTimer() ;
                std::vector<long unsigned int> hits;
                local_restricted_unique_range_values(myinitrange.first, myinitrange.second, 
                                                     c <= delta ? 1 : c - delta, // if subtracting results in less than 1, use 1
                                                     c + delta,
                                                     hits);
                double start2time = CSA::readTimer() ;
                // std::set<long unsigned int> hits2 = this->index.array_restricted_unique_range_values(myinitrange.first, myinitrange.second, 
                //                                                                                  c <= delta ? 1 : c - delta, // if subtracting results in less than 1, use 1
                //                                                                                  c + delta);
                double finishtime = CSA::readTimer() ;
//                std::cout << "wt search = " << start2time - starttime << "; array scan = " << finishtime - start2time << "; interval = " << myinitrange.second - myinitrange.first << "; query size = " << c << std::endl;
                // ensure new methods gives equivalent results
                // for(std::set<long unsigned int>::iterator h2i = hits2.begin(); h2i != hits2.end(); ++h2i) {
                //     if (std::binary_search(hits.begin(), hits.end(), *h2i)) {
                //         ;//std::cout << "Found " << needle << '\n';
                //     } else {
                //         std::cout << "Alert! Array search found the following element not found by the wt: " << *h2i << std::endl;
                //     }

                // }

                // if (hits.size() != hits2.size() ) {
                //     std::cout << "Alert! wt is of size " << hits.size() << " and array results are of size " << hits2.size() << std::endl;
                // }

                std::set<work_t > exhausted_nodes;
                    
                for(std::vector<long unsigned int>::iterator hit_itr = hits.begin(); hit_itr != hits.end(); ++hit_itr) {
                    pair_type myrange = /*this->index.getSARange()*/ this->index.getCharRange(*hit_itr);
                    myrange.second += 1;
                    long unsigned int subst_frag = *hit_itr;
                    if (VERBOSE >= 2) std::cout << "bootstrap range for initial query symbol candidate " <<subst_frag<< " is [" << myrange.first << ".." << myrange.second << "]" << std::endl;

                    int deviation = abs(subst_frag - c);
                    float chi_squared = std::pow((float)deviation / (float)get_stddev(uint_max(subst_frag, c)), 2);

                    std::vector<long unsigned int> hitvec;
                    hitvec.push_back(subst_frag);

                    std::vector<std::pair<long unsigned int, int> > qvec;
                    std::pair<long unsigned int, int> q;
                    q.first = c;
                    q.second = 0;
                    qvec.push_back(q);

                    std::vector<pair_type> ranges;


                    global_pattern = pattern;
                    global_exhausted_nodes=exhausted_nodes;
                    this->mybackwardSearch(//pattern, // pattern to search for
                                           pattern.size() - 1 - actv_la , // index of next symbol to search for
                                           myrange, // SA interval
                                           chi_squared, // chi**2 sum
                                           1, // matched count
                                           actv_la); // missed count
                                           //occurrences,
                                           //exhausted_nodes); 
                
                
                }
                std::cout << "find_one_dir (actv_la = " << actv_la << ") completed in " <<   CSA::readTimer() - search_start << " sec." << std::endl;
            }
        return pair_type(1,0);
    }



    const double MAX_CHISQUARED_CDF = .9;
    const int MAX_LOOKAHEAD = 2;
    const int MIN_MATCH_LEN = 10;
    const float LAMBDA = 1.2;
    const float NU = 0.9;
    //TODO: convert this to recursive call
    //TODO: use sigma*length as stddev

    void report_occurrence(unsigned int val) const {
     
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

    bool mybackwardSearch(/*const std::vector<usint>& pattern,*/  const unsigned int &pat_cursor, const pair_type &range, const double &chi_squared_sum, const unsigned int &matched_count, const unsigned int &missed_count/*, std::set<usint> &occurrence_set,                     std::set<work_t > &exhausted_nodes*/) const {
        // handle pat_cursor=0 to prevent underrun in the other branch
        if (pat_cursor == 0   || (matched_count >= MIN_MATCH_LEN && NU * matched_count - LAMBDA * missed_count >= 8.0)) { // stop the recurrsion
            float t_score = NU * matched_count - LAMBDA * missed_count;
            if (t_score < 8.0) return false;
            boost::math::chi_squared cs(matched_count );
            double chisqcdf = boost::math::cdf(cs, chi_squared_sum);

            std::vector<usint>* occurrences = this->index.locateRange(range);
            for (std::vector<usint>::iterator mi = occurrences->begin(); mi != occurrences->end(); ++mi) {
                //std::cout << "t-score: " << t_score <<  " chisqcdf: " << chisqcdf << " match found at: ";
                //report_occurrence(*mi);
                global_occurrences.insert(*mi);
            }
            delete occurrences;
        } else {

            int lookahead = MAX_LOOKAHEAD;
            if ((int)pat_cursor - 1 - lookahead < 0) {
                lookahead = pat_cursor - 1; //trim lookahead to max remaining, preventing pattern underrun

            }

            for (int actv_la = 0; actv_la <= lookahead; ++actv_la) {

                // compute the sum of the next lookahead fragments
                unsigned int c = 0;
                for (int j = 0; j <= actv_la; ++j) {
                    int index = pat_cursor - 1 - j;
                    assert(index >= 0);
                    c += global_pattern[index];
                }

                //wt stuff
                unsigned int delta = STDDEV_MULT * get_stddev(c) ;
                unsigned long long interval_size = range.second - range.first;
                double start2time = 0, finishtime = 0;//, starttime = CSA::readTimer() ;
                const int WT_MIN = 750;
                
                std::vector<long unsigned int> hits;
                std::set<long unsigned int> hits2;

                    // hits = this->index.restricted_unique_range_values(range.first, range.second, 
                    //                                                   c <= delta ? 1 : c - delta,  // if subtracting results in less than 1, use 1
                    //                                                   c + delta);

                
                if (interval_size > WT_MIN) {
                    hits = this->index.restricted_unique_range_values(range.first, range.second, 
                                                                      c <= delta ? 1 : c - delta,  // if subtracting results in less than 1, use 1
                                                                      c + delta);
                    //FIXME: find a way not to copy the WT results into a set! This is the uncommon case, and it's all on the stack, but still wasteful
                    for(std::vector<long unsigned int>::iterator h2i = hits.begin(); h2i != hits.end(); ++h2i)
                        hits2.insert(*h2i);
                } else {
                    
                    //start2time = CSA::readTimer() ;
                    this->index.array_restricted_unique_range_values(range.first, range.second, 
                                                                             c <= delta ? 1 : c - delta, // if subtracting results in less than 1, use 1
                                                                     c + delta, hits2);
                }
                //finishtime = CSA::readTimer() ;
//                if (interval_size > MIN_BOTH) std::cout << "wt search = " << start2time - starttime << "; array scan = " << finishtime - start2time << "; interval = " << interval_size << "; query size = " << c << std::endl;



                // for(std::set<long unsigned int>::iterator h2i = hits2.begin(); h2i != hits2.end(); ++h2i) {
                //     if (std::binary_search(hits.begin(), hits.end(), *h2i)) {
                //         ;//std::cout << "Found " << needle << '\n';
                //     } else {
                //         std::cout << "Alert! Array search found the following element not found by the wt: " << *h2i << std::endl;
                //     }

                // }

                // if (hits.size() != hits2.size() ) {
                //     std::cout << "Alert! wt is of size " << hits.size() << " and array results are of size " << hits2.size() << std::endl;
                // }

                
                for(std::set<long unsigned int>::iterator hit_itr = hits2.begin(); hit_itr != hits2.end(); ++hit_itr) {
                    // compute chi^2 score for putative substitute fragment in target
                    long unsigned int subst_frag = *hit_itr;
                    int deviation = abs(subst_frag - c);
                    float chi_squared = std::pow((float)deviation / (float)get_stddev(uint_max(c, subst_frag)), 2);
                    //boost::math::chi_squared cs(matched_count + 1);
                    boost::math::chi_squared_distribution<float, boost::math::policies::policy<boost::math::policies::digits10<3> >> cs(matched_count + 1);
                    double chisqcdf = boost::math::cdf(cs, chi_squared_sum + chi_squared);

                    if (chisqcdf <=   MAX_CHISQUARED_CDF) {
                        pair_type new_range = this->index.LF(range, subst_frag); 
                        unsigned int next_pat_cursor = pat_cursor - 1 - actv_la;

                        bool off_backbone_penalty = 0;
                        // if (new_range.first == new_range.second) {
                        //     std::vector<usint>* occurrences = this->index.locateRange(new_range);
                        //     if (*occurrences->begin() > this->index.getBackbone()->getSize()-1) {
                        //         off_backbone_penalty++; //FIXME can we figure out order 2 skip nodes here?
                        //     }

                        // }

                        work_t work(next_pat_cursor, new_range);

                        if(global_exhausted_nodes.count(work) == 0) {
                            this->mybackwardSearch(/*pattern,*/ next_pat_cursor, new_range, chi_squared_sum + chi_squared, matched_count + 1, missed_count + actv_la + off_backbone_penalty/*, occurrence_set, exhausted_nodes*/);
                            global_exhausted_nodes.insert(work);

                        }
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
