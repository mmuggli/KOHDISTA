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

// valuev stuff
// #include "./../om_set1/msfl.cpp"
// #include "./../om_set1/m_read.cpp"
// #include "./../om_set1/scoring.cpp"
// #include "./../om_set1/alignment.cpp"
#include "./../om_set1/scoring.h"

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
    const float lenwise_s_cutoffs[] = {-0.951583, 0.0793565, 0.43825000000000003, 0.5934475, 0.712202, 0.7915783333333333, 0.8361042857142857, 0.87324875, 0.8987633333333334, 0.9227270000000001, 0.9556363636363637, 0.9755583333333333, 0.9923230769230769};
    const unsigned int expected_s_lut_size = 13;
    
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
    // std::vector<usint> global_pattern;
    // std::set<usint> global_occurrences;
    // std::set<work_t > global_exhausted_nodes;

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
        const int MAX_DEPTH = 100;
        unsigned long long branch_fact_sum[MAX_DEPTH];
        unsigned long long branch_fact_count[MAX_DEPTH];
        for (int i = 0; i < MAX_DEPTH; ++i) {
            branch_fact_sum[i] = 0;
            branch_fact_count[i] = 0;
        }
        scoring_params sp(.2,1.2,.9,3,17.43,0.58, 0.0015, 0.8, 1, 3);//fixme: instantiate just once

        std::cout << "backward searching matching orientation" << std::endl;
        std::set<usint> occurrences;
        for (unsigned int skip = 0; skip < 3; ++skip) {
            std::cout << "skipping " << skip << " query symbols." << std::endl;
            

            std::vector<usint> pat;
            for (usint i = 0; i < pattern.size() - skip; ++i) {
                pat.push_back(pattern[i]);
            }
            assert(pat.size() > 0);
            find_one_dir(pat, occurrences, branch_fact_sum, branch_fact_count, sp);
        }
        if (occurrences.size()) {
            std::cout << "Found " << occurrences.size() << ":" << std::endl;
            for(std::set<usint>::iterator oi = occurrences.begin(); oi != occurrences.end(); ++oi) {
                report_occurrence(*oi);
            }
        }



        std::cout << "backward searching reverse orientation" << std::endl;
        std::set<usint> revoccurrences;
        for (unsigned int skip = 0; skip < 3; ++skip) {
            std::cout << "skipping " << skip << " query symbols." << std::endl;
        

            std::vector<usint> revpat;
            for (long long int i = pattern.size()  - 1; i >= skip ; --i) {
                revpat.push_back(pattern[i]);
            }
            assert(revpat.size() > 0);
            find_one_dir(revpat, revoccurrences, branch_fact_sum, branch_fact_count, sp);        
        }

        if (revoccurrences.size()) {
            std::cout << "Found " << revoccurrences.size() << ":" << std::endl;
            for(std::set<usint>::iterator oi = revoccurrences.begin(); oi != revoccurrences.end(); ++oi) {
                report_occurrence(*oi);
            }
        }

        std::cout << "Branching factor stats:" << std::endl;
        for (int i = 0; i < MAX_DEPTH and branch_fact_count[i]; ++i) {
            std::cout << "depth: " << i << " mean: " << (float)branch_fact_sum[i] / (float)branch_fact_count[i] << " count: " << branch_fact_count[i] << std::endl;
        }
        
     
    }
    const float STDDEV_MULT = 3;

    static inline unsigned int uint_max(unsigned int a, unsigned int b) {
        return a < b ? b : a;
    }

    void local_restricted_unique_range_values(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long d, std::vector<long unsigned int> &ret) const {
        ret = this->index.restricted_unique_range_values(a,b,c,d);
    }


    //fixme: abstract away the differences between the initial and subsequent frag matches, rational: don't repeat yourself
    pair_type find_one_dir(const std::vector<usint>& pattern, std::set<usint> &occurrences, unsigned long long branch_fact_sum[], unsigned long long branch_fact_count[], scoring_params &sp) const {

        int pat_cursor = pattern.size() - 1;

            int lookahead = MAX_LOOKAHEAD;
            if ((int)pat_cursor - 1 - lookahead < 0) {
                lookahead = pat_cursor - 1; //trim lookahead to max remaining, preventing pattern underrun

            }

            for (int actv_la = 0; actv_la <= lookahead; ++actv_la) {
                double search_start = CSA::readTimer();
                std::vector<usint> query_match_frag;

                // compute the sum of the next lookahead fragments
                unsigned int c = 0;
                for (int j = 0; j <= actv_la; ++j) {
                    int index = pat_cursor  - j;
                    assert(index >= 0);
                    c += pattern[index];
                    query_match_frag.push_back(pattern[index]);
                }

//                unsigned int c = (pat[pat.size() - 1] );
        



                pair_type myinitrange = this->index.getSARange();
                unsigned int delta = get_stddev(c) * STDDEV_MULT;
                //double starttime = CSA::readTimer() ;
                std::vector<long unsigned int> hits;
                local_restricted_unique_range_values(myinitrange.first, myinitrange.second, 
                                                     c <= delta ? 1 : c - delta, // if subtracting results in less than 1, use 1
                                                     c + delta,
                                                     hits);
                //double start2time = CSA::readTimer() ;
                // std::set<long unsigned int> hits2 = this->index.array_restricted_unique_range_values(myinitrange.first, myinitrange.second, 
                //                                                                                  c <= delta ? 1 : c - delta, // if subtracting results in less than 1, use 1
                //                                                                                  c + delta);
                //double finishtime = CSA::readTimer() ;
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
                unsigned int depth = 0;
                branch_fact_sum[depth] += hits.size();
                branch_fact_count[depth] += 1;


                for(std::vector<long unsigned int>::iterator hit_itr = hits.begin(); hit_itr != hits.end(); ++hit_itr) {
                    pair_type myrange = /*this->index.getSARange()*/ this->index.getCharRange(*hit_itr);
                    myrange.second += 1; //fixme: what is this for?  is it in original algo?
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


                    std::vector<usint> target_match_frags;
                    target_match_frags.push_back(*hit_itr);
                    std::vector<std::vector<usint>> query_match_frags;
                    query_match_frags.push_back(query_match_frag);
                    std::vector<pair_type> target_match_ranges;
                    target_match_ranges.push_back(myrange);
                    
                    float s_score = sp.opt_size_score((double)subst_frag/1000, (double)c/1000, 1, 1);
                    float expected_avg_s_score = lenwise_s_cutoffs[0];
//                    if (chisqcdf <=   MAX_CHISQUARED_CDF) {
                    if (s_score >= expected_avg_s_score) {
                        
                        this->mybackwardSearch(pattern, // pattern to search for
                                               pattern.size() - 1 - actv_la , // index of next symbol to search for
                                               myrange, // SA interval
                                               chi_squared, // chi**2 sum
                                               1, // matched count
                                               actv_la, // missed count
                                               occurrences,
                                               exhausted_nodes, branch_fact_sum, branch_fact_count, depth + 1/*depth*/,
                                               target_match_frags,
                                               query_match_frags,
                                               target_match_ranges, sp, s_score);
                    }


                }
                std::cout << "find_one_dir (actv_la = " << actv_la << ") completed in " <<   CSA::readTimer() - search_start << " sec." << std::endl;
            }
        return pair_type(1,0);
    }



    const double MAX_CHISQUARED_CDF = .1;
    const int MAX_LOOKAHEAD = 2;
    const unsigned int MIN_MATCH_LEN = 10;
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
        std::cout << "<(rmap #" << rmap_num << ")" << frag2rmap[rmap_num].second << "+" <<offset << ">" << std::endl;;
        
    }


    // [ 10:24.782 ]->[ 5:14.06, 6:20.276 ] s: 7.62454
    bool mybackwardSearch(const std::vector<usint>& pattern, const unsigned int &pat_cursor, const pair_type &range, const double &chi_squared_sum, const unsigned int &matched_count, const unsigned int &missed_count, std::set<usint> &occurrence_set,                     std::set<work_t > &exhausted_nodes, unsigned long long branch_fact_sum[], unsigned long long branch_fact_count[], unsigned int depth,
                          std::vector<usint> &target_match_frags,
                          std::vector<std::vector<usint> > &query_match_frags,
                          std::vector<pair_type> &target_match_ranges, scoring_params &sp, const float &s_score_sum) const {
        // handle pat_cursor=0 to prevent underrun in the other branch
        if (pat_cursor == 0   || (matched_count >= MIN_MATCH_LEN && NU * matched_count - LAMBDA * missed_count >= 8.0)) { // stop the recurrsion
            float t_score = NU * matched_count - LAMBDA * missed_count;
            if (t_score < 8.0) return false;
            boost::math::chi_squared cs(matched_count );
            double chisqcdf = boost::math::cdf(cs, chi_squared_sum);

            std::vector<usint>* occurrences = this->index.locateRange(range);
            if (occurrences->size() > 0) {
                std::cout << "query -> target alignment: " << std::endl;
                std::vector<pair_type>::iterator ri; // range iterator
                std::vector<usint>::iterator fi; // frag iterator
                std::vector<std::vector<usint> >::iterator qi; // query iterator
                float s_tot = 0.0;
                for (qi = query_match_frags.begin(), fi = target_match_frags.begin(), ri = target_match_ranges.begin();
                     qi != query_match_frags.end() &&  fi != target_match_frags.end() && ri != target_match_ranges.end();
                     ++qi, ++fi, ++ri) {
                    std::cout << "[ ";
                    usint querytotal = 0;
                    for (std::vector<usint>::iterator qfpi = qi->begin(); qfpi != qi->end(); ++ qfpi) {
                        std::cout << *qfpi;
                        querytotal += *qfpi;
                        if (qfpi + 1 != qi->end()) {
                            std::cout << ", ";
                        }
                    }
                    std::cout << " ]->[ ";
                    std::cout << *fi << " ";

                    for (usint sai = ri->first; sai <= ri->second && ri->second - ri->first < 6; ++sai) {
                            std::cout << ( (this->index.getBackbone()->contains(sai)) ? "C" : ".");
                            std::cout << ( (this->index.getBackbone()->originalContains(sai)) ? "O" : ".");
//                            std::cout << " " << sai << " ";
                    }

                    std::cout << " :: " ;

                    std::vector<usint>* mr_occurrences = this->index.locateRange(*ri); //match (sa) range
                    const int MAX_RANGES = 5;
                    if (mr_occurrences->size() <= 5) {
                        std::cout << " ";
                        for (std::vector<usint>::iterator mi = mr_occurrences->begin(); mi != mr_occurrences->end(); ++mi) {
                            // std::cout << ( (this->index.getBackbone()->contains(*mi)) ? "C" : ".");
                            // std::cout << ( (this->index.getBackbone()->originalContains(*mi)) ? "O" : ".");
                            std::cout << " " << *mi << " ";
                        }
                    } else {
                        std::cout << " (|SA| = " << mr_occurrences->size() << " > " << MAX_RANGES << ") ";
                    }
                    
                    delete mr_occurrences;

                    std::cout << " ]" ;

                    float s_score = sp.opt_size_score((double)querytotal/1000, (double)*fi/1000, (int)qi->size(), 1);
                    s_tot += s_score;
                    std::cout << " s: " << s_tot << " incr_s: " << s_score << std::endl;

                }
                // std::cout << std::endl;
                // std::cout << "Matched frag sequence in target: " << std::endl;
                // for (std::vector<usint>::iterator fi = target_match_frags.begin(); fi != target_match_frags.end(); ++fi) {
                //     std::cout << *fi << " ";
                // }
                // std::cout << std::endl;
                std::cout << "s-score:" << s_tot << std::endl;
                std::cout << "t-score: " << t_score << std::endl;
                std::cout <<  "chisqcdf: " << chisqcdf << " matches found at: " << std::endl;
                
                for (std::vector<usint>::iterator mi = occurrences->begin(); mi != occurrences->end(); ++mi) {
                    // usint val = *mi;
                    // CSA::DeltaVector::Iterator rmap_iter(rmap_starts);
                    // unsigned int rmap_num = rmap_iter.rank(val) - 1;
                    // unsigned int offset = val - rmap_iter.select(rmap_num );
                    //std::cout << "<(rmap #" << rmap_num << ")" << frag2rmap[rmap_num].second << "+" <<offset << ">" << std::endl;;
                
                    report_occurrence(*mi);


//                    occurrence_set.insert(*mi);
                }
            }
            delete occurrences;
        } else {

            int lookahead = MAX_LOOKAHEAD;
            if ((int)pat_cursor - 1 - lookahead < 0) {
                lookahead = pat_cursor - 1; //trim lookahead to max remaining, preventing pattern underrun

            }

            for (int actv_la = 0; actv_la <= lookahead; ++actv_la) {
                std::vector<usint> query_match_frag;
                // compute the sum of the next lookahead fragments
                unsigned int c = 0;
                for (int j = 0; j <= actv_la; ++j) {
                    int index = pat_cursor - 1 - j;
                    assert(index >= 0);
                    c += pattern[index];
                    query_match_frag.push_back(pattern[index]);
                }

                //wt stuff
                unsigned int delta = STDDEV_MULT * get_stddev(c) ;
                unsigned long long interval_size = range.second - range.first;
                //double start2time = 0, finishtime = 0;//, starttime = CSA::readTimer() ;
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

                branch_fact_sum[depth] += hits2.size();
                branch_fact_count[depth] += 1;
                
                for(std::set<long unsigned int>::iterator hit_itr = hits2.begin(); hit_itr != hits2.end(); ++hit_itr) {
                    // compute chi^2 score for putative substitute fragment in target
                    long unsigned int subst_frag = *hit_itr;

                    int deviation = abs(subst_frag - c);
                    float chi_squared = std::pow((float)deviation / (float)get_stddev(uint_max(c, subst_frag)), 2);
                    //boost::math::chi_squared cs(matched_count + 1);
                    boost::math::chi_squared_distribution<float, boost::math::policies::policy<boost::math::policies::digits10<3> >> cs(matched_count + 1);
                    double chisqcdf = boost::math::cdf(cs, chi_squared_sum + chi_squared);

                    float s_score = sp.opt_size_score((double)c/1000, (double)subst_frag/1000, 1, 1);
                    float avg_s_score = (s_score_sum + s_score) / (depth + 1);
                    float expected_avg_s_score = 1.0;
                    if (depth + 1 <= expected_s_lut_size) {
                        expected_avg_s_score = lenwise_s_cutoffs[depth + 1 /* for this piece */ - 1 /* b/c 0 based indexed */];
                    }
                    if (chisqcdf <=   MAX_CHISQUARED_CDF) {
//                    if (avg_s_score >= expected_avg_s_score) {
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

                        if(exhausted_nodes.count(work) == 0) {
                            target_match_frags.push_back(subst_frag);
                            target_match_ranges.push_back(new_range);
                            query_match_frags.push_back(query_match_frag);
                            this->mybackwardSearch(pattern, next_pat_cursor, new_range, chi_squared_sum + chi_squared, matched_count + 1, missed_count + actv_la + off_backbone_penalty, occurrence_set, exhausted_nodes, branch_fact_sum, branch_fact_count, depth + 1, target_match_frags, query_match_frags, target_match_ranges, sp, s_score_sum + s_score);
                            target_match_frags.pop_back();
                            query_match_frags.pop_back();
                            target_match_ranges.pop_back();
                            exhausted_nodes.insert(work);

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
    //valuev stuff


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
