#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <set>
#include <misc/utils.h>
#include "getRSS.c"

#include "gcsa.h"
#include "bwasearch.h"
#include "parameter_handler.h"
#include "pattern_classifier.h"
#include <bits/deltavector.h>

//using namespace CSA;
typedef CSA::usint usint;
typedef CSA::pair_type pair_type;

int main(int argc, char** argv)
{
  double appstart = CSA::readTimer();

  std::cout << "GCSA test" << std::endl;
  std::cout << std::endl;

  GCSA::ParameterHandler handler(argc, argv, false, "Usage: gcsa_test [options] base_name [patterns]");
  if(!(handler.ok)) { handler.printUsage(); return 1; }
  handler.printOptions();


  CSA::DeltaEncoder *rmap_startse = new CSA::DeltaEncoder(1024);
  unsigned int pos = 0;
  // read the frag2rmap file
  std::vector<std::pair<unsigned int, std::string> > frag2rmap;
  std::string f2rm_fname;
  f2rm_fname += handler.index_name; 
  f2rm_fname +=  ".frag2rmap";
  std::cout << "Loading GCSA node to rmap#+offset mapping data from " << f2rm_fname << std::endl;
  std::ifstream f2rm_file(f2rm_fname.c_str());
  unsigned int rmapnum, fragnum;
  std::string rmapname;
  while(f2rm_file >> rmapnum >> fragnum >> rmapname) {
      std::pair<unsigned int, std::string> t;
      t.first = fragnum;
      t.second = rmapname;
      frag2rmap.push_back(t);
      rmap_startse->setBit(fragnum);
      pos += fragnum;
  }
  f2rm_file.close();
  CSA::DeltaVector* rmap_starts = new CSA::DeltaVector(*rmap_startse, pos);  
  CSA::DeltaVector::Iterator itr(*rmap_starts);
// for(int jk = 0; jk < pos; ++jk) {
// if (itr.isSet(jk))  std::cout << "rmap_starts[" << jk << "] is set" << std::endl;
// }

  const GCSA::GCSA gcsa(handler.index_name);
  if(!gcsa.isOk()) { return 2; }
  gcsa.reportSize(true);


  if(handler.patterns_name == 0) { return 0; }
  GCSA::BWASearch<GCSA::GCSA> bwasearch(gcsa, *rmap_starts, frag2rmap, handler);

  std::ifstream patterns(handler.patterns_name, std::ios_base::binary);
  std::cout << "Reading patterns from file: " << handler.patterns_name << std::endl;
  if(!patterns)
  {
    std::cerr << "Error opening pattern file!" << std::endl;
    return 3;
  }
  GCSA::PatternClassifier classifier(handler.write ? handler.patterns_name : "");

  std::vector<std::vector<usint> > rows;
  CSA::readPatternRows(patterns, rows, true, handler.binary_patterns);

  usint total = 0, n = rows.size();
  usint found = 0, forward = 0, reverse = 0;
  usint forward_occurrences = 0, reverse_occurrences = 0;

  usint* edit_distances = new usint[handler.k + 1];
  for(usint i = 0; i <= handler.k; i++) { edit_distances[i] = 0; }
  usint* match_counts = new usint[n];
  for(usint i = 0; i < n; i++) { match_counts[i] = 0; }

  double start = CSA::readTimer();
  std::cout << "Data load time:         " <<  start - appstart  << " seconds" << std::endl;

//  for(usint i = 0/*0*/; i < 2/*n*/; i++)
  if (handler.end != 0) {
      if (handler.end < n) {
          n = handler.end;
      }
  }
  for(usint i = handler.begin; i < n; i++)
  {
    total += rows[i].size();
    //bool match = false;

    // Always start with exact matching.
    double row_start = CSA::readTimer();
    if (rows[i].size() < 10 ) {
        std::cout << "Skipping row " << i << " because size() < 10" << std::endl;
    } else {
        std::string rmap_name = frag2rmap[i*2/*two entries in the automaton/map per every sequence, introduced in valuev2bin.py*/].second;
        std::cout << "### Finding row " << i << "(" << rmap_name<< "): ###" <<std::endl;
        for (std::vector<usint>::iterator ri = rows[i].begin(); ri != rows[i].end(); ++ri) {
            std::cout << *ri << ", ";
        }
        std::cout << std::endl;

        
        bwasearch.find(rows[i], rmap_name);

    
        std::cout << "Find (row = " << i << ") completed in " <<   CSA::readTimer() - row_start << " seconds." << std::endl;
    }
  }
    //disabled as we now run find in backward search
    //    usint temp = bwasearch.handleOccurrences(result, handler.locate, handler.max_matches);
  //   if(temp > 0)
  //   {
  //     forward_occurrences += temp; match_counts[i] += temp;
  //     found++; forward++; match = true;
  //     edit_distances[0]++;
  //     if(handler.write) { classifier.forward(rows[i]); }
  //   }
  //   if(handler.reverse_complement)
  //   {
  //     pair_type reverse_result = bwasearch.find(rows[i], true, handler.skip);
  //     temp = bwasearch.handleOccurrences(reverse_result, handler.locate, handler.max_matches);
  //     if(temp > 0)
  //     {
  //       if(!match)         { found++; edit_distances[0]++; }
  //       if(handler.k == 0) { reverse_occurrences += temp; }
  //       else               { forward_occurrences += temp; }
  //       match_counts[i] += temp;
  //       if(handler.max_matches > 0 && match_counts[i] > handler.max_matches)
  //       {
  //         match_counts[i] = handler.max_matches + 1;
  //       }
  //       reverse++; match = true;
  //       if(handler.write)  { classifier.reverse(rows[i]); }
  //     }
  //   }
  //   if(handler.write && !match) { classifier.notfound(rows[i]); }

  //   // Do approximate matching only if there are no exact matches.
  //   if(handler.k > 0 && !match)
  //   {
  //     std::vector<GCSA::MatchInfo*>* results = bwasearch.find(rows[i], handler.k, handler.indels, handler.penalties);
  //     temp = bwasearch.handleOccurrences(results, handler.locate, handler.max_matches);
  //     if(temp > 0)
  //     {
  //       forward_occurrences += temp; match_counts[i] += temp;
  //       found++;
  //       edit_distances[(*(results->begin()))->errors]++;
  //     }
  //     bwasearch.deleteResults(results);
  //   }

  double time = CSA::readTimer() - start;
  double megabases = total / (double)CSA::MILLION;

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
  std::cout << " (peak RSS " << getPeakRSS() << ")" <<std::endl;

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

