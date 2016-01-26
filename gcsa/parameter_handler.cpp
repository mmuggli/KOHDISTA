#include <iostream>

#include "parameter_handler.h"


//using namespace CSA;


GCSA::ParameterHandler::ParameterHandler(int argc, char** argv, bool _rlcsa, std::string _usage) :
  ok(true), rlcsa(_rlcsa), usage(_usage),
  indels(false), locate(false), penalties(false), reverse_complement(false), verbose(false), write(false), binary_patterns(false), detailed(false), bounded_stddev(false),
  index_name(0), patterns_name(0),
  k(0), skip(0), max_matches(0)
{
  for(int i = 1; i < argc; i++)
  {
    if(argv[i][0] == '-')
    {
      switch(argv[i][1])
      {
      case 'd':
          this->detailed = true; break;
      case 't':
          this->trim = 1; break;
          
        case 'i':
          this->indels = true; break;
        case 'k':
          this->k = atoi(argv[i] + 2); break;
        case 'l':
          this->locate = true; break;
        case 'm':
          this->max_matches = atoi(argv[i] + 2); break;
        case 'p':
          this->penalties = true; break;
        case 'r':
          this->reverse_complement = true; break;
        case 'v':
          this->verbose = true; break;
        case 'w':
          this->write = true; break;
        case 'b':
          this->binary_patterns = true; break;
        case 's':
            if(this->rlcsa) { this->skip = atoi(argv[i] + 2); break; } // Fall through.
      case 'B':
          this->bounded_stddev = true; break;
            
      case 'S':
          this->begin = atoi(argv[i] + 2);
           break; 
      case 'E':
          this->end = atoi(argv[i] + 2);
           break; 
      case 'C':
          this->chi2cdf_thresh = atof(argv[i] + 2);
           break; 
      case 'O':
          this->min_overlap = atoi(argv[i] + 2);
           break; 
      case 'T':
          this->min_t_score = atof(argv[i] + 2);
           break; 
      case 'Z':  // Z for sigma because S and s were taken and Z looks a little like a capital sigma, and you can say zigma
          this->sigma_kbp = atof(argv[i] + 2);
          break;
        default:
          std::cout << "Invalid option: " << argv[i] << std::endl << std::endl;
          this->ok = false;
          return;
      }
    }
    else if(index_name == 0)
    {
      this->index_name = argv[i];
    }
    else if(patterns_name == 0)
    {
      this->patterns_name = argv[i];
    }
  }

  // Sanity checks.
  if(this->k > 0)
  {
    this->reverse_complement = true;
    this->write = false;
    this->skip = 0;
  }
  else
  {
    this->indels = false;
    this->penalties = false;
  }

  if(!(this->locate)) { this->max_matches = 0; }

  if(index_name == 0) { this->ok = false; }
}

void
GCSA::ParameterHandler::printUsage() const
{
  std::cout << this->usage << std::endl;
  std::cout << "  -i   Allow indels (requires -k)." << std::endl;
  std::cout << "  -k#  Allow up to # mismatches (implies -r)." << std::endl;
  std::cout << "  -l   Locate the occurrences." << std::endl;
  std::cout << "  -m#  Try to locate at most # matches (requires -l)." << std::endl;
  std::cout << "  -p   Use different penalties for different edit operations (requires -k)." << std::endl;
  std::cout << "  -r   Search also for the reverse complements of the patterns." << std::endl;
  std::cout << "  -S#  Start pattern to search" << std::endl;
  std::cout << "  -E#  End pattern to search" << std::endl;
  std::cout << "  -C#.#  Chi^2 CDF threshold" << std::endl;
  std::cout << "  -O#  minimum overlap" << std::endl;  
  if(this->rlcsa)
  {
    std::cout << "  -s#  Skip first # characters of the pattern (requires k = 0)." << std::endl;
  }
  std::cout << "  -v   Verbose output." << std::endl;
  std::cout << "  -w   Write found patterns into different files (does not work with -k)." << std::endl;
  std::cout << "  -b   Expect pattern file in binary format." << std::endl;
  std::cout << std::endl;
}

void
GCSA::ParameterHandler::printOptions() const
{
  std::cout << "Options:";
  if(this->indels) { std::cout << " indels"; }
  if(this->locate) { std::cout << " locate"; }
  if(this->max_matches > 0) { std::cout << " max_matches=" << this->max_matches; }
  if(this->k > 0) { std::cout << " mismatches=" << this->k; }
  if(this->penalties) { std::cout << " penalties"; }
  if(this->reverse_complement) { std::cout << " reverse_complement"; }
  if(this->skip > 0) { std::cout << " skip=" << this->skip; }
  if(this->verbose) { std::cout << " verbose"; }
  if(this->write) { std::cout << " write"; }
  std::cout << " min_overlap=" << this->min_overlap;
  std::cout << " chi2cdf_thresh=" << this->chi2cdf_thresh;
  std::cout << " min_t_score=" << this->min_t_score;
  if(this->begin) std::cout << " begin=" << this->begin;
  if(this->end) std::cout << " end=" << this->end;
  std::cout << std::endl;

  std::cout << "Index: " << this->index_name << std::endl;
  if(this->patterns_name != 0)
  {
    std::cout << "Patterns: " << this->patterns_name << std::endl;
  }
  std::cout << std::endl;
}
