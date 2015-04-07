#ifndef PATTERN_CLASSIFIER_H
#define PATTERN_CLASSIFIER_H

#include <fstream>
#include <vector>

#include <misc/definitions.h>


namespace GCSA
{


class PatternClassifier
{
  public:
    PatternClassifier(const std::string& base_name) :
      classify(false),
      forward_file(0), reverse_file(0), notfound_file(0)
    {
      if(base_name.length() == 0) { return; }

      std::string forward_name = base_name + ".forward";
      this->forward_file = new std::ofstream(forward_name.c_str(), std::ios_base::binary);
      if(!*(this->forward_file)) { return; }

      std::string reverse_name = base_name + ".reverse";
      this->reverse_file = new std::ofstream(reverse_name.c_str(), std::ios_base::binary);
      if(!*(this->reverse_file)) { return; }

      std::string notfound_name = base_name + ".notfound";
      this->notfound_file = new std::ofstream(notfound_name.c_str(), std::ios_base::binary);
      if(!*(this->notfound_file)) { return; }

      this->classify = true;
    }

    ~PatternClassifier()
    {
      delete this->forward_file; this->forward_file = 0;
      delete this->reverse_file; this->reverse_file = 0;
      delete this->notfound_file; this->notfound_file = 0;
    }

    std::string stringify(std::vector<usint>& v) {
        std::stringstream s;
        for(std::vector<usint>::iterator itr = v.begin(); itr != v.end(); ++itr)
            s << *itr << ",";
        return s.str();
    }
    void forward(std::vector<usint>& pattern)
    {
      if(this->classify)
      {
          *(this->forward_file) << stringify(pattern) << std::endl;
      }
    }

    void reverse(std::vector<usint>& pattern)
    {
      if(this->classify)
      {
          *(this->reverse_file) << stringify(pattern) << std::endl;
      }
    }

    void notfound(std::vector<usint>& pattern)
    {
      if(this->classify)
      {
          *(this->notfound_file) << stringify(pattern) << std::endl;
      }
    }

  private:
    bool classify;
    std::ofstream* forward_file;
    std::ofstream* reverse_file;
    std::ofstream* notfound_file;
};


} // namespace GCSA


#endif // GCSA_H
