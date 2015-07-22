#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

bool sequence_sites_lt(const pair<string, vector<int> >& a, const pair<string, vector<int> >& b) 
{

  return a.second.size() > b.second.size();
}
#define OPTMAP_SIZE 323093

// struct {
//   bool operator()(, int b)
//   {   
//     return a < b;
//   }   
// } sufcmp;

int main() 
{
  char optmap[OPTMAP_SIZE + 2];
  int fd = open("../optmap_bin", 0);
  if (fd == -1)
    printf("ERROR: couldn't open file\n");
  int bytesread = read(fd, optmap, OPTMAP_SIZE);
  printf("Loaded %d bytes from optical map", bytesread);

  std::vector<char*> sa;
  for (int i =0; i<OPTMAP_SIZE; ++i) {
    sa.push_back(optmap[i]);
  }
  sort(sa.begin(), sa.end(), strcmp);
  sdsl::csa_wt<> fm_index;
  sdsl::construct(fm_index, "../optmap_bin", 1);
  //sdsl::store_to_file(fm_index,"fm_index-file.sdsl");
  int placedcount = 0;

  ifstream sequence_file("/s/oak/a/nobackup/muggli/01.GENOME/scaffold/1000contigs.silico");

  map<string, vector<int> > sequence2sites; // in-silico sequence map container
  map<string, int> sequence2sizes;


  string id; int length; int num_sites;
  while(sequence_file >> id >> length >> num_sites) {

    sequence2sizes[id] = length;

    for (int i = 0; i < num_sites; i++) {

      int site; sequence_file >> site;
      sequence2sites[id].push_back(site);
    }
  }
  sequence_file.close();

  vector<pair<string, vector<int> > > sequence_sites; // another in-silico sequence map container
  for (map<string, vector<int> >::iterator it = sequence2sites.begin();
       it != sequence2sites.end(); it++) 
    sequence_sites.push_back(*it);

  sort(sequence_sites.begin(), sequence_sites.end(), sequence_sites_lt); // sort? maybe this makes it run faster


  for (vector<pair<string, vector<int> > >::iterator it = sequence_sites.begin(); it != sequence_sites.end(); it++) {

    sort(it->second.begin(), it->second.end()); //sort insilico frags, WTF?

    // convert offsets of restriction sites to in-silico fragment sizes, only when larger than threshold into sequence_pieces
    vector<int> sequence_pieces; 
    sequence_pieces.push_back(it->second[0]);
    int small_sequence_size = 700;
    int cnt = 0;
    char query[1000];
    for (vector<int>::iterator it2 = it->second.begin()+1; it2 != it->second.end(); it2++) {
      if (*it2-*(it2-1) >= small_sequence_size) {
        sequence_pieces.push_back(*it2-*(it2-1));
        query[cnt] = *it2-*(it2-1) * ((2*2*2*2*2*2*2*2 - 1) / 175) + 1;
        cnt++;
      }
    }

    sequence_pieces.push_back(sequence2sizes[it->first] - it->second[it->second.size()-1]); // put the last piece in
    query[cnt] = 0;

    if (sdsl::count(fm_index,query) > 0) placedcount++;

  }
  std::cout << "placed: " << placedcount << std::endl;
  //  sdsl::write_structure<HTML_FORMAT>(fm_index,"fm_index-file.sdsl.html");
}
