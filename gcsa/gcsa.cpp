#include <algorithm>
#include <fstream>
#include <iostream>
#include <stack>

#include <map>
#include <misc/utils.h>

#include "gcsa.h"


//TODO find every malloc(CHARS*CHARS) and add a corresponding free()
namespace GCSA
{
const usint CHARS = CSA::CHARS;
const usint MEGABYTE = CSA::MEGABYTE;
const pair_type EMPTY_PAIR = CSA::EMPTY_PAIR;
typedef CSA::uchar uchar;
typedef CSA::pair_type pair_type;
typedef CSA::sint sint;
//--------------------------------------------------------------------------

GCSA::GCSA(const std::string& base_name) :
  ok(false), support_locate(false),
  outgoing(0),
  sampled_positions(0), samples(0),
  alphabet(0),
  backbone(0)
{
//  this->array = (CSA::BitVector**)malloc(CHARS*sizeof(CSA::SDSLVector *));
    //for(usint i = 0; i < CHARS; i++) { this->array[i] = 0; }

  std::string index_name = base_name + GCSA_EXTENSION;
  std::ifstream input(index_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error opening input file (" << index_name << ")!" << std::endl;
    return;
  }
  std::cout << "Loading alphabet from " << index_name << std::endl;
  this->alphabet = new CSA::Alphabet(input);
//  for(usint i = 1; i < 256 /*FIXME:CHARS*/; i++)
  std::cout << "Loading BWT from " << index_name << std::endl;
  for( std::map<usint, pair_type>::const_iterator itr = this->alphabet->begin(); itr != this->alphabet->end(); ++itr)
  {
      usint i = itr->first;
      if (i == 0) continue;
      //if (i % (CHARS/256) == 0) std::cout << "gcsa: processing symbol " << i << " -- (this->array[i] = new CSA::SDSLVector(input);)" << std::endl;
      if(this->alphabet->hasChar(i)) { this->array.populate(i, new CSA::DeltaVector(input)); }

    //else { this->array[i] = 0; }
  }

  array.load(input);
  array.syncFMIndex();
  std::cout << "Loading M from " << index_name << std::endl;
  this->outgoing = new CSA::RLEVector(input);
  this->node_count = this->outgoing->getNumberOfItems();
  std::cout << "Loading sampled_positions for locate from " << index_name << std::endl;
  this->sampled_positions = new CSA::DeltaVector(input);
  usint sample_bits = 0;
  std::cout << "Loading sample_bits for locate from " << index_name << std::endl;
  input.read((char*)&sample_bits, sizeof(usint));
  if(sample_bits > 0)
  {
      std::cout << "Loading samples for locate from " << index_name << std::endl;
    this->samples = new CSA::ReadBuffer(input, this->sampled_positions->getNumberOfItems(), sample_bits);
    this->support_locate = true;
  }
  std::cout << "support_locate is " << this->support_locate << std::endl;
  this->ok = true;

    if(true)
    {
        // std::cout << "Nodes:           " << graph.node_count << std::endl;
        // std::cout << "Edges:           " << graph.edge_count << std::endl;
        // std::cout << "Automata:        " << graph.automata << std::endl;
        std::cout << "Samples:         " << this->sampled_positions->getNumberOfItems() << std::endl;
        std::cout << std::endl;

        this->reportSize(true);
    }


  input.close();
}

GCSA::GCSA(PathGraph& graph, Graph& parent, bool print) :
  ok(false), support_locate(false),
  outgoing(0),
  sampled_positions(0), samples(0),
  alphabet(0),
  backbone(0)
{
    if(graph.status != PathGraph::sorted || !(parent.ok)) { return; }
    //this->array = (CSA::BitVector**)malloc(CHARS*sizeof(CSA::DeltaVector *));
//  DeltaEncoder* array_encoders[CHARS];
    std::map<usint, CSA::DeltaEncoder*> array_encoders;
    //DeltaEncoder** array_encoders = (DeltaEncoder**)malloc(CHARS);
    CSA::RLEEncoder outedges(OUTGOING_BLOCK_SIZE, CSA::MEGABYTE);

    //usint *counts = (usint*)malloc(CHARS*sizeof(usint *));
    std::map<usint, usint> counts;
//  usint counts[CHARS*CHARS];
    // std::cout << "initializing array encoders..." << std::endl;
    // for(usint i = 0; i < CHARS; i++)
    // {
    //     if (i % 1000 == 0) std::cout << i/1000 << " thousand initialized" << std::endl;
    //     this->array[i] = 0;
    //     counts[i] = 0;
    //     array_encoders[i] = new DeltaEncoder(ARRAY_BLOCK_SIZE); // FIXME this uses a lot of memory
    // }
    // std::cout << "done initializing array encoders" << std::endl;
    if(print) { std::cout << "Generating edges... "; std::cout.flush(); }
    if(!graph.generateEdges(parent)) { return; }
    if(print) { std::cout << graph.edge_count << " edges." << std::endl; }


    usint offset = 0, edge_offset = 0, incomingedge_offset = 0;
        
    std::cout << "Writing BWT and M..." << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;

    sdsl::int_vector<> wt_data;
    sdsl::int_vector<1u> inedgetest; 
    inedgetest.resize(graph.edges.size() + 2 /*for the last node?*/);
//    sdsl::bit_vector  = sdsl::bit_vector(graph.edges.size(), 0);
    wt_data.resize(graph.edges.size() + 2);
    std::string inedgebv;

    //FIXME: build up a superflous map for debugging printing
    // std::map<int, int> num2lab;
    // for (std::vector<PathEdge>::iterator edge = graph.edges.begin(); edge != graph.edges.end(); ++edge)
    //     num2lab[edge->from] = edge->label;

    int nodecntr = 0;
    for(std::vector<PathNode>::iterator node = graph.nodes.begin(); node != graph.nodes.end(); ++node)
    {
        nodecntr++;
        if (nodecntr % 1000000 == 0) std::cout << "Processed " << nodecntr << " nodes." << std::endl;

        std::string inedgebvcontr;
        // Write BWT.
        // std::cout << "F[" << offset << "] = " << num2lab[offset] << " node->to: " << node->to << "\tnode->from: " << node->from 
        //      << "\tnode->key.first: "<<node->key.first <<"\tnode->key.second: "<<node->key.second 
        //            <<  "\tL(bwt)[" << offset << "] = " ;
        pair_type edge_range = graph.getEdges(node - graph.nodes.begin(), false);


        for(usint i = edge_range.first; i <= edge_range.second; i++)
        {
            uint label = graph.edges[i].label;
            // std::cout <<"\tedge->label: " << label << "\tedge->from:" << graph.edges[i].from << "\tedge->rank: " << graph.edges[i].rank;
            //if (label == 257) std::cout << "found larger label "<<label <<  std::endl;
            counts[label]++;
            if (array_encoders.find(label) == array_encoders.end()) {
                array_encoders[label] = new CSA::DeltaEncoder(ARRAY_BLOCK_SIZE); // FIXME this uses a lot of memory
            }
            array_encoders[label]->setBit(offset);


            // add WT tracking stuff
            if (i == edge_range.first) {
                inedgebvcontr+="1";  
                inedgetest[incomingedge_offset] = 1;
            } else {

                inedgebvcontr+="0";
                inedgetest[incomingedge_offset] = 0;
            }
            wt_data[incomingedge_offset] = label;
//            std::cout << " setting wt[" << incomingedge_offset << "] = " << label << " ";

            incomingedge_offset++;
        }
        //special case for initial symbol in automaton with no real incoming edges.  We actually add two nodes worth since we want to be able to query one past the final one and then back off one
        if (edge_range.second < edge_range.first) {
            for (int ijk = 0; ijk < 2; ++ijk) {
                inedgebvcontr+="1";  
                inedgetest[incomingedge_offset] = 1;
                //              std::cout << " setting wt[" << incomingedge_offset << "] = " << 0 << " ";
                wt_data[incomingedge_offset] = 0;
                incomingedge_offset++;
            }
        }

//        std::cout << "\t\t M = 1";
        offset++;

        // Write M
        outedges.addBit(edge_offset);
        unsigned addend = std::max((usint)1, (*node).outdegree());
        //for (unsigned ii = 0; ii < addend - 1; ++ii) std::cout <<"0";
        //      std::cout <<" F = " << inedgebvcontr << std::endl;
        edge_offset += addend;
        inedgebv += inedgebvcontr;
    }
//    std::cout << " F = " << inedgebv << std::endl;
    std::cout << "Done.   array_encoders has " << array_encoders.size() << " elements." << std::endl;
    counts[0] = graph.automata; //FIXME: figure out WTF is going on here, c++ static type checking doesn't catch this potential bug
    std::cout << "gcsa: Constructing Alphabet" << std::endl;
    CSA::Alphabet *thealphabet = new CSA::Alphabet(counts);
    this->alphabet = thealphabet;
    std::cout << "gcsa: Constructing array encoders and vectors" << std::endl;
//    for(usint i = 1; i < 256/*FIXME:CHARS*/; i++) //     for(std::map<usint,  pair_type>::iterator mapiter = this->alphabet.begin(); mapiter != this->alphabet.end(); ++mapiter)     
    int symnumcntr = 0;
    for( std::map<usint, pair_type>::const_iterator itr = this->alphabet->begin(); itr != this->alphabet->end(); ++itr)
    {
        symnumcntr++;
        if (symnumcntr % 1000 == 0) std::cout << "Constructed " << symnumcntr << " bit vectors for alphabet symbols." << std::endl;
        usint i = itr->first;
        if (i == 0) continue;
        //if (i % (CHARS/256) == 0) std::cout << "gcsa: processing symbol " << i << " -- (this->array[i] = new CSA::DeltaVector(*(array_encoders[i]), offset))" << std::endl;
        if(this->alphabet->hasChar(i)) {
            if (array_encoders.find(i) == array_encoders.end()) {
                std::cout << "alphabet has " << i << " but array_encoders does not!" << std::endl;
                array_encoders[i] = new CSA::DeltaEncoder(ARRAY_BLOCK_SIZE); // FIXME this uses a lot of memory
            }
            //this->array[i] = new CSA::DeltaVector(*(array_encoders[i]), offset);
            array.populate(i, array_encoders[i], offset);
        }
    }
    std::cout << "gcsa: constructing wavelet tree" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    array.setwt(wt_data);
    std::cout << "flushing outedges" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    outedges.flush();
    std::cout << "gcsa: Constructing outgoing" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    this->outgoing = new CSA::RLEVector(outedges, edge_offset);
    std::cout << "gcsa: Constructing F" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    this->array.constructF(inedgetest);
    std::cout << "gcsa: done Constructing F" << std::endl;
    std::cout << "Elapsed time: " << CSA::readTimer() - CSA::start_time  << std::endl;
    this->node_count = this->outgoing->getNumberOfItems();

//    size_t ones = sdsl::rank_support_v<1>(&inedgetest)(inedgetest.size());
    //sdsl::bit_vector b = inedgetest;
    //std::cout << "inedgetset " <<inedgetest << std::endl;
    sdsl::bit_vector::select_1_type b_sel(&inedgetest);
//    for (unsigned int i = 1; i <= graph.nodes.size(); ++i)
//        std::cout << "sdsl [" << i << "] = " << b_sel(i) << std::endl;


//    std::vector<usint> test = array.restricted_unique_range_values(0,10,1,10000);
    // std::cout << "wt before ser/des: ";
    // for (std::vector<usint>::iterator it = test.begin(); it!=test.end(); ++it) {
    //     std::cout << *it << " ";
    // }
    // std::cout << std::endl;


    // Create the backbone.
    if(parent.backbone != 0)
    {
        std::cout << "gcsa: Constructing backbone" << std::endl;
        this->backbone = new Backbone(*this, graph, parent, print);
    }


    // Sample the graph for locate().
    if(print) { std::cout << "Sampling the graph... "; std::cout.flush();  }
    usint max_sample = 0;
    std::vector<pair_type>* sample_pairs = graph.getSamples(SAMPLE_RATE, max_sample, parent);
    CSA::DeltaEncoder sample_encoder(SAMPLE_BLOCK_SIZE,  CSA::MEGABYTE);
    CSA::WriteBuffer sample_values(sample_pairs->size(), CSA::length(max_sample));
    CSA::parallelSort(sample_pairs->begin(), sample_pairs->end());

    for(usint i = 0; i < sample_pairs->size(); i++)
    {
        sample_encoder.addBit(sample_pairs->at(i).first);
        sample_values.writeItem(sample_pairs->at(i).second);
    }
    sample_encoder.flush();
    delete sample_pairs; sample_pairs = 0;
    std::cout << "gcsa: Constructing sampled positions" << std::endl;
    this->sampled_positions = new CSA::DeltaVector(sample_encoder, offset);
    this->samples = sample_values.getReadBuffer();
    this->support_locate = true;
    if(print)
    {
        std::cout << this->sampled_positions->getNumberOfItems() << " samples." << std::endl;
        std::cout << std::endl;
    }


    if(print)
    {
        std::cout << "Nodes:           " << graph.node_count << std::endl;
        std::cout << "Edges:           " << graph.edge_count << std::endl;
        std::cout << "Automata:        " << graph.automata << std::endl;
        std::cout << "Samples:         " << this->sampled_positions->getNumberOfItems() << std::endl;
        std::cout << std::endl;

        this->reportSize(true);
    }
    this->ok = true;
    //free(array_encoders); //FIXME: switch to hash version
    //free(counts);
}

GCSA::~GCSA()
{
    // for(std::map<usint,  CSA::BitVector*>::iterator mapiter = this->array.begin(); mapiter != this->array.end(); ++mapiter) {
    //     delete mapiter->second; 
    //     mapiter->second = 0;
    // }

  delete this->outgoing; this->outgoing = 0;
  delete this->sampled_positions; this->sampled_positions = 0;
  delete this->samples; this->samples = 0;
  delete this->alphabet; this->alphabet = 0;
  delete this->backbone; this->backbone = 0;
//  free(array);
}

void
GCSA::writeTo(const std::string& base_name) const
{
  std::string index_name = base_name + GCSA_EXTENSION;
  std::ofstream output(index_name.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error opening output file (" << index_name << ")!" << std::endl;
    return;
  }
  std::cout << "Writing alphabet to " << index_name << std::endl;
  this->alphabet->writeTo(output);
  std::cout << "Writing BWT to " << index_name << std::endl;
  this->array.writeTo(output);
  // for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
  // {
  //     if(this->alphabet->hasChar(i)) { this->array.at(i)->writeTo(output); }
  // }
  std::cout << "Writing M to " << index_name << std::endl;
  this->outgoing->writeTo(output);
  std::cout << "Writing sampled_positions to " << index_name << std::endl;
  this->sampled_positions->writeTo(output);
  usint sample_bits = (this->support_locate ? this->samples->getItemSize() : 0);
  std::cout << "Writing sample_bits to " << index_name << std::endl;
  output.write((char*)&sample_bits, sizeof(usint));
  std::cout << "Writing samples to " << index_name << std::endl;
  if(this->support_locate) { this->samples->writeBuffer(output); }

  output.close();

  if(this->backbone != 0)
  {
    this->backbone->writeTo(base_name);
  }
}

bool
GCSA::isOk() const
{
  return this->ok;
}

//--------------------------------------------------------------------------

usint
GCSA::reportSize(bool print) const
{
    usint array_size = this->array.reportSize();
  // for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
  // {
  //   if(this->alphabet->hasChar(i)) { array_size += this->array.at(i)->reportSize(); }
  // }

  usint outedges = this->outgoing->reportSize();

  usint sa_samples = this->sampled_positions->reportSize();
  if(this->support_locate) { sa_samples += this->samples->reportSize(); }

  usint bytes = sizeof(*this) + array_size + outedges + sa_samples + this->alphabet->reportSize();
  if(print)
  {
    std::cout << "Array:           " << (array_size / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Outgoing edges:  " << (outedges / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Samples:         " << (sa_samples / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Total size:      " << (bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << std::endl;
  }

  if(this->backbone != 0)
  {
    if(print) { std::cout << "Backbone information:" << std::endl; }
    bytes += this->backbone->reportSize(print);
  }

  return bytes;
}

void
GCSA::printCounts() const
{
  for(usint i = 0; i < CHARS; i++)
  {
    if(this->alphabet->hasChar(i))
    {
      std::cout << "counts[" << i << " (" << ((uchar)i) << ")] = " << this->alphabet->countOf(i) << std::endl;
    }
  }
  std::cout << std::endl;
}

//--------------------------------------------------------------------------

pair_type
GCSA::find(const std::string& pattern) const
{
  if(pattern.length() == 0) { return this->getSARange(); }

  sint pos = pattern.length() - 1;
  pair_type range = this->getCharRange(pattern[pos]);
  if(CSA::isEmpty(range)) { return range; }

  for(pos--; pos >= 0; pos--)
  {
    range = this->LF(range, pattern[pos]);
    if(CSA::isEmpty(range)) { return range; }
  }
  return range;
}

//--------------------------------------------------------------------------

std::vector<usint>*
GCSA::locate(pair_type range) const
{
  if(!(this->support_locate)) { return 0; }
  std::vector<usint>* vec = new std::vector<usint>;

  this->locate(range, *vec);
  CSA::removeDuplicates(vec, false);

  return vec;
}

std::vector<usint>*
GCSA::locate(std::vector<pair_type>& ranges) const
{
  if(!(this->support_locate)) { return 0; }
  std::vector<usint>* vec = new std::vector<usint>;

  for(std::vector<pair_type>::iterator iter = ranges.begin(); iter != ranges.end(); ++iter)
  {
    this->locate(*iter, *vec);
  }
  CSA::removeDuplicates(vec, false);

  return vec;
}

bool
GCSA::isSampled(usint index) const
{
  if(!(this->support_locate) || index >= this->getSize())  { return false; }

  CSA::DeltaVector::Iterator sample_iter(*(this->sampled_positions));
  return sample_iter.isSet(index);
}

usint
GCSA::locate(usint index) const
{
  if(!(this->support_locate) || index >= this->getSize())  { return this->getSize(); }

  return this->locateUnsafe(index);
}

usint
GCSA::labelOf(usint index) const
{
  if(index >= this->getSize()) { return 0; }
  char iterbuf[this->outgoing->iterSize()];
  CSA::BitVector::Iterator* outgoing_iter = this->outgoing->newIterator(iterbuf);
  index = outgoing_iter->select(index);
  outgoing_iter->~Iterator();
  return this->alphabet->charAt(index);
}

//--------------------------------------------------------------------------

pair_type
GCSA::getSARange() const
{
  return pair_type(0, this->getSize() - 1);
}

pair_type
GCSA::getBWTRange() const
{
  return pair_type(0, this->getSize() - 1);
}

pair_type
GCSA::getCharRange(usint c) const
{
  pair_type range = this->alphabet->getRange(c);
  if(CSA::isEmpty(range)) { return range; }
  return this->convertToNodeRange(range);
}

void
GCSA::convertToBWTRange(pair_type& sa_range) const
{
}

void
GCSA::convertToSARange(pair_type& bwt_range) const
{
}

void
GCSA::convertToSARange(std::vector<pair_type>& bwt_ranges) const
{
}


    std::vector<usint> GCSA::restricted_unique_range_values(usint l, usint r, usint min, usint max) const
    {
        return array.restricted_unique_range_values(l, r, min, max);
    }

pair_type
GCSA::LF(pair_type range, usint c) const
{
    if (VERBOSE >= 2) std::cout << "LF(<" <<range.first << "," << range.second << ">, " << c <<") = <";
    
    if(!(this->alphabet->hasChar(c))) { return EMPTY_PAIR; }

    // Follow edges backward using BWT.
    char iterbufc[this->array.iterSize(c)];
    CSA::CharVector::Iterator* array_iter = this->array.newIterator(c, iterbufc);
  
    range.first = this->alphabet->cumulative(c) + array_iter->rank(range.first, true) - 1;

    range.second = this->alphabet->cumulative(c) + array_iter->rank(range.second) - 1;

    if (VERBOSE >= 2)   std::cout << range.first <<  "," << range.second << ">" << std::endl;
    if(CSA::isEmpty(range)) { return EMPTY_PAIR; }
    array_iter->~Iterator();
    pair_type ret_range = this->convertToNodeRange(range);
    if (VERBOSE >= 2)   std::cout << " --->  <" << ret_range.first <<  "," << ret_range.second << ">" << std::endl;
    return ret_range;
}

std::vector<usint>*
GCSA::locateRange(pair_type range) const
{
  return this->locate(range);
}

std::vector<usint>*
GCSA::locateRanges(std::vector<pair_type>& ranges) const
{
  return this->locate(ranges);
}

//--------------------------------------------------------------------------

std::vector<usint>*
GCSA::getSuccessors(usint index) const
{
  if(index >= this->getSize()) { return 0; }

  std::vector<usint>* result = new std::vector<usint>;
  if(index == 0) { return result; } // Final node.
  char iterbuf[this->outgoing->iterSize()];
  CSA::BitVector::Iterator* outgoing_iter = this->outgoing->newIterator(iterbuf);
  index = outgoing_iter->select(index);
  usint successors = outgoing_iter->selectNext() - index;
  usint c = this->alphabet->charAt(index);

  // Find the corresponding incoming edges using BWT.
  char iterbufc[this->array.iterSize(c)];
  CSA::CharVector::Iterator* array_iter = this->array.newIterator(c, iterbufc);
  result->push_back(array_iter->select(index - this->alphabet->cumulative(c)));
  for(usint i = 1; i < successors; i++) { result->push_back(array_iter->selectNext()); }
  outgoing_iter->~Iterator();
  array_iter->~Iterator();
  return result;
}

CSA::CharVector::Iterator*
GCSA::getIterator(usint c, char* placement) const
{
  if(c < CHARS && c > 0 && this->alphabet->hasChar(c)) { return this->array.newIterator(c, placement); }
  return 0;
}

CSA::BitVector::Iterator*
GCSA::getEdgeIterator(char *placement) const
{
    return this->outgoing->newIterator(placement);
}

 size_t GCSA::edgeIterSize() const
 {
     return this->outgoing->iterSize();
 }

//--------------------------------------------------------------------------

pair_type
GCSA::convertToNodeRange(pair_type edge_range) const
{
  char iterbuf[this->outgoing->iterSize()];
  CSA::BitVector::Iterator* outgoing_iter = this->outgoing->newIterator(iterbuf);
  edge_range.first = outgoing_iter->rank(edge_range.first) - 1;
  edge_range.second = outgoing_iter->rank(edge_range.second) - 1;
  outgoing_iter->~Iterator();
  return edge_range;
}

usint
GCSA::Psi(usint index) const
{
  if(index == 0) { return this->getSize() - 1; } // Final node.

  char iterbuf[this->outgoing->iterSize()];
  CSA::BitVector::Iterator* outgoing_iter = this->outgoing->newIterator(iterbuf);
  index = outgoing_iter->select(index);
  usint c = this->alphabet->charAt(index);

  // Find the corresponding incoming edge using BWT.
  char iterbufc[this->array.iterSize(c)];
  CSA::CharVector::Iterator* array_iter = this->array.newIterator(c, iterbufc);
  index = array_iter->select(index - this->alphabet->cumulative(c));
  outgoing_iter->~Iterator();
  array_iter->~Iterator();
  return index;
}

usint
GCSA::LF(usint index, usint c) const
{
  char iterbufc[this->array.iterSize(c)];
	CSA::CharVector::Iterator* array_iter = this->array.newIterator(c, iterbufc);
  index = this->alphabet->cumulative(c) + array_iter->rank(index) - 1;
  array_iter->~Iterator();
  char iterbuf[this->outgoing->iterSize()];
  CSA::BitVector::Iterator* edge_iter = this->outgoing->newIterator(iterbuf);
  index = edge_iter->rank(index) - 1;
  edge_iter->~Iterator();
  return index;
}

//--------------------------------------------------------------------------

void
GCSA::locate(pair_type range, std::vector<usint>& vec) const
{
  range.second = std::min(range.second, this->getSize() - 1);
  if(CSA::isEmpty(range)) { return; }

  for(usint pos = range.first; pos <= range.second; pos++)
  {
    vec.push_back(this->locateUnsafe(pos));
  }
}

usint
GCSA::locateUnsafe(usint index) const
{
	CSA::DeltaVector::Iterator sample_iter(*(this->sampled_positions));
  usint temp = index, steps = 0;

  while(!(sample_iter.isSet(temp))) { temp = this->Psi(temp); steps++; }
  if (VERBOSE >= 3) std::cout <<"locateUnsafe(" << index << ") steps: " << steps << "rank: " << sample_iter.rank(temp) << std::endl;
  return this->samples->readItem(sample_iter.rank(temp) - 1) - steps;
}

//--------------------------------------------------------------------------


Backbone::Backbone(const GCSA& _gcsa, PathGraph& graph, Graph& parent, bool print) :
  gcsa(_gcsa),
  nodes(0), edges(0), original(0),
  size(0), ok(false)
{
  if(graph.status != PathGraph::ready || parent.backbone == 0) { return; }

  if(print)
  {
    std::cout << "Compressing backbone... "; std::cout.flush();
  }

  CSA::SuccinctEncoder original_encoder(NODE_BLOCK_SIZE, CSA::MEGABYTE);
  CSA::SuccinctVector::Iterator iter(*(parent.backbone));
  for(usint i = 0; i < graph.node_count; i++)
  {
    if(iter.isSet(graph.nodes[i].from)) { original_encoder.addBit(i); }
  }
  original_encoder.flush();
  this->original = new CSA::SuccinctVector(original_encoder, graph.node_count);
  if(print)
  {
    std::cout << "original(" << this->original->getNumberOfItems() << ") ";
    std::cout.flush();
  }



  // Use a kind of backward searching to find the backbone in PathGraph.
  usint bb_nodes = 0;
  const unsigned int num_automata = this->gcsa.getNumberOfAutomata();
  for(usint j = 0; j < num_automata; j++) {
      std::cout << "processing automata " << j << " of " <<  num_automata  << std::endl;
      usint curr_node = j;  // Final node of automaton j.
      graph.nodes[curr_node].setBackbone(); 
      bb_nodes++;
      while(curr_node < this->gcsa.getSize() - num_automata) {// Not the initial node.
          //std::cout << "\tcurr_node: " << curr_node << " of " << this->gcsa.getSize() - num_automata<< std::endl;
          bool found = false;
          pair_type edge_range = graph.getEdges(curr_node, false);
          for(usint i = edge_range.first; i <= edge_range.second; i++) {
              usint prev = graph.edges[i].from;
              //std::cout << "\t\tedge " << i << " of " << edge_range.second << " (prev " << prev << ")" <<std::endl;
              //std::cout << "\t\tpredicate: " << iter.isSet(graph.nodes[prev].from) << " " << graph.nodes[prev].value() << " " <<  graph.nodes[curr_node].value() - 1 << std::endl;
              if(iter.isSet(graph.nodes[prev].from) && graph.nodes[prev].value() == graph.nodes[curr_node].value() - 1) {
                  graph.nodes[prev].setBackbone(); 
                  bb_nodes++;
                  curr_node = prev;
                  found = true; break;
              }
          }
          if(!found) {
              std::cerr << "Error: Cannot find previous backbone node!" << std::endl;
              return;
          }
      }
  }
  if(print)
  {
    std::cout << "found(" << bb_nodes << ") "; std::cout.flush();
  }


  // Scan the graph forward and build backbone information.
  CSA::SuccinctEncoder node_encoder(NODE_BLOCK_SIZE, CSA::MEGABYTE);
  CSA::RLEEncoder edge_encoder(EDGE_BLOCK_SIZE, CSA::MEGABYTE);
  char iterbuf[this->gcsa.outgoing->iterSize()];
  std::cout << "edge iter buffer size is " << this->gcsa.outgoing->iterSize() << std::endl;
  CSA::BitVector::Iterator* edge_iter = this->gcsa.outgoing->newIterator(iterbuf);
  usint offset = edge_iter->select(0);

  graph.sortEdges(true, true);
  for(usint i = 0; i < graph.node_count; i++)
  {
    if(graph.nodes[i].isBackbone())
    {
      node_encoder.addBit(i); this->size++;
      if(i >= num_automata)
      {
        bool found = false;
        pair_type successors = graph.getEdges(i, true);
        for(usint j = successors.first; j <= successors.second; j++)
        {
          usint to = graph.edges[j].rank;
          if(graph.nodes[to].isBackbone() && graph.nodes[to].value() == graph.nodes[i].value() + 1)
          {
            offset += j - successors.first;
            found = true; break;
          }
        }
        if(!found)
        {
          std::cerr << "Error: Cannot find next backbone node!" << std::endl;
          std::cerr << "       i = " << i << ", from = " << graph.nodes[i].from << ", value = "
                    << graph.nodes[i].value() << std::endl;
          return;
        }
      }
    }
    edge_encoder.addBit(offset);
    offset = edge_iter->selectNext();
  }

  node_encoder.flush(); edge_encoder.flush();
  this->nodes = new CSA::SuccinctVector(node_encoder, graph.node_count);
  this->edges = new CSA::RLEVector(edge_encoder, this->gcsa.outgoing->getSize());


  if(print)
  {
    std::cout << "done." << std::endl;
  }
  this->ok = true;
  std::cout << "edge_iter addr " << edge_iter << std::endl;
  edge_iter->~Iterator(); //FIXME: leak, but avoid crash
}

Backbone::Backbone(const std::string& base_name, const GCSA& _gcsa) :
  gcsa(_gcsa),
  nodes(0), edges(0), original(0),
  ok(false)
{
  std::string backbone_name = base_name + BACKBONE_EXTENSION;
  std::ifstream input(backbone_name.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "Error opening input file (" << backbone_name << ")!" << std::endl;
    return;
  }

  input.read((char*)&(this->size), sizeof(this->size));
  this->nodes = new CSA::SuccinctVector(input);
  this->edges = new CSA::RLEVector(input);
  this->original = new CSA::SuccinctVector(input);

  this->ok = true;
  input.close();
}

Backbone::~Backbone()
{
  delete this->nodes; this->nodes = 0;
  delete this->edges; this->edges = 0;
  delete this->original; this->original = 0;
}

void
Backbone::writeTo(const std::string& base_name) const
{
  std::string backbone_name = base_name + BACKBONE_EXTENSION;
  std::ofstream output(backbone_name.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "Error opening output file (" << backbone_name << ")!" << std::endl;
    return;
  }

  output.write((char*)&(this->size), sizeof(this->size));
  this->nodes->writeTo(output);
  this->edges->writeTo(output);
  this->original->writeTo(output);

  output.close();
}

usint
Backbone::getSize() const
{
  return this->size;
}

bool
Backbone::isOk() const
{
  return this->ok;
}

usint
Backbone::reportSize(bool print) const
{
  usint nodes_size = this->nodes->reportSize();
  usint edges_size = this->edges->reportSize();
  usint original_size = this->original->reportSize();

  usint bytes = sizeof(*this) + nodes_size + edges_size + original_size;
  if(print)
  {
    std::cout << "Nodes:           " << (nodes_size / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Edges:           " << (edges_size / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Original:        " << (original_size / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << "Total size:      " << (bytes / (double)MEGABYTE) << " MB" << std::endl;
    std::cout << std::endl;
  }

  return bytes;
}

//--------------------------------------------------------------------------

bool
Backbone::contains(usint index) const
{
	CSA::SuccinctVector::Iterator iter(*(this->nodes));
  return iter.isSet(index);
}

bool
Backbone::originalContains(usint index) const
{
	CSA::SuccinctVector::Iterator iter(*(this->original));
  return iter.isSet(index);
}

usint
Backbone::next(usint index) const
{
  if(index < this->gcsa.getNumberOfAutomata())  // From final node to initial node.
  {
    return this->gcsa.getSize() - this->gcsa.getNumberOfAutomata() + index;
  }

  CSA::RLEVector::Iterator iter(*(this->edges));
  index = iter.select(index);
  usint c = this->gcsa.alphabet->charAt(index);

  // Find the corresponding incoming edge using BWT.
  char iterbufc[this->gcsa.array.iterSize(c)];
  CSA::CharVector::Iterator* array_iter = this->gcsa.array.newIterator(c, iterbufc);
  index = array_iter->select(index - this->gcsa.alphabet->cumulative(c));
  array_iter->~Iterator();
  return index;
}

usint
Backbone::previous(usint index) const
{
  if(index >= this->gcsa.getSize() - this->gcsa.getNumberOfAutomata())  // From initial node to final node.
  {
    return index - (this->gcsa.getSize() - this->gcsa.getNumberOfAutomata());
  }

  CSA::Alphabet* alpha = this->gcsa.alphabet;
  for(usint i = 0; i < alpha->getAlphabetSize(); i++)
  {
    usint c = alpha->getTextChar(i);

    // If BWT[index] contains c, follow the corresponding edge backward.
  char iterbufc[this->gcsa.array.iterSize(c)];
    CSA::CharVector::Iterator* array_iter = this->gcsa.array.newIterator(c, iterbufc);
    if(!(array_iter->isSet(index))) { continue; }
    index = array_iter->rank(index) - 1;

    // If we followed a backbone edge, return the node we ended up in.
    CSA::RLEVector::Iterator edge_iter(*(this->edges));
    if(!(edge_iter.isSet(index))) { continue; }
    array_iter->~Iterator();
    return edge_iter.rank(index) - 1;
  }

  return 0;
}

//--------------------------------------------------------------------------

} // namespace GCSA
