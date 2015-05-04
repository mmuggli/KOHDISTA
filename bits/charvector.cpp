#include "charvector.h"
#include "deltavector.h"
#include <iostream>

namespace CSA
{

    usint CharVector::maxlength() const
    {
        usint maxl = 0;
//        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        for(std::map<usint, CSA::BitVector*>::const_iterator itr = array.begin(); itr != array.end(); ++itr)
        {
            usint i = itr->first;
            if (i == 0) continue;
            if(array.count(i) && array.at(i)->getSize() > maxl) { maxl =  array.at(i)->getSize();}
        }
        return maxl;
    }

    std::vector<usint> CharVector::restricted_unique_range_values(usint l, usint r, usint min, usint max) const //FIXME: don't copy vector in return
    {
        std::vector<long unsigned int> hits;
        // const sdsl::csa_wt<sdsl::wt_int<>, 64, 64, sdsl::sa_order_sa_sampling<>, 
        //                    sdsl::int_vector<>, 
        //                    sdsl::int_alphabet<>> * const 
        //     fm = dynamic_cast<const sdsl::csa_wt<sdsl::wt_int<>, 64, 64, sdsl::sa_order_sa_sampling<>, sdsl::int_vector<>, sdsl::int_alphabet<>> * const>(&fm_index);
        



        //select is sometimes one based it seems, so the first is considered the 1st, not the 0th
        sdsl::bit_vector::select_1_type b_sel(&inedgetest);
        usint mapped_l = b_sel(l+1) ;//lincoming_itr->select(l)
        usint mapped_r = b_sel(r+2)-1; //we want the subsequent node's 1 and then pack off one wt position //lincoming_itr->select(r);


//        std::cout << "mapping " << l << "," << r << " to " << mapped_l <<"," << mapped_r <<std::endl;
        hits = sdsl::restricted_unique_range_values(*wt /*fm->wavelet_tree*/, /*l, fm->wavelet_tree.size(),*/ mapped_l, mapped_r, min, max);
        if (VERBOSE >= 2) {
            std::cout <<"-- wavelet tree query in SA interval [" << l << ".."<<r << "] (mapped to " << mapped_l <<"," << mapped_r << ") "
                      << " has the following symbols within " << (max - min)/2 << " alphabet symbols of " << min + (max-min)/2 << ": " ;

            for(std::vector<long unsigned int>::iterator itr = hits.begin(); itr != hits.end(); ++itr) {
                std::cout << *itr << " ";
            }
            std::cout << std::endl;
        }

        return hits;

    }
    void CharVector::populate(usint c, DeltaVector::Encoder* encoder, usint offset)
    {

        array[c] = new CSA::DeltaVector(*encoder, offset);
//        std::cout << "populating element " << c << " with " << array.at(c)->reportSize() << " elements. Total: " << reportSize() << std::endl;
    }

    void CharVector::populate(usint c, DeltaVector *dv)
    {
        array[c] = dv;
        //      std::cout << "populating element " << c << " with " << array.at(c)->reportSize() << " elements. Total: " << reportSize() << std::endl;
    }

    void CharVector::constructF(    sdsl::int_vector<1u> &a_inedgetest)
    {
     
        inedgetest = a_inedgetest;

        }

    void CharVector::setwt(sdsl::int_vector<> &v)
    {
        wt = new sdsl::wt_int<>; //FIXME: leaking memory here
        construct_im(*wt, v);
    }

    void CharVector::syncFMIndex()
    {
//        std::cout << "Syncing FMIndex/WT to charwise bit vectors..." << std::endl;
//         sdsl::int_vector<> temp;
//         usint m = maxlength();
//         temp.resize(m);
//         for(std::map<usint, CSA::BitVector*>::iterator itr = array.begin(); itr != array.end(); ++itr)
// //        for(usint c = 1; c < 256/*FIXME:CHARS*/; c++)
//         {
//             usint c = itr->first;
//             if (c == 0) continue;
//             if(array.count(c))  {  
//                 Iterator *itr = newIterator(c);
//                 for(usint i = 0; i < array.at(c)->getSize(); ++i) {
//                     if (itr->isSet(i)) {
//                         temp[i] = c;
//                     }
//                 }
//             }
//         }
//         for (usint i = 0; i < m; ++i) {
//             if (temp[i] == 0) {
//                 std::cout << "replacing 0 element at position " << i << std::endl;
//                 temp[i] = 1; //FIXME: hack to deal with "Error: File "@24135_0" contains zero symbol."
//             }
//         }

//        sdsl::construct_im(fm_index, temp);
        //      std::cout << "WT size is " << fm_index.size() << std::endl;

    }
    void CharVector::writeTo(std::ofstream& file) const
    {
        for(std::map<usint, CSA::BitVector*>::const_iterator itr = array.begin(); itr != array.end(); ++itr)
//        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        {

            usint i = itr->first;
            if (i == 0) continue;
            if(array.count(i)) { array.at(i)->writeTo(file); }

        }

        std::cout << "Writing F to file" << std::endl;
        inedgetest.serialize(file);
        std::cout << "Writing wavelet tree to file"  << std::endl;
        wt->serialize(file);

    }
    void CharVector::load(std::ifstream &file)
    {

        wt = new sdsl::wt_int<>;
        std::cout << "Loading F from file "  << std::endl;
        inedgetest.load(file);
        std::cout << "Loading wavelet tree from file "  << std::endl;
        wt->load(file);
        std::cout << "Adding select support..." ;

        std::cout << "done.";
    }

    void CharVector::writeTo(FILE* file) const
    {
//        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        for(std::map<usint, CSA::BitVector*>::const_iterator itr = array.begin(); itr != array.end(); ++itr)
        {
            usint i = itr->first;
            if (i == 0) continue;
            if(array.count(i)) { array.at(i)->writeTo(file); }
        }

        assert(false); //assert till we find the right way to do this:
     //sdsl::store_to_file(inedgetest, file);
     //sdsl::store_to_file(wt, file);



    }
    usint CharVector::reportSize() const
    {

//        std::cout << "maxxlength = " << maxlength() << std::endl;
        usint array_size = 0;
//        for(usint i = 1; i < 256/*FIXME:CHARS*/; i++)
        for(std::map<usint, CSA::BitVector*>::const_iterator itr = array.begin(); itr != array.end(); ++itr)
        {
            usint i = itr->first;
            if (i == 0) continue;
            if(array.count(i)) { array_size += array.at(i)->reportSize(); }
        }
        return array_size;


    }
//     usint CharVector::rank(usint c, usint num, bool at_least) const {
//     // BitVector::Iterator* array_iter = this->array.at(c)->newIterator();
//     // usint r = array_iter->rank(range.first, at_least);
//     // delete array_iter;
//     // return r;
//     return 1;
// }

    usint CharVector::Iterator::select(usint i){return itr->select(i);}
    usint  CharVector::Iterator::selectNext(){return itr->selectNext();}
    usint  CharVector::Iterator::rank(usint i, bool at_least){return itr->rank(i, at_least);}
    bool  CharVector::Iterator::isSet(usint i){return itr->isSet(i);}


    CharVector::Iterator *CharVector::newIterator(usint c) const{return new Iterator(array.at(c)->newIterator());}

    CharVector::Iterator::Iterator(BitVector::Iterator *itr) : itr(itr)
    {
    }
    CharVector::Iterator::~Iterator() 
    {
        delete itr;
    }
        
} //namespace CSA
