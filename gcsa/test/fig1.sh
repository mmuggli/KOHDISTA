#!/bin/sh -xe
rm -f fig1_base.automaton  fig1_cleaned	 fig1_base.gcsa	     fig1_cleaned.automaton fig1_base.backbone 
../clean_alignment fig1.ex fig1_cleaned
../build_automaton -b fig1_cleaned 0
../determinize -b fig1_cleaned.automaton fig1_base
../build_index -b fig1_base
../gcsa_test fig1_base fig1.patterns



