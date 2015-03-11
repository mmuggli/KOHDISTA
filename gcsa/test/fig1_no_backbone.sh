#!/bin/sh -xe
rm -f fig1_base.automaton  fig1_cleaned	 fig1_base.gcsa	     fig1_cleaned.automaton fig1_base.backbone 
../clean_alignment fig1.ex fig1_cleaned
../build_automaton  fig1_cleaned 0
../determinize  fig1_cleaned.automaton fig1_base
../build_index  fig1_base
../gcsa_test fig1_base fig1.patterns



