#!/bin/sh -xe
rm -f fig1_base.automaton  fig1_cleaned	 fig1_base.gcsa	     fig1_cleaned.automaton fig1_base.backbone 
/bin/time -v ../clean_alignment fig1.ex fig1_cleaned
/bin/time -v ../build_automaton  fig1_cleaned 0
/bin/time -v ../determinize  fig1_cleaned.automaton fig1_base
/bin/time -v ../build_index  fig1_base
/bin/time -v ../gcsa_test fig1_base fig1.patterns



