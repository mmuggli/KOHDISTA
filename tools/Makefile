

all: om2automaton auto2dot

om2automaton: om2automaton.cpp
	g++ -O3  om2automaton.cpp -o om2automaton

auto2dot: auto2dot.cpp
	g++ auto2dot.cpp -o auto2dot -g

clean:
	rm -f om2automaton auto2dot
