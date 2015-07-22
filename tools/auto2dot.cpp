// auto2dot.cpp <gcsa format graph>
// generates a graphviz dot text file description of a graph


#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include <list>

class Edge;
class Node {
public:
    Node(unsigned int _label, unsigned int _value, unsigned int _pos) { 
        label = _label; 
        value = _value; 
        


        pos = _pos;}
    std::vector<Edge *> in_edges;
    std::vector<Edge *> out_edges;
    unsigned int label;
    unsigned int value;
    unsigned int pos;
};

class Edge {
public:
    Edge(Node *_from, Node *_to) { from = _from; to = _to;}
    Node *from;
    Node *to;
};


std::vector<Node *> nodes;
std::vector<Edge *> edges;
int main(int argc, char** argv)
{
    int fd = 0;

    fd = open(argv[1], O_RDONLY);
    unsigned int num_nodes = 0;
    unsigned int num_edges = 0;
    read(fd, &num_nodes, 4);
    read(fd, &num_edges, 4);

    unsigned int min = (unsigned int)-1;
    unsigned int max = 0;
    Node * min_node;
    Node * max_node;

    for (int i = 0; i < num_nodes; ++i) {
        unsigned int label;
        unsigned int value;
        read(fd, &label, 4);
        read(fd, &value, 4);
        Node * n = new Node(label, value, i);
        if (label > max) {
            max = label;
            max_node = n;
        }
        if (label < min) {
            min = label;
            min_node = n;
        }
        nodes.push_back(n);
    }
    std::cout << "// Nodes: " << num_nodes << " Edges: " << num_edges << " Max label: " << max << " Min label: " << min << std::endl;
    for (int i = 0; i < num_edges; ++i) {    
        unsigned int from;
        unsigned int to;
        read(fd, &from, 4);
        read(fd, &to, 4);
        Edge * e = new Edge(nodes[from], nodes[to]);
        edges.push_back(e);
        nodes[to]->in_edges.push_back(e);
        nodes[from]->out_edges.push_back(e);
    }

    // loading should be done here.

    std::cout << "digraph G {" << std::endl;
    for (std::vector<Node *>::iterator ni = nodes.begin(); ni != nodes.end(); ++ni) {
        std::cout << "    n" << (*ni)->pos << " [label=\"" << (*ni)->label << "\"];" << std::endl;
    }

    for (std::vector<Edge *>::iterator ei = edges.begin(); ei != edges.end(); ++ei) {
        std::cout << "    n" <<  (*ei)->from->pos << " -> n" << (*ei)->to->pos << ";" << std::endl;
    }
    std::cout << "}" << std::endl;
}
