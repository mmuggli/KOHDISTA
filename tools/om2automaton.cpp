// om2automaton <binary optical map> <gcsa format graph file>



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
#include <fstream>
#include <sys/stat.h>
#include <stdlib.h>

class Node {
public:
    Node(unsigned int _label, unsigned int _value, unsigned int _pos) { 

        label = _label; 
        value = _value; 
        
        if (_label == 0) {
            int ijk = 0;
        }



        pos = _pos;}
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


unsigned int remap(unsigned int l)
{
    return l;
    if (l == 0) return 0;
    if (l == (unsigned int)-1) return 255;
    //FIXME: add log stuff in here
    float max = 1364039;//(unsigned int)-1;
    float thelog = log(l);
    float ret =  thelog * 253 / log(max)  + 1;

    return ret;

}

unsigned int mytoupper(unsigned int lab)
{
        if (lab == 0) return lab;
        if (lab & 0x1) lab -= 1; // mark it as a backbone "uppercase" node
        return lab;
}

unsigned int mytolower(unsigned int lab)
{
        lab |= 0x1; // mark it as a nonbackbone "lowercase" node        
        return lab;
}

bool myisupper(unsigned int lab)
{
    if (lab & 0x1) return false;
    return true;
}

bool mislower(unsigned int lab)
{
    if (lab & 0x1) return true;
    return false;
}

int bin_size = 1;
const int DESORPTION_THRESH = 1000;

unsigned int quantize(unsigned int val)
{
    unsigned int new_val = 0;
    if (val % bin_size < bin_size / 2.0)
        new_val = val - val % bin_size;
    else
        new_val = val - val % bin_size + bin_size;
    return new_val;
}


int main(int argc, char** argv)
{


    if (argc < 5) {
        printf("Usage: %s <binary optical map> <gcsa format graph> <quantization bin size> <file prefix>\n", argv[0]);
        exit(1);
    }
    bin_size = atol(argv[3]);
    unsigned long long requested_elems = atoll(argv[4]);
    printf("Quantizing with bin size %d\n", bin_size);
    std::map<unsigned int, unsigned int> counts;

    

    int ifd = 0;
    char *addr = 0;
    char *ifname = argv[1];
    struct stat stat_struct;
    fstat(ifd, &stat_struct);
    int file_size = stat_struct.st_size;
    unsigned int elems = 0;
    if (requested_elems == 0) {
        elems = file_size / 4; //364876384 / 4;
    }else{
        elems = requested_elems;
    }
    ifd = open(ifname, O_RDONLY);
//printf("Attempting to mmap %d elements from %s\n", (elems - 1) * 4, fname);
    // // void *mmap(void *addr, size_t length, int prot, int flags,
    // //            int ifd, off_t offset);
    addr = (char*)malloc(elems * 4);
    off_t  pa_offset = 0;
    
    printf("Reading file %s\n", ifname);
    int r = read(ifd, addr, elems * 4);
    printf("read %d bytes\n", r);
    if (r != 4*elems) { printf("Incomplete read ERROR!\n");}

    int j = 0;
    Node *start = new Node(-2, 0, 0);
    ++j;
    Node *current = start;
    std::vector<Node *> nodes;
    std::vector<Edge *> edges;
    nodes.push_back(start);

    printf("building baseline\n");
    int i;

    for (i = 0; i < r / 4; ++i) {

        unsigned int lab = mytoupper(quantize(((unsigned int *)addr)[i]));

        Node *n = new Node(lab, j, j);
        ++j;
        nodes.push_back(n);
        Edge *e = new Edge(current, n);
        edges.push_back(e);
        current = n;
    }


    Node *end = new Node(0, j, j);
    int end_value = j;
    std::cout <<"end node value is " <<end_value << std::endl;
    j++;

    nodes.push_back(end);
    std::cout << "creating edge from " << current->pos << " to " << end->pos << std::endl;
    Edge *e = new Edge(current, end);
    edges.push_back(e);

    std::vector<Node *> skip_nodes;
    std::vector<Edge *> skip_edges;


    std::cout << "Adding skip edges and nodes" << std::endl;;
    for (i = 0; i < r / 4 - 1; ++i) {

        // add order 1 skip nodes
        unsigned int lab = mytolower(nodes[i+1]->label + nodes[i+2]->label);

        Node *sumnode = new Node(lab, nodes[i+1]->value + end_value, j); //j+1 means non-backbone none and no alignment which we can hopefully detect somehow 
        //std::cout << "o1 skip node lab=" << lab << " val_and_pos=" << j << std::endl;
        ++j;
        nodes.push_back(sumnode);
        Edge *e1 = new Edge(nodes[i], sumnode);
        edges.push_back(e1);
        Edge *e2 = new Edge(sumnode, nodes[i+3]);
        edges.push_back(e2);


        //add order 2 skip nodes
        if (i  < r / 4 - 4) {
            unsigned int lab = mytolower(nodes[i+1]->label + nodes[i+2]->label + nodes[i+3]->label);

            Node *sumnode = new Node(lab, nodes[i+1]->value + end_value, j); //j+1 means non-backbone none and no alignment
            ++j;
            nodes.push_back(sumnode);
            Edge *e1 = new Edge(nodes[i], sumnode);
            edges.push_back(e1);
            Edge *e2 = new Edge(sumnode, nodes[i+4]);
            edges.push_back(e2);
        }

        // add skip edges
        if (nodes[i+1]->value < DESORPTION_THRESH) {
            Edge *eskip = new Edge(nodes[i], nodes[i+2]);
            edges.push_back(eskip);
        }

        // add double skip edges
        if (nodes[i+1]->value < DESORPTION_THRESH && nodes[i+2]->value < DESORPTION_THRESH) {
            Edge *eskip = new Edge(nodes[i], nodes[i+3]);
            edges.push_back(eskip);
        }


    }
    char *ofname = argv[2];
    printf("writing file %s\n", ofname);
    //int ofd = open(ofname, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    std::ofstream ofd(ofname, std::ofstream::binary);
    //if (ofd == -1) printf("problem opening file to write ERROR!\n");
    unsigned int numnodes = nodes.size();
    unsigned int numedges = edges.size();
    ofd.write((char*)&numnodes, 4);
    ofd.write((char*)&numedges, 4);
    for (std::vector<Node *>::iterator ni = nodes.begin(); ni != nodes.end(); ++ni) {
        unsigned int lab = ((*ni)->label);
        if (lab == 0) std::cout << "0 lab at node " << (*ni)->value <<  " pos " << (*ni)->pos << std::endl;

        // quantize, preserving backboniness
        if (myisupper(lab)) {
            lab = mytoupper((lab));
        } else {
            lab = mytolower((lab));
        }
//        lab = remap(lab);
        
        counts[lab] += 1;
        ofd.write((char*)&lab, 4);
        ofd.write((char*)&((*ni)->value), 4);
    }

    std::cout << "Number of nodes: " << numnodes << std::endl;
    std::cout << "largest node pos: " << (*(nodes.end() - 1))->pos << std::endl;
    for (std::vector<Edge *>::iterator ei = edges.begin(); ei != edges.end(); ++ei) {
        ofd.write((char*)&((*ei)->from->pos), 4);
        ofd.write((char*)&((*ei)->to->pos), 4);
    }

    printf("size of active alphabet: %d\n", counts.size());

     // for(std::map<unsigned int, unsigned int>::iterator ci = counts.begin(); ci != counts.end(); ++ci) {
     //     std::cout << "counts[" << ci->first << "] = " << ci->second << std::endl;
     // }

    ofd.close();
}
