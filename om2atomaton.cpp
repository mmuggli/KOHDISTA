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

int main(int argc, char** argv)
{
    std::map<unsigned int, unsigned int> counts;
    unsigned int elems = 364876384 / 4;
    int fd = 0;
    char *addr = 0;
    fd = open(argv[1], O_RDONLY);
//printf("Attempting to mmap %d elements from %s\n", (elems - 1) * 4, fname);
    // // void *mmap(void *addr, size_t length, int prot, int flags,
    // //            int fd, off_t offset);
    addr = (char*)malloc(elems * 4);
    off_t  pa_offset = 0;
    printf("Reading file\n");
    int r = read(fd, addr, elems * 4);
    printf("read %d bytes\n", r);
    if (r != 4*elems) { printf("Incomplete read ERROR!\n");}

    int j = 0;
    Node *start = new Node(-1, 0, 0);
    ++j;
    Node *current = start;
    std::vector<Node *> nodes;
    std::vector<Edge *> edges;
    nodes.push_back(start);

    printf("building baseline\n");
    int i;

    for (i = 0; i < r / 4; ++i) {

        unsigned int val = ((unsigned int *)addr)[i];
        if (val == 0) continue;
        Node *n = new Node(val, j, j);
        ++j;
        nodes.push_back(n);
        Edge *e = new Edge(current, n);
        edges.push_back(e);
        current = n;
    }


    Node *end = new Node(0, 0, j);
    j++;
    nodes.push_back(end);
    std::cout << "creating edge from " << current->pos << " to " << end->pos << std::endl;
    Edge *e = new Edge(current, end);
    edges.push_back(e);

    std::vector<Node *> skip_nodes;
    std::vector<Edge *> skip_edges;


    std::cout << "Adding skip edges" << std::endl;;
    for (i = 0; i < r / 4 - 1; ++i) {
        Node *sumnode = new Node(nodes[i+1]->label + nodes[i+2]->label, nodes[i+1]->value, j);
        ++j;
        nodes.push_back(sumnode);
        Edge *e1 = new Edge(nodes[i], sumnode);
        edges.push_back(e1);
        Edge *e2 = new Edge(sumnode, nodes[i+3]);
        edges.push_back(e2);
    }
    printf("writing file\n");
    int fd2 = open(argv[2], O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
    if (fd2 == -1) printf("problem opening file to write ERROR!\n");
    unsigned int numnodes = nodes.size();
    unsigned int numedges = edges.size();
    write(fd2, &numnodes, 4);
    write(fd2, &numedges, 4);
    for (std::vector<Node *>::iterator ni = nodes.begin(); ni != nodes.end(); ++ni) {
        unsigned int lab = ((*ni)->label);
        if (lab == 0) std::cout << "0 lab at node " << (*ni)->value <<  " pos " << (*ni)->pos << std::endl;
        lab = remap(lab);
        
        counts[lab] += 1;
        write(fd2, &lab, 4);
        write(fd2, &((*ni)->value), 4);
    }

    std::cout << "Number of nodes: " << numnodes << std::endl;
    std::cout << "largest node pos: " << (*(nodes.end() - 1))->pos << std::endl;
    for (std::vector<Edge *>::iterator ei = edges.begin(); ei != edges.end(); ++ei) {
        write(fd2, &((*ei)->from->pos), 4);
        write(fd2, &((*ei)->to->pos), 4);
    }

    printf("size of active alphabet: %d\n", counts.size());

    // for(std::map<unsigned int, unsigned int>::iterator ci = counts.begin(); ci != counts.end(); ++ci) {
    //     std::cout << "counts[" << ci->first << "] = " << ci->second << std::endl;
    // }

    close(fd2);
}
