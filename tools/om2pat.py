#!/usr/bin/python

import sys,re, struct
item = struct.Struct("i")

def load_optmap(fname):
    """Return a list of the frags for a binary opt map file fname"""
    f = open(fname, "rb")
    frags = []
    word = f.read(4)
    while len(word) == 4:
        frag = item.unpack(word)
        frags.append(frag[0])
        word = f.read(4)
    f.close()
    return frags


opt_map = load_optmap(sys.argv[1])

binfile = open(sys.argv[2], "wb")

print("len(opt_map):", len(opt_map))
num_pats = 5
pat_len = 20
mult = 50

def emit(frag):
    binfile.write(item.pack(frag))
    print(frag,' ', end='')

for i in range(num_pats):
    # each pattern starts with the number of elements
    binfile.write(item.pack(pat_len)) 
    # followed by the elements
    print("pat = opt_map[",i*mult,":",i*mult+pat_len,"] = ", end='')
    for fnum, frag in enumerate(opt_map[i*mult:i*mult+pat_len]):
        if fnum == 1:
            print("(",frag,opt_map[i*mult+fnum + 1],frag + opt_map[i*mult+fnum + 1],")",end='')
            emit(frag + opt_map[i*mult+fnum + 1])
            continue
        if fnum == 2:
            continue
        if fnum == 5:
            emit(int(frag*.4))
            emit(int(frag*.6))
            continue
        emit(frag + 5)

    print()
# 0 elements means we're done reading
binfile.write(item.pack(0))

binfile.close()

print (list(enumerate(opt_map)))
