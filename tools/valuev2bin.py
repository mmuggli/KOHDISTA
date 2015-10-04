#!/usr/bin/python
import sys

import struct

if len(sys.argv) < 4:
    print("Usage: valuev2bin.py map.valuev map.bin pat.bin")
    sys.exit(1)


def write_files(frags, rmap, pos, name):
    for frag in frags:
        fragsize = int(float(frag)*1000)
        target.write(item.pack(fragsize))
        targettext.write(str(fragsize) + "\n")
    fragsize = int(1000000000)
    target.write(item.pack(fragsize))
    targettext.write(str(fragsize) + "\n")
    targetrmaps.write(str(rmap) + "\t" + str(pos) + "\t" + name + "\n")
    

    
item = struct.Struct("i")

target = open(sys.argv[2], "wb")
targetrmaps = open(sys.argv[2] + ".frag2rmap", "w")
query = open(sys.argv[3], "wb")
targettext = open(sys.argv[2] + ".text", "w")
pos = 0
rmap = 0
for i,line in enumerate(open(sys.argv[1])):
    if i % 100000 == 0: print("processing lineno:",i)
    # parse name line
    if i % 3 == 0:
        name = line.strip()
        
    # parse the frag line and write the target binary file
    if i % 3 == 1:
        fields = line.strip().split("\t")
        frags = fields[2:]
        write_files(frags, rmap, pos, name)
        rmap += 1
        pos += len(frags) + 1
        frags.reverse()
        write_files(frags, rmap, pos, name + "_reversed")
        rmap += 1
        pos += len(frags) + 1
        
        # write the query file
        query.write(item.pack(len(frags)))
        #print("writting",len(fields) - 2,"frags: ",end=" ")
        for frag in frags:
            #print(int(float(frag)*1000), end=" ")
            query.write(item.pack(int(float(frag)*1000)))
        #print()
print("Processed aprox",i/3,"rmaps")               

target.close()
query.write(item.pack(0))
query.close()


