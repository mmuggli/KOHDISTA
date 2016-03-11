#!/usr/bin/python
import sys

import struct

if len(sys.argv) < 5:
    print("Usage: valouev2bin.py map.valuev map.bin pat.bin min_desorption_thresh")
    sys.exit(1)


def write_files(frags, rmap, pos, name):

    for frag in frags:

        fragsize = int(float(frag)*1000)
#        if fragsize < 50 : print("found 0 frag while converting", (frags, rmap, pos, name))
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
min_desorption_thresh = float(sys.argv[4])
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
        fields = line.strip().split()
        frags = fields[2:]
#        print("found",len(frags),"frags in line")
        write_files(frags, rmap, pos, name)
        rmap += 1
        pos += len(frags) + 1
        frags.reverse()
        write_files(frags, rmap, pos, name + "_reversed")
        rmap += 1
        pos += len(frags) + 1
        
        # write the query file

        qfrags = [frag for frag in frags if float(frag) > min_desorption_thresh]
        query.write(item.pack(len(qfrags)))
#        print("writting",len(qfrags) ,"frags: ",end=" ")
        for qfrag in qfrags:
#            print(int(float(qfrag)*1000))
            query.write(item.pack(int(float(qfrag)*1000)))
#        print()
print("Processed aprox",i/3,"rmaps")               

target.close()
query.write(item.pack(0))
query.close()


