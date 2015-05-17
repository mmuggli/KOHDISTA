#!/usr/bin/python
import sys

import struct
item = struct.Struct("i")

target = open(sys.argv[2], "wb")
targetrmaps = open(sys.argv[2] + ".frag2rmap", "w")
query = open(sys.argv[3], "wb")
pos = 0
rmap = 0
for i,line in enumerate(open(sys.argv[1])):
    fields = line.strip().split("\t")
    if "map" in line: name = line.strip()
    # write the target binary file
    if len(fields) > 2 and i % 3 == 1:
        frags = fields[2:]
        for frag in frags:
            target.write(item.pack(int(float(frag)*1000)))
        target.write(item.pack(int(1000000000)))
        targetrmaps.write(str(rmap) + "\t" + str(pos) + "\t" + name + "\n")
        rmap += 1
        pos += len(frags) + 1
    # write the query file
    if len(fields) > 3 and i % 3 == 1:
        query.write(item.pack(len(fields) - 3))
        for frag in fields[3:]:
            query.write(item.pack(int(float(frag)*1000)))

print("Processed aprox",i,"lines")               

target.close()
query.write(item.pack(0))
query.close()


