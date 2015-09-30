#!/usr/bin/python

import sys
from math import ceil,sqrt

if len(sys.argv) < 4:
    print("Usage: ", sys.argv[0], "<input.valuev> <STDDEVs per bin> <STDDEV per kb> <output.valuev>")
    sys.exit(1)

sigma = float(sys.argv[3])
m = float(sys.argv[2]) / 2  # divide by two just because I solved the algebra in terms of a delta defined as the bin boundary to the center

# U(n) = 1/4 m^2 sigma^2 (4 (c_1+n)^2-1)  (c_1 is an arbitrary parameter)
# let c_1 = 0
# position = 1/4 * m**2 * sigma**2 * (4 * n**2 - 1)
# formulate a function that gives a real number for bin from a position, inverse of U(n)


# p = 1/4 * m**2 * sigma**2 * (4 * n**2 - 1)
# p * 4 = m**2 * sigma**2 * (4 * n**2 - 1)
# p * 4 / m**2 = sigma**2 * (4 * n**2 - 1)
# p * 4 / m**2 / sigma**2 = 4 * n**2 - 1
# (p * 4 / m**2 / sigma**2) + 1 = 4 * n**2
# ((p * 4 / m**2 / sigma**2) + 1) / 4 =  n**2
# sqrt(((p * 4 / m**2 / sigma**2) + 1) / 4), 2) = n

c_1 = 0

def bin_pos(n):
    "Return the bin position"
    return     1/4 * m**2 * sigma**2 * (4 * (c_1 + n)**2 - 1)

def inv_bin_pos(p):
    "Return the inverse bin position"
    return sqrt(((p * 4 / m**2 / sigma**2) + 1) / 4) - c_1

#    return log((((p * 4 / m**2 / sigma**2) + 1) / 4), 2)
#    return log( ((4* p / m**2 / sigma**2) + 1 ) / 4 ,2)

def quantize(n):
    estimated = inv_bin_pos(n)
    bin = ceil(estimated)
    upper = bin_pos(bin)
    lower = bin_pos(bin-1)
    # return the mid point
    mid_point = lower + (upper-lower)/2
    if lower > n or upper < n:
        print(n,estimated, bin,lower,mid_point,upper)
    return mid_point

output = open(sys.argv[4], "w")
for lno,line in enumerate(open(sys.argv[1])):

    if lno % 3 == 0 or lno % 3 == 2:
        output.write(line)
        continue
    # otherwise, is a second line of three line sets
    fields = line.strip().split("\t")
    newfields = fields[:2] + [str(quantize(float(field))) for field in fields[2:]]
    output.write("\t" + "\t".join(newfields) + "\n")
    
    
output.close()

# for i in range(1,100):
#     print(quantize(i))

# print("xxxxxxxxxxxxxxxxxxxxxxxxxxx")
# for i in range(1,100):
#     print(i, inv_bin_pos(i), bin_pos(inv_bin_pos(i)))
