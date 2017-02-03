# Takes a file in the format used by Valouev and extracts k-mers as new maps
import sys
import numpy as np
k = 12 # window size                                                                                                                                                                                                                                                           
for lno,line in enumerate(open(sys.argv[1])):
    if lno % 3 == 0:
        label = line.strip()
    if lno % 3 == 1:
        parts = line.strip().split()
        enz, enzabr = parts[:2]
        frags =  parts[2:]
        windows = len(frags) - k
        if windows >= 0:
            for window in range(windows+1):
                print(label + "-" + str(window))
                print("\t".join(parts[:2] + frags[window:window+k]))
                print()
