#!/usr/bin/python

import os
import sys
import optparse
import tempfile
import subprocess
import shutil
p = optparse.OptionParser()

p.add_option("--query", action="store", dest="query")
p.add_option("--target", action="store", dest="target")
p.add_option("--detailed", action="store_true", dest="detailed")
p.set_defaults(detailed=False)

opts,args = p.parse_args()
#query = opts.query



# install directory
dopp_instdir = os.path.dirname(os.path.realpath(sys.argv[0]))

tempdir = tempfile.mkdtemp()
print("Storing temporary files in directory", tempdir)

# python /s/chopin/l/grad/muggli/git/rlcsa/tools/valouev2bin.py  ../ecoli_verif_100x_experimental.valuev sim_ecoli_XhoI_rev.bin sim_ecoli_XhoI_rev_pat.bin
print("*** build the automaton for the target")
print("***   step 1.) convert target rmaps file to binary file")
p = subprocess.Popen(["/usr/bin/python", dopp_instdir +  "/tools/valouev2bin.py", opts.target, tempdir + "/target.bin", tempdir + "/target_pat.bin"])
sout, serr = p.communicate()
ret = p.returncode
print((sout), (serr), ret)
 
# ~/git/rlcsa/tools/om2automaton sim_ecoli_XhoI_rev.bin sim_ecoli_XhoI_rev.automaton 100 0
print("***   step 2.) build the automaton from the binary file")                     
p = subprocess.Popen([dopp_instdir + "/tools/om2automaton", tempdir + "/target.bin",  tempdir + "/target.automaton",
                      "100", # quantization bin size
                      "0" # prefix of the whole set to use (for testing with smaller substring of a large set of rmaps concatenated
])
sout, serr = p.communicate()
ret = p.returncode
print((sout), (serr), ret)


print("*** build the automaton for the target")
print("***   convert query rmaps file to binary file (for the side effect of producing the patterns file")
p = subprocess.Popen(["/usr/bin/python", dopp_instdir +  "/tools/valouev2bin.py", opts.query, tempdir + "/query.bin", tempdir + "/query_pat.bin"])
sout, serr = p.communicate()
ret = p.returncode
print((sout), (serr), ret)


# ~/git/rlcsa/gcsa/determinize -b sim_ecoli_XhoI_rev.automaton sim_ecoli_XhoI_rev_base
print("*** Determinizing the automaton")                     
p = subprocess.Popen([dopp_instdir + "/gcsa/determinize", "-b", tempdir + "/target.automaton", tempdir + "/target_base"])
sout, serr = p.communicate()
ret = p.returncode
print((sout), (serr), ret)

# ~/git/rlcsa/gcsa/build_index -b  sim_ecoli_XhoI_rev_base
print("*** Building the GCSA data structure")                     
p = subprocess.Popen([dopp_instdir + "/gcsa/build_index", "-b", tempdir + "/target_base"])
sout, serr = p.communicate()
ret = p.returncode
print((sout), (serr), ret)


print("*** Moving rmap index file returned", shutil.move(tempdir + "/target.bin.frag2rmap", tempdir + "/target_base.frag2rmap"))



# cp sim_ecoli_XhoI_rev.bin.frag2rmap sim_ecoli_XhoI_rev_base.frag2rmap #FIXME this step is convoluted
# echo "arguments" "$@"
# /bin/time -v ~/git/rlcsa/gcsa/gcsa_test sim_ecoli_XhoI_rev_base sim_ecoli_XhoI_rev_pat.bin  -b -l "$@"
print("*** Performing alignment")                     
p = subprocess.Popen([dopp_instdir + "/gcsa/gcsa_test", "-b", "-l", tempdir + "/target_base", tempdir + "/query_pat.bin"])
sout, serr = p.communicate()
ret = p.returncode
print((sout), (serr), ret)

print("*** Removing temporary directory", tempdir)
shutil.rmtree(tempdir)
