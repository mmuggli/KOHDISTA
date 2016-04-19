#!/usr/bin/python

import os
import re
import sys
import optparse
import tempfile
import subprocess
import shutil
p = optparse.OptionParser()

p.add_option("--query", action="store", dest="query", help="query for incrementally longer matches from the right end of these Rmaps, up to 'min-overlap' value.  Repeats this for the reverse of each of these Rmaps.")
p.add_option("--target", action="store", dest="target", help="build an atomaton for the concatenation of these Rmaps")
p.add_option("--bin-size", action="store", dest="bin_size", help="round input fragments to this bin size for quantization")
p.add_option("--min-overlap", action="store", dest="min_overlap", help="minimum number of compound fragments for an alignment")
p.add_option("--detailed", action="store_true", dest="detailed", help="report a Valouev et al. style list of fragment groupings for the alignment and their accompanying s-scores and t-score")
p.add_option("--keep-tempdir", action="store_true", dest="keep_tempdir", help="don't delete the temporary directory (useful for debuggin)")
p.add_option("--chi-squared-cdf-thresh", action="store", dest="chi2cdf_thresh", default=".1", help="Chi^2 CDF threshold used for pruning the partial alignments which are extended to full alignments (default=%default)")
p.add_option("--max-desorption-thresh", action="store", dest="max_desorption_thresh", help="value below which skip edges are added to the target automaton")
p.add_option("--min-desorption-thresh", action="store", dest="min_desorption_thresh", help="value below which fragments are discarded from the query")
p.add_option("--min-t-score", action="store", dest="min_t_score", help="minimum t-score required to report an alignment")
p.add_option("--sigma", action="store", dest="sigma_kbp", help="per Kbp standard deviation of fragment size estimation error (default=%default)", default=".58")
p.add_option("--constant-sigma", action="store_true", dest="bounded_sigma", help="Assume all fragments have sigma as the upper bound of their stddev, instead of calculated from fragment length.")
p.add_option("--trim-query-ends", action="store_true", dest="trim_query_ends", help="Remove one fragment off each end of query rmaps, for fragments that are only cleaved by an enzyme on one end")
p.add_option("--query-order", action="store", dest="query_order", help="How many combinations of successive missed sites in target [0..n) to accomodate")
p.add_option("--one-sided-error", action="store_true", dest="one_sided_error", help="Assume only one sequence has a sizing error (e.g. in silico digested contigs aligned to whole genome optical map)")
p.add_option("--single_ended", action="store_true", dest="single_ended", help="Only search from one end of a query.  Suitable for 'fit' alignments but not for 'overlap' alignments")
p.add_option("--pretend", action="store_true", dest="pretend", help="Don't run anything, just show what would be done")
p.set_defaults(pretend=False)
p.set_defaults(single_ended=False)
p.set_defaults(chi2cdf_thresh=".1")
p.set_defaults(one_sided_error=False)
p.set_defaults(detailed=False)
p.set_defaults(bounded_sigma=False)
p.set_defaults(keep_tempdir=False)
p.set_defaults(trim_query_ends=False)
p.set_defaults(bin_size="100")
p.set_defaults(min_t_score="1")
p.set_defaults(max_desorption_thresh="1000")
p.set_defaults(min_desorption_thresh="500")
p.set_defaults(min_overlap="10")
p.set_defaults(query_order="2")

opts,args = p.parse_args()
#query = opts.query
print("pretend is", opts.pretend)


# install directory
dopp_instdir = os.path.dirname(os.path.realpath(sys.argv[0]))

tempdir = tempfile.mkdtemp()
print("# Storing temporary files in directory", tempdir)

# python /s/chopin/l/grad/muggli/git/rlcsa/tools/valouev2bin.py  ../ecoli_verif_100x_experimental.valuev sim_ecoli_XhoI_rev.bin sim_ecoli_XhoI_rev_pat.bin
print("# *** build the automaton for the target")
cmd = ["/usr/bin/python", dopp_instdir +  "/tools/valouev2bin.py", opts.target, tempdir + "/target.bin", tempdir + "/target_pat.bin", str(float(opts.min_desorption_thresh)/1000.0)]
print("# ***   step 1.) convert target rmaps file to binary file with command:")
print(" ".join(cmd))
if not opts.pretend:
    p = subprocess.Popen(cmd)
    sout, serr = p.communicate()
    ret = p.returncode
    print((sout), (serr), ret)
 
# ~/git/rlcsa/tools/om2automaton sim_ecoli_XhoI_rev.bin sim_ecoli_XhoI_rev.automaton 100 0
cmd = [dopp_instdir + "/tools/om2automaton", tempdir + "/target.bin",  tempdir + "/target.automaton",
                      opts.bin_size, # quantization bin size
                      "0",# prefix of the whole set to use (for testing with smaller substring of a large set of rmaps concatenated
                      opts.max_desorption_thresh] 

print("# ***   step 2.) build the automaton from the binary file with command:")
print(" ".join(cmd))                     
if not opts.pretend:
    p = subprocess.Popen(cmd)
    sout, serr = p.communicate()
    ret = p.returncode
    print((sout), (serr), ret)


# cmd = ["/usr/bin/python", dopp_instdir +  "/tools/valouev2bin.py", opts.query, tempdir + "/query.bin", tempdir + "/query_pat.bin", str(float(opts.min_desorption_thresh)/1000.0)]
# print("*** convert query rmaps file to binary file (for the side effect of producing the patterns file with command:", " ".join(cmd))

# p = subprocess.Popen(cmd)
# sout, serr = p.communicate()
# ret = p.returncode
# print((sout), (serr), ret)

# guard against something Omfio parser class can't handle
id_re = re.compile("^\W*(\w+)")
if not opts.pretend:
    for lno, line in enumerate(open(opts.query)):
        if lno % 3 == 0:
            if ord("0") >= ord(id_re.search(line).group(1)[0]) >= ord("9"):
                print("Error: ", opts.query, "line", lno, "query rmap id's must not begin with a digit")
                sys.exit(1)


# ~/git/rlcsa/gcsa/determinize -b sim_ecoli_XhoI_rev.automaton sim_ecoli_XhoI_rev_base
cmd = [dopp_instdir + "/gcsa/determinize", "-b", tempdir + "/target.automaton", tempdir + "/target_base"]
print("# *** Determinizing the automaton with command:")
print(" ".join(cmd))                     
if not opts.pretend:
    p = subprocess.Popen(cmd)
    sout, serr = p.communicate()
    ret = p.returncode
    print((sout), (serr), ret)

# ~/git/rlcsa/gcsa/build_index -b  sim_ecoli_XhoI_rev_base
cmd = [dopp_instdir + "/gcsa/build_index", "-b", tempdir + "/target_base"]
print("# *** Building the GCSA data structure with command:")
print( " ".join(cmd))                     
if not opts.pretend:
    p = subprocess.Popen(cmd)
    sout, serr = p.communicate()
    ret = p.returncode
    print((sout), (serr), ret)


print("# *** Moving rmap index file returned")
print("mv", tempdir + "/target.bin.frag2rmap", tempdir + "/target_base.frag2rmap")
if not opts.pretend:
    shutil.move(tempdir + "/target.bin.frag2rmap", tempdir + "/target_base.frag2rmap")



# cp sim_ecoli_XhoI_rev.bin.frag2rmap sim_ecoli_XhoI_rev_base.frag2rmap #FIXME this step is convoluted
# echo "arguments" "$@"
# /bin/time -v ~/git/rlcsa/gcsa/gcsa_test sim_ecoli_XhoI_rev_base sim_ecoli_XhoI_rev_pat.bin  -b -l "$@"
cmd = [dopp_instdir + "/gcsa/gcsa_test", "-b", "-l", tempdir + "/target_base", opts.query]
if opts.detailed:
    cmd.append("-d")
if opts.bounded_sigma:
    cmd.append("-B")
if opts.trim_query_ends:
    cmd.append("-t")
if opts.one_sided_error:
    cmd.append("-1")
if opts.single_ended:
    cmd.append("-M")
cmd.append("-Q" + opts.query_order)    
cmd.append("-O" + opts.min_overlap)
cmd.append("-C" + opts.chi2cdf_thresh)
cmd.append("-T" + opts.min_t_score)
cmd.append("-Z" + opts.sigma_kbp)
print("# *** Performing alignment with command:")
print(" ".join(cmd))                     
if not opts.pretend:
    p = subprocess.Popen(cmd)
    sout, serr = p.communicate()
    ret = p.returncode
    print((sout), (serr), ret)

if not opts.keep_tempdir:
    print("# *** Removing temporary directory", tempdir)
    shutil.rmtree(tempdir)
