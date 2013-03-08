#!/usr/bin/env python
# encoding: utf-8
"""
From read names, create metadata table.
"""
import sys
import itertools
from toolshed import nopen
from parsers import read_fasta

def main(args):
    # fields from issake
    fields = "contig_id length reads_needed coverage seed_name v_region j_region".split()
    # the only fields i believe make any sense to keep
    out = "contig_id v_region j_region length coverage".split()
    print "#" + "\t".join(out)
    with nopen(args.fasta) as fasta:
        for name, seq in read_fasta(fasta):
            # remove some text from iSSAKE output
            name = name.replace("size","").replace("cov","").replace("read","").replace("seed:","")
            d = dict(zip(fields, name.split("|")))
            print "\t".join([d[o] for o in out])

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("in", help="fasta with detailed read names")
    p.add_argument("out", help="fasta with unique, truncated read names")
    p.add_argument("meta", help="metadata for use in topiary explorer")
    main(p.parse_args())