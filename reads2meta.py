#!/usr/bin/env python
# encoding: utf-8
"""
From read names, create metadata table.
"""
import sys
import itertools
from toolshed import nopen
from parsers import read_fasta, write_fasta

def main(args):
    # fields from issake
    fields = "contig_id length reads avg_coverage seed v_region j_region".split()
    # the only fields i believe make any sense to keep
    out = "id v_region j_region length reads avg_coverage percent_of_total".split()
    # total reads used in assembly
    total = 0.
    with nopen(args.fasta_in) as fasta:
        for name, seq in read_fasta(fasta):
            name = name.replace("size","").replace("cov","").replace("read","").replace("seed:","")
            d = dict(zip(fields, name.split("|")))
            total += int(d['reads'])
    with nopen(args.fasta_in) as fasta, open(args.fasta_out, 'wb') as fasta_out, \
            open(args.meta, 'wb') as meta:
        # print header
        meta.write("\t".join(out) + "\n")
        for i, retvals in enumerate(read_fasta(fasta)):
            name, seq = retvals
            # remove some text from iSSAKE output
            name = name.replace("size","").replace("cov","").replace("read","").replace("seed:","")
            d = dict(zip(fields, name.split("|")))
            # want to shorten the read names
            d['id'] = "contig_%d" % i
            d['percent_of_total'] = 100 * (int(d['reads']) / total)
            meta.write("\t".join(map(str, [d[o] for o in out])) + "\n")
            write_fasta(fasta_out, d['id'], seq.upper())

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta_in", help="fasta with detailed read names")
    p.add_argument("fasta_out", help="fasta with unique, truncated read names")
    p.add_argument("meta", help="metadata for use in topiary explorer")
    main(p.parse_args())