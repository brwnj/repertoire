#!/usr/bin/env python
# encoding: utf-8
"""
Parses the read name down to only include the necessary gene label.
"""
import re
import sys
import itertools
from toolshed import nopen
from parsers import read_fasta

def main(args):
    for name, seq in read_fasta(nopen(args.fasta)):
        try:
            # rename from imgt
            name = re.findall(r'(%s[^\|]+)' % args.gene.upper(), name)[0]
            print ">%s\n%s" % (name, seq.upper())
        except IndexError:
            sys.stderr.write(">> unable to parse: %s\n>> for gene: %s\n" \
                                % (name, args.gene))
            pass

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fasta')
    req = p.add_argument_group('required arguments')
    req.add_argument('-g', '--gene', required=True,
            help="gene name, eg. TRAJ or TRBV")
    main(p.parse_args())