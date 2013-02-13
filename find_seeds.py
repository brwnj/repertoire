#! /usr/bin/env python
# encoding: utf-8
"""
Given sequence tags, finds reads containing that sequence and creates a fasta.
"""
import sys
import itertools
import os.path as op
from toolshed import nopen
from parsers import read_fastq, read_fasta

def print_record(tags, f_id, f_seq):
    for tcr_name, tcr_seq in tags.iteritems():
        if f_seq.find(tcr_seq) != -1:
            print ">%s|%s\n%s" % (f_id, tcr_name, f_seq)

def main(args):
    tags = {}
    if args.verbose:
        sys.stderr.write(">> reading in tag sequences...\n")
    for name, seq in read_fasta(args.tags):
        tags[name] = seq
    i = 0
    for fx in args.reads:
        if args.verbose:
            sys.stderr.write(">> processing %s...\n" % op.basename(fx))
        # process either fasta or fastq.
        if ".fasta" in fx or ".fa" in fx:
            for f_id, f_seq in read_fasta(fx):
                i += 1
                if i%1000000 == 0 and args.verbose:
                    sys.stderr.write(">> processed %d reads...\n" % i)
                print_record(tags, f_id, f_seq)
        else:
            for f_id, f_seq, f_qual in read_fastq(fx):
                i += 1
                if i%1000000 == 0 and args.verbose:
                    sys.stderr.write(">> processed %d reads...\n" % i)
                print_record(tags, f_id, f_seq)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            usage="%(prog)s [options] tags fastq",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tags", help="unique tag sequences as fasta")
    p.add_argument("reads", nargs="+", help="read pool to create unique seeds")
    p.add_argument("-v", "--verbose", action="store_true",
            help="maximum verbosity")
    main(p.parse_args())