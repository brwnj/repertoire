#!/usr/bin/env python
# encoding: utf-8
"""
Parse the output of SeqPrep joined reads when using trim_adapter.py's output.
"""
from toolshed import nopen
from parsers import read_fastx
from collections import Counter

def main(args):
    for fastq in args.fastq:
        print fastq
        meta = Counter()
        with nopen(fastq) as fh:
            for name, seq, qual in read_fastx(fh):
                name, cregion, fwork = name.split()[0].split(":")
                meta.update(["%s:%s" % (cregion, fwork)])
        for combination, count in meta.iteritems():
            print "%s\t%d" % (combination, count)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fastq', nargs="+", help="SeqPrep joined reads.")
    args = p.parse_args()
    main(args)
