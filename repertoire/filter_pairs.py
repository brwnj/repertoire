#!/usr/bin/env python
# encoding: utf-8
"""
Remove reads from R1 and R2 that are missing a mate. Writes output in same order
for each file.
"""
import sys
import gzip
import os.path as op
from parsers import read_fastx
from toolshed import nopen

def fq_to_set(fq):
    fq_l = []
    fq_d = {}
    with nopen(fq) as fh:
        for name, seq, qual in read_fastx(fh):
            # split name on whitespace, leaving off "1" or "2"
            read_id = name.split()[0]
            fq_l.append(read_id)
            # a million reads ~ 1 GB
            fq_d[read_id] = "@%s\n%s\n+\n%s\n" % (name, seq, qual)
    return set(fq_l), fq_d

def get_file_name(file_name):
    path, ext = op.splitext(file_name)
    if ext == ".gz":
        path, ext = op.splitext(path)    
    return "%s.filtered%s" % (path, ext)

def main(args):
    # finding intersection of read headers
    print >> sys.stderr, ">> building sequence dictionaries"
    r1_set, r1_d = fq_to_set(args.r1)
    r2_set, r2_d = fq_to_set(args.r2)
    intersect = (r1_set & r2_set)
    # new file names
    filtered1 = get_file_name(args.r1)
    filtered2 = get_file_name(args.r2)
    # write new fastqs
    print >> sys.stderr, ">> writing output to %s and %s" % (filtered1, filtered2)
    with open(filtered1, 'wb') as f1, open(filtered2, 'wb') as f2:
        for i, name in enumerate(intersect):
            if i % 100000 == 0:
                print >> sys.stderr, ">> processed %d reads" % i
            f1.write(r1_d[name])
            f2.write(r2_d[name])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("r1", help="R1 fastq")
    p.add_argument("r2", help="R2 fastq")
    main(p.parse_args())