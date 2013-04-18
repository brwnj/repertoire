#!/usr/bin/env python
# encoding: utf-8
"""
Find unique reads by primer among all joined reads.
"""
import sys
import editdist as ed
from toolshed import nopen
from parsers import read_fastx
from collections import OrderedDict

def fastq_to_dict(fastq):
    """docstring for fastq_to_dict"""
    d = {}
    with nopen(fastq) as fh:
        for name, seq, qual in read_fastx(fh):
            d[name] = {'seq':seq,'qual':qual}
    return d

def distance(a, b):
    """Calculate Levenshtein distance accounting for length differences between
    the two strings. Returns int.
    >>> distance('CATGGGTGGTTCAGTGGTAGAATTCTCGCCTGCC', 'GTGCTGTAGGCATT')
    2
    """
    return ed.distance(a, b) - abs(len(a) - len(b))

def main(args):
    fd = fastq_to_dict(args.fastq)
    # longest sequence first
    fd = OrderedDict(sorted(fd.items(), key=lambda (k, v): len(v['seq']), reverse=True))
    seen = set()
    ignore = set()
    for i, (t_name, target) in enumerate(fd.iteritems(), start=1):
        # skip reads that have already been determined as subsequences
        if t_name in ignore: continue
        seen.add(t_name)
        if i % 10 == 0:
            print >> sys.stderr, ">> processed %d reads..." % i
        t_id, t_cregion, t_fwork = t_name.split(":")
        for q_name, query in fd.iteritems():
            # seen - tested in first loop
            # ignore - already found to be subsequence of a target
            if q_name in seen or q_name in ignore: continue
            q_id, q_cregion, q_fwork = q_name.split(":")
            # only attempt to collapse things of the same c-region and framework
            if t_cregion != q_cregion and t_fwork != q_fwork: continue
            # compare the sequences
            if distance(target['seq'], query['seq']) < args.mismatches:
                # mark subsequences to be ignored
                ignore.add(q_name)
        print "@%s\n%s\n+\n%s" % (t_name, target['seq'], target['qual'])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fastq', help="reads to collapse to unique")
    p.add_argument('-m', '--mismatches', type=int, default=0,
            help="mismatches to allow during mapping [%(default)s]")
    main(p.parse_args())