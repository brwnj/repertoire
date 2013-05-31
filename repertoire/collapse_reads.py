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
    """find best edit distance between two strings of potentially uneven length.
    """
    la, lb = len(a), len(b)
    if la < lb:
        return distance(b, a)
    if la == lb:
        return editdist.distance(a, b)
    else:
        dists = []
        for i in xrange(0, la-lb+1):
            dists.append(editdist.distance(a[i:i+lb], b))
        return min(dists)

def main(args):
    fd = fastq_to_dict(args.fastq)
    # longest sequence first
    fd = OrderedDict(sorted(fd.items(), key=lambda (k, v): len(v['seq']), reverse=True))
    seen = set()
    ignore = set()
    i = 1
    for t_name, target in fd.iteritems():
        if i % 10 == 0:
            print >> sys.stderr, ">> processed %d reads..." % i
        seen.add(t_name)
        i += 1
        # skip reads that have already been determined as subsequences
        if t_name in ignore: continue
        t_id, t_cregion, t_fwork = t_name.split(":")
        for q_name, query in fd.iteritems():
            if q_name in seen or q_name in ignore: continue
            q_id, q_cregion, q_fwork = q_name.split(":")
            # only attempt to collapse things of the same c-region and framework
            if t_cregion != q_cregion and t_fwork != q_fwork: continue
            if distance(target['seq'], query['seq']) < args.mismatches:
                ignore.add(q_name)
                i += 1
        print "@%s\n%s\n+\n%s" % (t_name, target['seq'], target['qual'])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fastq', help="reads to collapse to unique")
    p.add_argument('-m', '--mismatches', type=int, default=1,
            help="mismatches to allow during mapping [%(default)s]")
    args = p.parse_args()
    main(args)
