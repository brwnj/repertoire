#!/usr/bin/env python
# encoding: utf-8
"""
matches = pairwise2.align.localms(target, query, 1, -1, -3, -2)
try:
    # highest scoring match first
    return int(matches[0][3])
except IndexError:
"""
import sys
from toolshed import nopen
from parsers import read_fastx
from Bio import pairwise2
from collections import OrderedDict

def fastq_to_dict(fastq):
    """docstring for fastq_to_dict"""
    d = {}
    with nopen(fastq) as fh:
        for name, seq, qual in read_fastx(fh):
            d[name] = {'seq':seq,'qual':qual}
    return d

def main(args):
    fd = fastq_to_dict(args.fastq)
    # convert to ordered dictionary
    fd = OrderedDict(sorted(fd.items(), key=lambda (k, v): len(v['seq'])))
    seen = {}
    for i, (name, query) in enumerate(fd.iteritems(), start=1):
        if i % 1000 == 0:
            print >> sys.stderr, ">> processed %d reads..." % i
        subseq = False
        q_id, q_cregion, q_fwork = name.split(":")
        expected_score = len(query['seq']) - args.mismatches
        # maps onto same length or longer seqs
        for t_name, target in fd.iteritems():
            if t_name == name: continue
            # skipping reads we've already mapped
            if seen.has_key(t_name): continue
            t_id, t_cregion, t_fwork = t_name.split(":")
            # only attempt to collapse things of the same c-region and framework
            if q_cregion != t_cregion and q_fwork != t_fwork: continue
            # locally align using smith-waterman
            matches = pairwise2.align.localms(target['seq'], query['seq'], 1, -1, -1, -1)
            high_score = matches[0][2]
            if high_score == expected_score:
                subseq = True
                break
        if not subseq:
            # print fastq record
            print "@%s\n%s\n+\n%s" % (name, query['seq'], query['qual'])
        seen[name] = ""

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fastq', help="reads to collapse to unique")
    p.add_argument('-m', '--mismatches', type=int, default=0,
            help="mismatches to allow during mapping [ %(default)s ]")
    main(p.parse_args())