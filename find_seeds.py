#! /usr/bin/env python
# encoding: utf-8
"""
given sequence tags, finds reads containing that sequence and creates a fasta.
"""
import sys
import os.path as op
from toolshed import nopen

class read_fastq(object):
    """Yields name, seq, qual."""
    def __init__(self, fastq):
        self.fastq = nopen(fastq)

    def __iter__(self):
        fq = self.fastq
        while True:
            id1 = fq.next().strip()
            seq = fq.next().strip()
            id2 = fq.next().strip()
            qual = fq.next().strip()
            if qual == "":
                if id1 != "":
                    sys.stderr.write(">> Incomplete fastq record... skipping.\n")
                break
            yield id1[1:], seq, qual

def read_fasta(fa):
    """yields name and seq from fasta."""
    name, seq = None, []
    for record in nopen(fa):
        record = record.rstrip()
        if record.startswith('>'):
            if name: yield (name, ''.join(seq))
            name, seq = record.lstrip('>'), []
        else:
            seq.append(record.upper())
    if name: yield (name, ''.join(seq))

def main(args):
    tags = {}
    if args.verbose:
        sys.stderr.write(">> reading in tag sequences...\n")
    for name, seq in read_fasta(args.tags):
        tags[name] = seq
    i = 0
    for fq in args.fastq:
        if args.verbose:
            sys.stderr.write(">> processing %s...\n" % op.basename(fq))
        for fq_id, fq_seq, fq_qual in read_fastq(fq):
            i += 1
            if i%1000000 == 0 and args.verbose:
                sys.stderr.write(">> processed %d reads...\n" % i)
            for tcr_name, tcr_seq in tags.iteritems():
                if fq_seq.find(tcr_seq) != -1:
                    print ">%s|%s\n%s" % (fq_id, tcr_name, fq_seq)
                    break

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            usage="%(prog)s [options] tags fastq",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tags", help="unique tag sequences as fasta")
    p.add_argument("fastq", nargs="+", help="read pool to create unique seeds")
    p.add_argument("-v", "--verbose", action="store_true",
            help="maximum verbosity")
    main(p.parse_args())