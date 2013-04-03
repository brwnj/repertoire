#!/usr/bin/env python
# encoding: utf-8
"""
Trim primer adapters from forward and reverse reads and update read name with
identifier.
"""
import sys
import numpy as np
import editdist as ed
from toolshed import nopen
from itertools import izip
from parsers import read_fastx
from Bio import pairwise2

BASE_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'S':'S',
                    'R':'R', 'Y':'Y', 'M':'M', 'K':'K', 'W':'W', 'N':'N',
                    'H':'H', 'B':'B', 'V':'V', 'D':'D'}

def rev_comp(seq):
    """return reverse complement of seq."""
    return ''.join([BASE_COMPLEMENT[b] for b in seq[::-1]])

def fasta_to_dict(fasta):
    """intention is to save the primer sequences."""
    d = {}
    with nopen(fasta) as fa:
        for name, seq, qual in read_fastx(fa):
            d[name] = seq
    return d

def get_primer(query, primers, mismatches):
    """return the primer name and the length to trim."""
    match_name = ""
    distance = 99
    for name, target in primers.iteritems():
        d = ed.distance(query[:len(target)], target)
        if d < distance:
            distance = d
            match_name = name
    # return False if distance to length ratio is terrible
    if distance < mismatches:
        return match_name, len(primers[match_name])
    else:
        return False, False

def get_trim_loc(target, query):
    """docstring for get_trim
    return the trim location if possible or the sequence length
    match=1, mismatch=-1, gapopen=-5, gapextend=-3
    """
    matches = pairwise2.align.localms(target, query, 1, -1, -3, -2)
    try:
        # highest scoring match first
        return int(matches[0][3])
    except IndexError:
        # no match
        return len(target)

def main(args):
    # read 1 and 2 specific primers
    r1_primers = fasta_to_dict(args.r1primer)
    r2_primers = fasta_to_dict(args.r2primer)

    with nopen(args.r1) as r1, nopen(args.r2) as r2,\
            open(args.r1out, 'wb') as r1out, open(args.r2out, 'wb') as r2out:

        for i, ((r1name, r1seq, r1qual), (r2name, r2seq, r2qual)) in \
                enumerate(izip(read_fastx(r1), read_fastx(r2)), start=1):
            if i % 100000 == 0:
                print >> sys.stderr, ">> processed %d reads" % i

            # this should never happen after `filter_pairs.py`
            if r1name.split()[0] != r2name.split()[0]:
                print >> sys.stderr, ">> r1 and r2 are out of sync."
                sys.exit(1)

            # determine primer being used
            r1_pname, r1_left_trim = get_primer(r1seq, r1_primers, args.mismatches)
            r2_pname, r2_left_trim = get_primer(r2seq, r2_primers, args.mismatches)

            # unable to determine primer; skip this read
            if not r1_pname or not r2_pname:
                # maybe dump into badprimer.fastq
                # print >> sys.stderr, ">> undetermined primer for %s" % r1name
                continue

            # find start of RC of primer in opposing sequence
            r1_right_trim = get_trim_loc(r1seq[r1_left_trim:], rev_comp(r2_primers[r2_pname]))
            r2_right_trim = get_trim_loc(r2seq[r2_left_trim:], rev_comp(r1_primers[r1_pname]))
            
            # read count:c-region:vh-framework
            r1name = "read_%d:%s:%s 1" % (i, r1_pname, r2_pname)
            r2name = "read_%d:%s:%s 2" % (i, r1_pname, r2_pname)
            r1seq = r1seq[r1_left_trim:r1_right_trim + r1_left_trim]
            r1qual = r1qual[r1_left_trim:r1_right_trim + r1_left_trim]
            r2seq = r2seq[r2_left_trim:r2_right_trim + r2_left_trim]
            r2qual = r2qual[r2_left_trim:r2_right_trim + r2_left_trim]
            if len(r1seq) < args.minlength or len(r2seq) < args.minlength:
                # count number discarded due to length
                continue
            # write the records
            r1out.write("@%s\n%s\n+\n%s\n" % (r1name, r1seq, r1qual))
            r2out.write("@%s\n%s\n+\n%s\n" % (r2name, r2seq, r2qual))
            
if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("r1", help="read 1 fastq")
    p.add_argument("r2", help="read 2 fastq")
    p.add_argument("r1primer", help="expected R1 primer fasta")
    p.add_argument("r2primer", help="expected R2 primer fasta")
    p.add_argument("r1out", help="read 1 output fastq")
    p.add_argument("r2out", help="read 2 output fastq")
    p.add_argument("--mismatches", type=int, default=6, help="number of \
            mismatches to allow while matching primers of the 5' end \
            [ %(default)s ]")
    p.add_argument("--minlength", type=int, default=10, help="minimum \
            acceptable sequence length after trimming [ %(default)s ]")
    # lao = p.add_argument_group("local alignment options")
    # lao.add_argument("--match", type=int, default=1, help="match score [ %(default)s ]")
    # lao.add_argument("--mismatch", type=int, default=-1, help="mismatch penalty [ %(default)s ]")
    # lao.add_argument("--gapopen", type=int, default=-5, help="gap open penalty [ %(default)s ]")
    # lao.add_argument("--gapextent", type=int, default=-3, help="gap extend penalty [ %(default)s ]")
    main(p.parse_args())