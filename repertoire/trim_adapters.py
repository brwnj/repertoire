#!/usr/bin/env python
# encoding: utf-8
"""
Trim primer adapters from forward and reverse reads and update read name with
identifier.
"""
import sys
import gzip
import string
import editdist as ed
from toolshed import nopen
from itertools import izip, groupby, islice

COMPLEMENT = string.maketrans('ACGTNSRYMKWHBVD','TGCANSRYMKWHBVD')

def readfq(fq):
    class Fastq(object):
        def __init__(self, args):
            self.name = args[0][1:]
            self.seq = args[1]
            self.qual = args[3]
            assert len(self.seq) == len(self.qual)
    
        def __repr__(self):
            return "Fastq({name})".format(name=self.name)
    
        def __str__(self):
            return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                    seq=self.seq, qual=self.qual)
    
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)

def readfa(fa):
    class Fasta(object):
        def __init__(self, name, seq):
            self.name = name
            self.seq = seq
    
        def __repr__(self):
            return "Fasta({name})".format(name=self.name)
    
        def __str__(self):
            return ">{name}\n{seq}".format(
                        name=self.name,
                        seq="\n".join([self.seq[i:i + 70] for i in \
                                range(0, len(self.seq), 70)]))
    
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield Fasta(name, seq)

def rev_comp(seq):
    """return reverse complement of seq."""
    return seq.translate(COMPLEMENT)[::-1]

def fasta_to_dict(fasta):
    """intention is to save the primer sequences."""
    d = {}
    for rd in readfa(fasta):
        d[rd.name] = rd.seq
    return d

def get_primer(query, primers, n):
    """return the primer name and the length to trim."""
    primer = ""
    distance = n + 1
    for name, target in primers.iteritems():
        d = ed.distance(query[:len(target)], target)
        if d < distance:
            distance = d
            primer = name
    if distance < n:
        return primer, len(primers[primer])
    else:
        return False, False

def trim_loc(a, b):
    """find best edit distance and return its index else return length of a"""
    la, lb = len(a), len(b)
    dists = []
    for i in xrange(0, la-lb+1):
        dists.append(ed.distance(a[i:i+lb], b))
    best = min(dists)
    # 20% mismatch okay for now
    return dists.index(best) if best < .2*lb else la

def get_name(name, insert):
    if ".fastq" in name:
        sample = name.split(".fastq")[0]
    else:
        sample = name.split(".fq")[0]
    return "{sample}.{insert}.fastq.gz".format(**locals())

def main(args):
    mmatch = args.mismatches
    minleng = args.minlength
    r1_primers = fasta_to_dict(args.r1primer)
    r2_primers = fasta_to_dict(args.r2primer)
    r1_out = gzip.open(get_name(args.r1, "rmadptr"), 'wb')
    r2_out = gzip.open(get_name(args.r2, "rmadptr"), 'wb')
    for i, (r1, r2) in enumerate(izip(readfq(args.r1), readfq(args.r2)), start=1):
        if i % 100000 == 0: print >>sys.stderr, ">> processed %d reads" % i
        assert r1.name.split()[0] == r2.name.split()[0]

        # determine primer being used, trim location
        p1, r1_left_trim = get_primer(r1.seq, r1_primers, mmatch)
        p2, r2_left_trim = get_primer(r2.seq, r2_primers, mmatch)
        if not p1 or not p2: continue

        # find start of RC of primer in opposing sequence
        r1_right_trim = trim_loc(r1.seq[r1_left_trim:], rev_comp(r2_primers[p2]))
        r2_right_trim = trim_loc(r2.seq[r2_left_trim:], rev_comp(r1_primers[p1]))
        
        r1.name = "{id}:{cregion}:{fwork} 1".format(id=r1.name.split()[0], cregion=p1, fwork=p2)
        r2.name = "{id}:{cregion}:{fwork} 2".format(id=r2.name.split()[0], cregion=p1, fwork=p2)
        
        # do the trimming of seq and qual
        r1_full_trim = r1_right_trim + r1_left_trim
        r1.seq = r1.seq[r1_left_trim:r1_full_trim]
        r1.qual = r1.qual[r1_left_trim:r1_full_trim]
        r2_full_trim = r2_right_trim + r2_left_trim
        r2.seq = r2.seq[r2_left_trim:r2_full_trim]
        r2.qual = r2.qual[r2_left_trim:r2_full_trim]
        if len(r1.seq) < minleng or len(r2.seq) < minleng: continue

        # write the records
        r1_out.write(r1.__str__() + "\n")
        r2_out.write(r2.__str__() + "\n")

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("r1", help="read 1 fastq")
    p.add_argument("r2", help="read 2 fastq")
    p.add_argument("r1primer", help="expected R1 primer fasta")
    p.add_argument("r2primer", help="expected R2 primer fasta")
    p.add_argument("--mismatches", type=int, default=3, help="number of \
            mismatches to allow while matching primers")
    p.add_argument("--minlength", type=int, default=30, help="minimum \
            acceptable sequence length after trimming")
    args = p.parse_args()
    main(args)
