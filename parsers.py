#! /usr/bin/env python
import itertools

def read_fastq(fq):
    while True:
        values = list(itertools.islice(fq, 4))
        if len(values) == 4:
            id1, seq, id2, qual = values
        elif len(values) == 0:
            raise StopIteration
        else:
            raise EOFError("unexpected end of file")
        assert id1.startswith('@'),\
                ">> Fastq out of sync at read:\n%s\n" % id1
        assert id2.startswith('+'),\
                ">> Fastq out of sync at read:\n%s\n" % id1
        assert len(seq) == len(qual),\
                ">> Sequence and Quality are not the same length \
                for read:\n%s\n" % id1
        yield id1[1:-1], seq[:-1], qual[:-1]

def read_fasta(fa):
    for header, group in itertools.groupby(fa, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq
