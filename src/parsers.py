#! /usr/bin/env python
import itertools

def read_fastq(fq):
    """parses fastq filehandle and returns name, sequence, and qualities."""
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
    """parses fasta filehandle and returns name and sequence."""
    for header, group in itertools.groupby(fa, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq

def write_fasta(fh, name, seq, wrap=70):
    """given filehandle, name, and sequence, writes fasta lines and wraps
    sequence lines at wrap setting.
    """
    fh.write(">%s\n" % name)
    # wrap the sequence lines at length `wrap`
    fh.write("\n".join([seq[i:i + wrap] for i in range(0, len(seq), wrap)]) + "\n")