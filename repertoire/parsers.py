#! /usr/bin/env python
import itertools

def read_fastx(fp):    
    """author: @lh3. parses fasta and fastq."""
    last = None
    while True:
        if not last:
            for l in fp:
                if l[0] in '>@':
                    last = l[:-1]
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp:
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':
            yield name, ''.join(seqs), None
            if not last: break
        else:
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):
                    last = None
                    yield name, seq, ''.join(seqs);
                    break
            if last:
                yield name, seq, None
                break

def write_fasta(fh, name, seq, wrap=70):
    """given filehandle, name, and sequence, writes fasta lines and wraps
    sequence lines at wrap setting.
    """
    fh.write(">%s\n" % name)
    # wrap the sequence lines at length `wrap`
    fh.write("\n".join([seq[i:i+wrap] for i in range(0,len(seq),wrap)])+"\n")
