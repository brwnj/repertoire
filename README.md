#iSSAKE

Profiling model T-cell metagenomes with short reads

#Using iSSAKE

1. Quality trim your fastq

    ``seqtk trimfq in.fq > trimmed.fq``

2. Convert fastq to fasta

    ``bioawk -c fastx '{print ">"$name"\n"$seq}' trimmed.fq > trimmed.fa``

3. Download TCRB predictions from IMGT

4. Create tags from IMGT

    ``python create_tags.py -v -l 40 trav.fa > trav.tags.fa``

5. Run iSSAKE