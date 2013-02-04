#iSSAKE

Profiling model T-cell metagenomes with short reads

#Using iSSAKE

1. Quality trim your fastq

    ``seqtk trimfq in.fq > trimmed.fq``

1. Convert fastq to fasta

    ``bioawk -c fastx '{print ">"$name"\n"$seq}' trimmed.fq > trimmed.fa``

1. Download TCRB predictions from IMGT

1. Create tags from IMGT

    ``python create_tags.py -v -l 35 trav.fa > trav.tags.fa``

1. Find seeds among your reads

    ``python find_seeds.py -v trav.tags.fa in.fq > seeds.fa``

1. Run iSSAKE