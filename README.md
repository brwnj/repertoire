#iSSAKE

Profiling model T-cell metagenomes with short reads

#Using iSSAKE

1. Quality trim your fastq

    seqtk trimfq in.fq > trimmed.fq

2. Convert fastq to fasta

    bioawk -c fastx '{print ">"$name"\n"$seq}' trimmed.fq > trimmed.fa

3. Download TCRB predictions from IMGT

4. Align reads using exonerate to find seeds

    exonerate --bestn 1 --score 1 --percent 0 -q trimmed.fa -t TCRB.fa 

#TODO
