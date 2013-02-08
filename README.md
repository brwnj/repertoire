#iSSAKE

Profiling model T-cell metagenomes with short reads

#Using iSSAKE

.gz is supported throughout the pipeline

1. Quality trim your fastq

    `seqtk trimfq in.fq > trimmed.fq`

1. Convert fastq to fasta

    `bioawk -c fastx '{print ">"$name"\n"$seq}' trimmed.fq > trimmed.fa`

1. Download TCRB predictions from IMGT

1. Create tags from IMGT regions

    `python create_tags.py -v -l 35 trav.fa > trav.tags.fa`
    
    * result must be a substring within the reads. too long and you won't find
    any reads. too short and you'll lose regions as a unique substring among them
    will not found.

1. Find seeds among your reads

    `python find_seeds.py -v trav.tags.fa in.fq > seeds.fa`

1. Run iSSAKE

    `iSSAKE -f trimmed.fa -s seeds.fa`

#Links

Bioawk: https://github.com/lh3/bioawk

Python dependency: ``pip install toolshed``

TRBV group: http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+TRBV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON

TRAV group: http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+TRAV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON