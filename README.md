#iSSAKE

Profiling model T-cell metagenomes with short reads

##Using iSSAKE

.gz is supported throughout the pipeline

Quality trim your fastq

```
seqtk trimfq in.fq > trimmed.fq
```

Convert fastq to fasta

```
bioawk -c fastx '{print ">"$name"\n"$seq}' trimmed.fq > trimmed.fa
```

Download TCRB predictions from IMGT ([TRAV][1] or [TRBV][2])

Create tags from IMGT regions

```
python create_tags.py -v -l 35 trav.fa > trav.tags.fa
```

Find seeds among your reads

```
python find_seeds.py -v trav.tags.fa in.fq > seeds.fa
```

Run iSSAKE

```
iSSAKE -f trimmed.fa -s seeds.fa -b sampleid
```

##Further analysis

Download J regions based on strand ([TRAJ][3] or [TRBJ][4]).

Rename fasta names

```
python renameIMGT.py --gene TRAJ imgt_traj.fa > traj.fa
```

Locally align J regions to assembled contigs

```
exonerate -q sampleid.contigs \
    -t traj.fa \
    --bestn 1 \
    --ryo ">%qi|%ti\n%qs" \
    --showalignment FALSE \
    --showvulgar FALSE \
    > sampleid.jregion.fa
```

That will add the J region name onto the read the name

Parse read names into data table

Align sequences into tree using [Muscle][5] 

Visualize data in [Topiary Explorer][6]

##Links

Bioawk: https://github.com/lh3/bioawk

Python dependency: ``pip install toolshed``

[1]: http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+TRAV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON
[2]: http://www.imgt.org/IMGT_GENE-DB/GENElect?query=8.1+TRBV&species=Homo+sapiens&IMGTlabel=L-PART1+V-EXON
[3]: http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.2+TRAJ&species=Homo+sapiens
[4]: http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.2+TRBJ&species=Homo+sapiens
[5]: http://www.ebi.ac.uk/Tools/msa/muscle/
[6]: https://github.com/qiime/Topiary-Explorer
