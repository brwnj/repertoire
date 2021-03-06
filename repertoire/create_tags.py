#! /usr/bin/env python
# encoding: utf-8
"""
Find unique sequence tags among regions (see IMGT for downloads). Discards
any region with non-unique tag.
"""
import sys
import itertools
from toolshed import nopen
from parsers import read_fasta

def main(args):
    full_seqs = {}
    tags = {}
    with nopen(args.fasta) as fasta:
        for name, seq in read_fasta(fasta):
            full_seqs[name] = seq.upper()
        
    # each tcr in original fasta
    for tcr, seq in full_seqs.iteritems():
        unique_expected = len(full_seqs) - 1
        # all possible tags
        for i in range(len(seq) - args.length, 0, -1):
            # tag matches favor 3' end
            tag = seq[i:args.length + i]
            # reached the end of the sequence
            if len(tag) < args.length: break
            unique_found = 0

            # tag not present in any other tcr
            for ss_tcr, ss_seq in full_seqs.iteritems():
                # the current tcr
                if ss_tcr == tcr: continue
                # finding unique tags
                if ss_seq.find(tag) == -1:
                    unique_found += 1

            if unique_found == unique_expected:
                tags[tcr] = tag
                # exit loop on first unique tag
                break
                
    # ensure this actually worked
    taglist = [tag for name, tag in tags.iteritems()]
    tagset = set(taglist)
    assert(len(taglist) == len(tagset))
    taglist = []
    
    # print results
    for tcr, tag in tags.iteritems():
        taglist.append(tcr)
        print ">%s\n%s" % (tcr, tag)

    # tags found stats
    if args.verbose:
        sys.stderr.write("Of %d regions, %d tags were found.\n" \
                            % (len(full_seqs), len(tagset)))
        alltcrs = []
        for tcr, seq in full_seqs.iteritems():
            alltcrs.append(tcr)
        alltcrs = set(alltcrs)
        taglist = set(taglist)
        diff = alltcrs - taglist
        if diff:
            sys.stderr.write("Unable to find a unique tag for:\n")
            sys.stderr.write("\n".join(diff) + "\n")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            usage="%(prog)s [options] fasta",
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", metavar="FASTA", help="regions fasta")
    p.add_argument("-l", "--length", default=35, type=int,
            help="desired tag length [ %(default)s ]")
    p.add_argument("-v", "--verbose", action="store_true",
            help="maximum verbosity")
    main(p.parse_args())