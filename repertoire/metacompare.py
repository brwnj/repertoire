#!/usr/bin/env python
# encoding: utf-8
"""
Compiles per V:J stats and fold changes between samples.
"""

import sys
import itertools
import numpy as np
import pandas as pd

def scale_factors(df):
    median = np.median(df.sum(0))
    factors = {}
    for name, total in df.sum(0).iteritems():
        if total == median:
            factors[name] = 1.
        else:
            factors[name] = median / total
    return factors

def fold_change(a, b, counts):
    if counts[a] > counts[b]:
        return np.divide(counts[a], counts[b])
    else:
        return np.divide(counts[b], counts[a])

def main(args):
    # build the dataframe
    meta_df = pd.DataFrame()
    usecols = ["v_region", "j_region", "reads"]
    index = ["v_region", "j_region"]
    names = args.names.split(",")
    for name, file in itertools.izip(names, args.meta_files):
        print >> sys.stderr, ">> reading sample %s: %s" % (name, file)
        if meta_df.empty:
            # read into dataframe
            meta_df = pd.io.parsers.read_table(file, usecols=usecols)
            # sum v and j regions in pivot table because v_region is not unique
            meta_df = pd.pivot_table(meta_df, values="reads", rows=index, aggfunc=np.sum)
            # convert summed data back into dataframe indexed on v and j
            meta_df = meta_df.reset_index().set_index(index)
            # rename "reads" to sample name
            meta_df = meta_df.rename(columns={"reads":name})
        else:
            # same as above but we have to join to existing
            temp = pd.io.parsers.read_table(file, usecols=usecols)
            temp = pd.pivot_table(temp, values="reads", rows=index, aggfunc=np.sum)
            temp = temp.reset_index().set_index(index)
            temp = temp.rename(columns={"reads":name})
            # append the column onto table
            meta_df = meta_df.join(temp)
    # scale the data to total read counts per sample
    factors = scale_factors(meta_df)
    for name, factor in factors.iteritems():
        print >> sys.stderr, ">> sample %s scaling factor: %f" % (name, factor)
        meta_df[name] = meta_df[name] * factor
    # reset the index to simplify adding fold change to dataframe
    meta_df = meta_df.reset_index()
    # adding fold change columns
    for k, counts in meta_df.iterrows():
        for (a, b) in itertools.combinations(names, 2):
            pair = "fold_change(%s:%s)" % (a, b)
            try:
                # add fold change into dataframe as larger / smaller
                meta_df[pair][int(k)] = fold_change(a, b, counts)
            except KeyError:
                # create the column
                meta_df[pair] = np.nan
                meta_df[pair][int(k)] = fold_change(a, b, counts)
    meta_df.to_csv(sys.stdout, sep="\t", na_rep="NaN", index=False)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("names",
            help="comma delimited list of sample names defining the metadata files")
    p.add_argument("meta_files", nargs="+",
            help="metadata of the samples to compare")
    args = p.parse_args()
    main(args)
