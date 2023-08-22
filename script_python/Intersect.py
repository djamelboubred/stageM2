#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Djamel Boubred, stage M2 IRD DIADE

"""
Intersect between 2 files by genes tags. In input A a tsv file with in firts columns the gene tags and in input B a output of SelectedKmerIntoGens.py

example commande lines :  python Intersect.py -wa ../../../stageM2/stageM2/input/diff_raws_counts_rename.tsv -wb test_DE_stats_GENE_INTERSECT.tsv
"""

import numpy as np
import csv
import pandas as pd
import argparse
from collections import Counter
import statistics as stat
from tqdm import tqdm
import random
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-wa", "--fileA", type=str, help="File with counts for each sample, with in the first colums the tags, then counts for each samples, in tsv format")
parser.add_argument("-wb", "--fileB", type=str, help="Minimum number of samples to support a k-mer [default: 20% of the size of the input condition].")
parser.add_argument("-o", "--output", type=str, help="output file name without extension, [default -o = counts_reccurence_threshold]")


args = parser.parse_args()
if args.fileA and args.fileB:
    if args.output:
        if "." in args.output:
            output=args.output
        else:
            output=args.output+".tsv"
    else:
        output="IntersectAnalysisA-B.tsv"

    fileA = pd.read_csv(args.fileA, sep='\t')
    fileB = pd.read_csv(args.fileB, sep='\t')

    # Effectuer la fusion sur les colonnes différentes
    merged_AB = pd.merge(fileA, fileB, left_on=fileA.columns[0], right_on=fileB.columns[0])

    print(merged_AB.head(5))
    print(f"The numbers of matching tags between {args.fileA} and {args.fileB} is : {len(merged_AB)}")
    merged_AB.to_csv(output, sep='\t', index=False)



    
    

