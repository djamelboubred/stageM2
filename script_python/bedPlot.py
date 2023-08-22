#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Djamel Boubred, stage M2 IRD DIADE

"""
PLOT statistics as number of gene (or kmer) by regions from a output bed_analysis or tsv file,
Input tsv file and tags which is a basename of chromosomes in samples analyze. if regions is differents of basename regions are classifield in tags PAN
"""


import matplotlib.pyplot as plt
import csv
import pandas as pd
import argparse
import plotly.graph_objs as go
import plotly.express as px
import plotly.io as pio
from tqdm import tqdm
import pdfkit
#import plotly.offline as pyo


parser = argparse.ArgumentParser()

parser.add_argument("-tsv", "--tsv", type=str, help="tsv file from bedAnalysis.py, the format of file is with only two colums:\n in first colums region and second colums reccurence")
parser.add_argument("-o", "--output", type=str, help="output file name without extension, [default -o = GENE_INTERSECT]")
#parser.add_argument("-nb_filter", "--nb_filter", type=int, help="Given the numbers filter for concidere a regions in analyse (numbers of gene or kmer in regions concidere)(execpt Chromosomes with baseChr)")
parser.add_argument("-colnames", "--colnames", nargs='+', help="Given the colnames of file\n ex: -colnames REGIONS nb_gene")
parser.add_argument("-baseChr", "--baseChr", type=str, help="Given the prefix of Chromosomes names in your samples")
parser.add_argument("-baseUN", "--baseUN", type=str, help="Given the prefix of Unknown Chromosomes names in your samples")
args = parser.parse_args()

if args.tsv and args.colnames and args.baseChr:
    if args.output:
        if "." in args.output:
            output=args.output
        else:
            output=args.output+".pdf"
            output=args.output+".pdf"
    else:
        output=args.colnames[1]+"_histogramme.pdf"
    if len(args.colnames) == 2:
        tsv=pd.read_csv(args.tsv, delimiter='\t')
        regions=tsv[args.colnames[0]]
        nb_occurence=tsv[args.colnames[1]]
        nb_pan=tsv[~tsv[args.colnames[0]].str.contains(args.baseChr)][args.colnames[1]].sum()
        chr="^"+args.baseChr+"[0-9]+"
        df_CHR=tsv[tsv[args.colnames[0]].str.contains(chr)]
        df_PAN=pd.DataFrame({args.colnames[0]:["PAN"],args.colnames[1]:[nb_pan]})
        if args.baseUN:
            nb_UN=tsv[tsv[args.colnames[0]].str.contains(args.baseUN)][args.colnames[1]].sum()
            df_UN=pd.DataFrame({args.colnames[0]:["UNKNOWN"],args.colnames[1]:[nb_UN]})
            fig = px.histogram(pd.concat([df_CHR,df_UN,df_PAN], ignore_index=True), x=args.colnames[0], y=args.colnames[1])
            fig.show()
        else:
            fig = px.histogram(pd.concat([df_CHR,df_PAN], ignore_index=True), x=args.colnames[0], y=args.colnames[1])
            
    else:
        print("****\n\n-colnames argument accept only 2 colums names (because the input file has constituate only 2 colums). Please Check your file or your argument.\n\n****")
        print(" (°~/°) ")


#plt.figure()
## Définition de la largeur des barres
#width = 0.4
## Extraire les clés et les valeurs du dictionnaire
#cles = list(dico_gene_in_region_count.keys())
#valeurs = list(dico_gene_in_region_count.values())
#    
#
## Créer le graphique à barres
#plt.bar(cles, valeurs, width=width)
#plt.xlabel("Localisation")
#plt.ylabel("Nb gene")
#title='Numbers of genes for each regions'
#plt.title(title)
#
#for i in range(len(cles)):
#    plt.text(cles[i], valeurs[i], str(valeurs[i]), ha='center', va='bottom')
## Sauvegader le graphique
#plt.savefig("NumbersOfGeneByRegions.pdf") 

