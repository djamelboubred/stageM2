#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Djamel Boubred, stage M2 IRD DIADE

"""
Calculate statistics from a bedfile after launch this samtools command line
bedtools intersect -abam aln_ref_vs_kmers.bam -b genes.gff -bed -wb -wa > intersection_btw_aligned_kmers_and_genes.bed
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import argparse
from collections import Counter
import statistics as stat
from tqdm import tqdm

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bed", type=str, help="Bed file from bedtools intersect, the format of tags should be 12354_tags_AATTATTT or AATTAAT\n If your tags is a only numerical tags this programs run")
parser.add_argument("-s", "--similar", action='store_true',help="Create a column called Similar if this word is found in bed file")
parser.add_argument("-o", "--output", type=str, help="output file name without extension, [default -o = GENE_INTERSECT]")
parser.add_argument("-mapq", "--mapq", type=int, help="mapq filter values, [default mapq = 15]")
parser.add_argument("-k", "--nbkmer", type=int, help=" numbers of kmer reccurence for concider a gene, [default k = 10]")
args = parser.parse_args()

if args.bed:
    if args.nbkmer:
        nb_kmers_filter=args.nbkmer
    else:
        nb_kmers_filter = 10
    print(f"Numbers ok kmer reccurence : {nb_kmers_filter}")
    if args.mapq:
        mapq=args.mapq
    else:
        mapq=15
    dico_gene_in_contigs={}
    liste_kmer_in_contigs=[]
    print(f"Mapping quality : {mapq}")
    bed_file=pd.read_csv(args.bed, sep='\t')
    #pandas storage of bed
    bed=pd.DataFrame(bed_file)
    print(f"\n\n*******************\n\nSTEP: \t\t NB KMERS TO GENE\n\n{args.bed} file is load\n\n*******************")

    dico_count_kmer_in_gene={}
    dico_gene_kmer= {}
    gene_identique=False # booléen dans le cas ou le gene est le même gène que la ligne d'avant on a la valeur True
    nb_kmer=0 #itération qui compte le nombre de kmer total qui sont pris en compte en fonction de nos kmer
    
    # iteration in bed
    for i in tqdm(range(len(bed))):
        # counting kmers by gene filtering by mapq
        if int(bed.iloc[i,4]) > mapq:
            nb_kmer+=1
            contigs=bed.iloc[i,0]
            gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]

            if "_tags_" in str(bed.iloc[i,3]):
                kmer=bed.iloc[i, 3].strip().split('_')[2]
            else:
                #kmer name recovery 
                kmer=bed.iloc[i,3]
            if gene_name in dico_gene_kmer:
                # La clé existe, ajouter le k-mer à la liste existante
                dico_gene_kmer[gene_name]['kmer'].append(kmer)
                
            else:
                # La clé n'existe pas, créer une nouvelle entrée avec la clé et le k-mer
                gene_start=int(bed.iloc[i, 15])
                gene_stop=int(bed.iloc[i, 16])
                length=gene_stop-gene_start

                # adding informations about annotation if available
                if args.similar:
                    if bed.iloc[i,20].strip().split(';')[2].split('=')[0]=="Note":
                        gene_similar = bed.iloc[i,20].strip().split(';')[2].split('=')[1].split(':')[0].split(' ')[2]            
                        note_function = bed.iloc[i,20].strip().split(';')[2].split('=')[1]
                    else:
                        gene_similar = bed.iloc[i,20].strip().split(';')[3].split('=')[1].split(':')[0].split(' ')[2]
                        note_function = bed.iloc[i,20].strip().split(';')[3].split('=')[1]
                    dico_gene_kmer[gene_name] = {'kmer': [kmer],'length':length, 'contigs':contigs,'START': gene_start, 'STOP': gene_stop, 'Similar':gene_similar,'Note': note_function}
                else:
                    dico_gene_kmer[gene_name] = {'kmer': [kmer],'length':length, 'contigs':contigs, 'START': gene_start, 'STOP': gene_stop}

                #on l'ajoute à une liste qui va contenir toute les contigs de chaque gène
            if gene_name in dico_gene_in_contigs:
                dico_gene_in_contigs[gene_name]['contigs'].append(contigs)
            else:
                dico_gene_in_contigs[gene_name]={'contigs': [contigs]}
            #on l'ajoute à une liste qui va contenir toute les contigs de chaque kmer
            liste_kmer_in_contigs.append(contigs)

        #CALCULATE THE NUMBERS OF KMER FOR EACH GENES
        dico_count_kmer_in_gene[gene_name] = len(dico_gene_kmer[gene_name]['kmer'])
    
    # Créer un dico qui va contenir en clé la région et en valeurs le nombre de gène
    dico_gene_in_contigs_count={}
    for gene_name in dico_gene_in_contigs:
        if len(dico_gene_in_contigs[gene_name]['contigs']) >=nb_kmers_filter:
            element = dico_gene_in_contigs[gene_name]['contigs']
            if element[0] in dico_gene_in_contigs_count:
                # La clé existe, ajouter le gene à la liste existante
                dico_gene_in_contigs_count[element[0]]['nb_gene']+=1        
            else:
                # La clé n'existe pas, créer une nouvelle entrée avec la clé et le gene
                dico_gene_in_contigs_count[element[0]]={'nb_gene':1}

        # Créer un dico qui va contenir en clé la région et en valeurs le nombre de kmer
    dico_kmer_in_contigs_count={}
    for element in liste_kmer_in_contigs:
        if element in dico_kmer_in_contigs_count:
            # La clé existe, ajouter le kmer à la liste existante
            dico_kmer_in_contigs_count[element]['nb_kmer']+=1        
        else:
            # La clé n'existe pas, créer une nouvelle entrée avec la clé et le kmer
            dico_kmer_in_contigs_count[element]={'nb_kmer':1}
    
    nb_kmer_tot=len(bed)

    ####
    #SAVE the number of gene and kmer by contigs
    ####

    if args.output:
        if "." in args.output:
            output=args.output
        else:
            output_gene=args.output+"_genesIntoContigs.tsv"
            output_kmer=args.output+"_kmersIntoContigs.tsv"
    else:
        output_gene="GenesIntoContis.tsv"
        output_kmer="KmersIntoContigs.tsv"

    # Conversion on dataframe to dictoniary
    dico_gene_contigs = pd.DataFrame.from_dict(dico_gene_in_rcontigs_count, orient='index')
    dico_kmer_contigs = pd.DataFrame.from_dict(dico_kmer_in_contigs_count, orient='index')

    # Afficher les 5 premières lignes du DataFrame
    print(f"\n\nTSV File {output_gene} Structure :\n\n")
    
    print(dico_gene_contigs.head(5))

    print(f"\n\nTSV File {output_kmer} Structure :\n\n")
    
    print(dico_kmer_contigs.head(5))
    
    
    # save DataFrame on tsv format
    dico_gene_contigs.to_csv(output_gene, sep='\t', index_label='CONTIGS')
    dico_kmer_contigs.to_csv(output_kmer, sep='\t', index_label='CONTIGS')
    print(f"****\n\n{output_gene} and {output_kmer} are sucessly created.\n\n****")
    

    print(f"\nNumber of gene is {len(dico_gene_kmer.keys())}\n")
    print(f"Kmer number a mapq > {mapq} is: {nb_kmer}/{nb_kmer_tot}\n")
    if args.similar: 
        dico_gene = {gene_name: {'contigs':dico_gene_kmer[gene_name]['contigs'],'START': dico_gene_kmer[gene_name]['START'],'STOP':dico_gene_kmer[gene_name]['STOP'],'length': dico_gene_kmer[gene_name]['length'],'nb_kmers': nb_kmers,'Similar':dico_gene_kmer[gene_name]['Similar'],'Note':dico_gene_kmer[gene_name]['Note']} 
                            for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= nb_kmers_filter} 
    else:
        dico_gene = {gene_name: {'contigs':dico_gene_kmer[gene_name]['contigs'],'START': dico_gene_kmer[gene_name]['START'],'STOP':dico_gene_kmer[gene_name]['STOP'],'length': dico_gene_kmer[gene_name]['length'],'nb_kmers': nb_kmers} 
                    for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= nb_kmers_filter}     

    #cleaning memory
    del bed_file
    del bed
    del dico_count_kmer_in_gene
    del dico_gene_kmer
    

    print(f"Genes number containing >= {nb_kmers_filter} kmers and a mapq > {mapq} : {len(dico_gene.keys())}\n\n")    


    if args.output:
        if "." in args.output:
            output=args.output
        else:
            output=args.output+".tsv"
    else:
        output=args.output+"GENE_INTERSECT.tsv"
    
    print(f"****\n\nTSV File {output} Structure :\n\n****")
    count=0
    for key in dico_gene.keys():
        if count==0:
            print("GENE |", end='')
            for value in dico_gene[key]:
                print(f"{value} |", end='')
            print("\n")    
        if count < 5: 
            print(f"{key} | ", end='')
            count+=1
            
            for value in dico_gene[key].values():
                print(f"{value} | ", end='')
            print("\n")
        else:
            break

    # Conversion du dictionnaire en DataFrame
    dico_gene = pd.DataFrame.from_dict(dico_gene, orient='index')

    # Enregistrement du DataFrame en tant que fichier TSV
    dico_gene.to_csv(output, sep='\t', index_label='Gene')
    
    
    print(f"****\n\n{output}.tsv was sucessly created.\n\n****")

else:
    print("\n\n******************\n\nPlease give bed file using -b option\n\n*******************")

