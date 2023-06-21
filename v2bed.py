import subprocess
from subprocess import PIPE
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import argparse
from collections import Counter
import statistics
from tqdm import tqdm
import time




gff=pd.read_csv(args.gff, sep='\t')
gff=pd.DataFrame(gff)

dico_kmer_pval_log={}

if args.mapq:
    mapq=args.mapq
else:
    mapq=15    

print(f"Le seuil de mapping est {mapq}")

if args.diff:
    diff_file=pd.read_csv(args.diff, sep='\t')
    diff=pd.DataFrame(diff_file)
    # Supprimer le fichier de la mémoire
    del diff_file

    print(f"\n\n*******************\n\nLancement de l'étape: \t\t PVALUE AND logFC\n\nLe fichier {args.diff} est chargé\n\n*******************")
    if args.index:
        # Affichage du premier élément de la liste
        print("Premier élément de la liste :", args.index[0])
        print("Deuxième élément de la liste :", args.index[1])
        print("Troisième élément de la liste :", args.index[2])
        for i in tqdm(range(len(diff))):
            kmer=diff.loc[i, args.index[0]]
            pvalue=diff.loc[i, args.index[1]]
            logFC=diff.loc[i, args.index[2]]
            dico_kmer_pval_log[kmer]={'pvalue': pvalue, 'logFC': logFC}
    else:
        print("Aucune liste d'indices fournie.")
        for i in tqdm(range(len(diff))):
            kmer=diff.loc[i, 'tag']
            pvalue=diff.loc[i, 'pvalue']
            logFC=diff.loc[i, 'log2FC']
            dico_kmer_pval_log[kmer]={'pvalue': pvalue, 'logFC': logFC}

    #Visualisation de la première clé du dico_kmer_pval_log
    premiere_cle = next(iter(dico_kmer_pval_log))
    print(dico_kmer_pval_log[premiere_cle])
    
    # Supprimer la variable diff
    del diff

if args.bed:
    bed=pd.read_csv(args.bed, sep='\t')
    bed=pd.DataFrame(bed)
    print(f"\n\n*******************\n\nLancement de l'étape: \t\t NB KMERS TO GENE\n\nLe fichier {args.bed} est chargé\n\n*******************")
    dico_read_length_gene={}
    dico_coverage={}
    dico_count_gene={}
    gene_name='first'
    thresholds = [5, 10, 25, 50, 75, 90]
    coverage_threshold = [0] * len(thresholds)
    dico_gene_kmer= {}
    gene_identique=False # booléen dans le cas ou le gene est le même gène que la ligne d'avant on a la valeur True
    nb_kmer=0 #itération qui compte le nombre de kmer total qui sont pris en compte en fonction de nos kmer
    for i in tqdm(range(len(bed))):
        if int(bed.iloc[i,4]) > mapq:
            nb_kmer+=1
            gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
            kmer=bed.iloc[i, 3].strip().split('_')[2]
            if gene_name in dico_gene_kmer:
                # La clé existe, ajouter le k-mer à la liste existante
                dico_gene_kmer[gene_name]['kmer'].append(kmer)
                kmer_start=int(bed.iloc[i, 1])
                kmer_stop=int(bed.iloc[i, 2])

                #Compare les valeur de la ligne (kmer) avec celle déjà enregistrer dans le dico (pour un gene)
                if kmer_start < dico_read_length_gene[gene_name]['kmer_start']:
                    dico_read_length_gene[gene_name]['kmer_start'] = kmer_start
                if kmer_stop > dico_read_length_gene[gene_name]['kmer_stop']:
                    dico_read_length_gene[gene_name]['kmer_stop'] = kmer_stop
                if kmer_start > dico_read_length_gene[gene_name]['kmer_stop']:
                    correction= kmer_start - dico_read_length_gene[gene_name]['kmer-stop']
                    if correction < dico_read_length_gene[gene_name]['correction']:
                        dico_read_length_gene[gene_name]['correction']=correction
            else:
                # La clé n'existe pas, créer une nouvelle entrée avec la clé et le k-mer
                dico_gene_kmer[gene_name] = {'kmer': [kmer]}
                #Analyse du coverage 
                kmer_start=int(bed.iloc[i, 1])
                kmer_stop=int(bed.iloc[i, 2])
                gene_start=int(bed.iloc[i, 15])
                gene_stop=int(bed.iloc[i, 16])
                correction=0
                dico_read_length_gene[gene_name]={'kmer_start':kmer_start, 'kmer_stop':kmer_stop,'gene_start':gene_start, 'gene_stop':gene_stop, 'correction':correction}

            #
            #Vérification
            #
            
            if gene_name in dico_gene_kmer:
                dico_count_gene[gene_name] = {'nb_kmers': len(dico_gene_kmer[gene_name])}
            else:
                print(f"La clé '{gene_name}' n'existe pas dans le dictionnaire 'dico_gene_kmer'. Sa valeur de mapQ est = '{int(bed.iloc[i,4])}")


    premiere_cle = next(iter(dico_gene_kmer))
    print(dico_gene_kmer[premiere_cle])

    print(f"Le Nombre de gène dans le fichier {args.bed} est {len(dico_gene_kmer.keys())}")
    print(f"Le nombre de kmer ayant une valeur de mapq > {mapq} est {nb_kmer}/{len(bed)}")

else:
    print("\n\n*******************\n\nLancement de l'étape: \t\t NB KMERS TO GENE\n\n*******************")
    print("\n\n******************\n\nFichier BED ABSENT\n\nVeuillez relancer l'analyse en vérifiant que l'arguement -bed est bien renseigner\n\n*******************")
