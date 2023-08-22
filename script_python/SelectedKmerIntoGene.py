#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @author Djamel Boubred, stage M2 IRD DIADE

"""
Calculate statistics DE from a bedfile after launch this samtools command line
bedtools intersect -abam aln_ref_vs_kmers.bam -b genes.gff -bed -wb -wa > intersection_btw_aligned_kmers_and_genes.bed
and diff file after lauch dekupl pipeline
"""

import collections
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import argparse
from collections import Counter
import statistics as stat
from tqdm import tqdm
import random
import seaborn as sns
#from qmplot import manhattanplot

parser = argparse.ArgumentParser()

parser.add_argument("-diff", "--diff", type=str, help="File with pvlaue and logFoldChange after DE analysis, the file should have pvalue and log2FC name of colums or only pvalue header for DE analysis")
parser.add_argument("-b", "--bed", type=str, help="Bed file from bedtools intersect, the format of tags should be 12354_tags_AATTATTT or AATTAAT")
parser.add_argument("-s", "--similar", action='store_true',help="Create a column called Similar if this word is found in bed file, [default: no activate]")
parser.add_argument("-o", "--output", type=str, help="output file name without extension, [default -o = GENE_INTERSECT]")
parser.add_argument("-mapq", "--mapq", type=int, help="mapq filter values, [default mapq = 15]")
parser.add_argument("-k", "--nbkmer", type=int, help=" numbers of kmer reccurence for concidere a gene, [default k = 10]")
parser.add_argument("-chunk", "--chunk", type=int, help="Nb of kmer DE are selected by group, [default 100 000]")
parser.add_argument("-baseChr", "--baseChr", type=str, help="Given the prefix of Chromosomes names in your samples")
parser.add_argument("-baseUN", "--baseUN", type=str, help="Given the prefix of Unknown Chromosomes names in your samples")

#parser.add_argument("-acp", "--ACP", type=int, help="Activate a ACP visualization if have aboundance in diff file for each samples in the last columns, mentionned the first colums with aboundance\nex: -acp 5 for the colums 6 (first index=0)")

args = parser.parse_args()

dico_kmer_pval_log={}
print("\nVerify the argument:\n")
if args.diff and args.bed:
    
    if args.baseChr and args.baseUN:
        baseChr = args.baseChr
        baseUN = args.baseUN
    else:
        print("baseCHR and baseUN not available, mahattan plot not available")
    if args.output:
        if "." in args.output:
            output=args.output
        else:
            output=args.output+"_GENE_INTERSECT.tsv"
         
    else:
        output="GENE_INTERSECT.tsv"
    if args.mapq:
        mapq=args.mapq
    else:
        mapq=15
    print(f"Mapping quality : {mapq}\n")
    if args.nbkmer:
        nb_kmers_filter=args.nbkmer
    else:
        nb_kmers_filter = 10
    print(f"Numbers ok kmer reccurence : {nb_kmers_filter}") 
    
    #INITIALIZE BALEAN VALUES FOR PRESENCE ABSENCE LOG2FC IN COLUMNS
    col_pvalue=False
    col_log2FC=False
    col_tag=True
    
    bed_file=pd.read_csv(args.bed, sep='\t')
    bed=pd.DataFrame(bed_file)
    
    print(f"\n\n{args.bed} FILE IS LOAD: \n\n")
    
    # Charge un fichier
    diff_file=pd.read_csv(args.diff, sep='\t')
    
    # Charger votre DataFrame à partir d'un fichier
    diff=pd.DataFrame(diff_file)
    
    print(f"\n\n{args.diff} FILE is LOAD: \n\n")
    print(f"\nVerify the file format of {args.diff}:\n")
    # Vérifier si au moins une colonne a pour en-tête "tag"
    if any(colonne == "tag" for colonne in diff.columns):
        print(f"\n\nColums | tag | Checked.\n\n")
        col_tag=True
    else:
        print(f"\n\nColums | tag |.\n\n**/ERREOR STOPPING ANALYZE\**\n\n")
        col_tag=False

    # Vérifier si au moins une colonne a pour en-tête "pvalue"
    if any(colonne == "pvalue" for colonne in diff.columns):
        print(f"\n\nColums | pvalue | Checked.\n\n")
        col_pvalue=True
    else:
        print(f"\n\nColums | pvalue |.\n\n **/ERREOR STOPPING ANALYZE\**\n\n")
        col_pvalue=False 
    
    # Vérifier si au moins une colonne a pour en-tête "log2FC"
    if any(colonne == "log2FC" for colonne in diff.columns):
        print(f"****\n\nColums | log2FC | Checked.\n\n\nVolcano plot in output.\n\n")
        col_log2FC=True
    else:
        print(f"\n\nColums | log2FC |.\n\nVolcano plot not in output.\n\n")
        col_log2FC=False
    
    # Supprimer le fichier de la mémoire
    del diff_file
    data_diff={}  

    if col_pvalue==True and col_tag==True:
        if col_log2FC==True:
            if args.chunk:
                # Vérifier que le nombre de lignes dans le fichier est suffisant pour la sélection aléatoire
                if len(len(diff)) < args.chunk:
                    print('File does not contain enough rows for random selection.')
                else:

                    # Sélectionner args.chunk valeurs d'indices uniques aléatoires
                    indices_aleatoires = random.sample(range(len(diff)), args.chunk)

                    # Sélectionner les lignes correspondant aux indices aléatoires
                    lignes_aleatoires = diff.loc[indices_aleatoires]
    
                    # Afficher les lignes sélectionnées
                    print(f'{args.chunk} lines selected:')
    
                data_diff_random=lignes_aleatoires   
            else:
                # Sélectionner args.chunk valeurs d'indices uniques aléatoires
                indices_aleatoires = random.sample(range(len(diff)), 100000)
                
                # Sélectionner les lignes correspondant aux indices aléatoires
                data_diff_random = diff.loc[indices_aleatoires, ['tag', 'pvalue', 'log2FC']]

                # Afficher les lignes sélectionnées
                print('100 000 lines selected :')
    
            # Ajouter une colonne "Expression" basée sur le seuil de p-value et le log2FC
            for index, row in data_diff_random.iterrows():
                
                tag_value = row['tag']
                pvalue_value = row['pvalue']
                log2FC_value = row['log2FC']

                if log2FC_value > 0 and pvalue_value < 0.05:
                    data_diff_random.at[index, 'Expression'] = 'Up-regulated'
                elif log2FC_value < 0 and pvalue_value < 0.05:
                    data_diff_random.at[index, 'Expression'] = 'Down-regulated'
                else:
                    data_diff_random.at[index, 'Expression'] = 'Unchanged'
            
            #définis les couleurs pour le volcano plot    
            colors = {'Up-regulated': 'red', 'Down-regulated': 'blue', 'Unchanged': 'grey'}

            
            data_diff_random['pvalue'] = np.where(data_diff_random['pvalue'] == 0, 1e-100, data_diff_random['pvalue'])

            # Calculer le logarithme base 10 en évitant la division par zéro
            pvalue_log10 = -1 * np.log10(data_diff_random['pvalue'])

            # Tracer le scatter plot avec les couleurs appropriées
            plt.scatter(data_diff_random['log2FC'], pvalue_log10, c=[colors[expression] for expression in data_diff_random['Expression']])

            # Ajouter les labels des axes et la légende
            plt.xlabel('log2 Fold Change')
            plt.ylabel('-log10 p-value')
            plt.legend(colors)

            # Sauvegader le graphique
            plt.savefig("volcano_diff_counts.pdf") 
        print(f"\n\n*******************\n\nSTEP: \t\t STAT KMER SELECTED FROM DIFF FILE\n\n*******************")

        # Copie du DataFrame avec uniquement les colonnes 'pvalue', 'log2FC' et 'tags'
        if col_log2FC==True: 
            data_diff = diff[['tag','pvalue', 'log2FC']].copy()
        else:
            data_diff = diff[['tag','pvalue']].copy()
        for index, row in tqdm(data_diff.iterrows()):
                
            kmer = row['tag']
            pvalue = row['pvalue']

            if col_log2FC==True:
                log2FC_value = row['log2FC']
                dico_kmer_pval_log[kmer]={'pvalue': pvalue, 'log2FC': log2FC_value}
            else:
                dico_kmer_pval_log[kmer]={'pvalue': pvalue}

        # Supprimer la variable diff
        del diff

        print(f"\n\n*******************\n\nSTEP: \t\t NB KMERS TO GENE\n\n*******************")
        
        dico_count_kmer_in_gene={}
        dico_gene_kmer= {}
        gene_identique=False # booléen dans le cas ou le gene est le même gène que la ligne d'avant on a la valeur True
        nb_kmer=0 #itération qui compte le nombre de kmer total qui sont pris en compte en fonction de nos kmer

        liste_gene_in_region=[]

        for i in tqdm(range(len(bed))):
            if int(bed.iloc[i,4]) > mapq:
                nb_kmer+=1
                region=bed.iloc[i,0]
                print(bed.iloc[i, 20].strip().split(';')[1].split('=')[1])
                gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
                kmer=bed.iloc[i, 3].strip().split('_')[2]
                position_kmer=bed.iloc[i,1]

                if gene_name in dico_gene_kmer:
                    # La clé existe, ajouter le k-mer à la liste existante
                    dico_gene_kmer[gene_name]['kmer'].append(kmer)

                else:
                    # La clé n'existe pas, créer une nouvelle entrée avec la clé et le k-mer
                    gene_start=int(bed.iloc[i, 15])
                    gene_stop=int(bed.iloc[i, 16])
                    length=gene_stop-gene_start

                    if args.similar:
                        if bed.iloc[i,20].strip().split(';')[2].split('=')[0]=="Note":
                            gene_similar = bed.iloc[i,20].strip().split(';')[2].split('=')[1].split(':')[0].split(' ')[2]
                            note_function = bed.iloc[i,20].strip().split(';')[2].split('=')[1]            
                        else:
                            gene_similar = bed.iloc[i,20].strip().split(';')[3].split('=')[1].split(':')[0].split(' ')[2]
                            note_function = bed.iloc[i,20].strip().split(';')[3].split('=')[1]

                        dico_gene_kmer[gene_name] = {'kmer': [kmer], 'region':region,'START': gene_start, 'STOP': gene_stop,'length':length,'Similar':gene_similar,'Note':note_function}
                    else:
                        dico_gene_kmer[gene_name] = {'kmer': [kmer], 'region':region, 'START': gene_start, 'STOP': gene_stop,'length':length}

                if gene_name in dico_gene_kmer:
                    dico_count_kmer_in_gene[gene_name] = len(dico_gene_kmer[gene_name]['kmer'])
                else:
                    print(f"La clé '{gene_name}' n'existe pas dans le dictionnaire 'dico_gene_kmer'. Sa valeur de mapQ est = '{int(bed.iloc[i,4])}")

        nb_kmer_tot=len(bed)

        print(f"\nNumber of gene with kmer mapq > {mapq} is {len(dico_gene_kmer.keys())}\n")
        print(f"Kmer number a mapq > {mapq} is: {nb_kmer}/{nb_kmer_tot}\n")

        if args.similar:
            dico_gene = {gene_name: {'region':dico_gene_kmer[gene_name]['region'],'START': dico_gene_kmer[gene_name]['START'],'STOP':dico_gene_kmer[gene_name]['STOP'], 'length': dico_gene_kmer[gene_name]['length'],'nb_kmers': nb_kmers,'Similar':dico_gene_kmer[gene_name]['Similar'], 'Note':dico_gene_kmer[gene_name]['Note'],'kmer': dico_gene_kmer[gene_name]['kmer']} 
                                for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= nb_kmers_filter} 
        else:
            dico_gene = {gene_name: {'region':dico_gene_kmer[gene_name]['region'],'START': dico_gene_kmer[gene_name]['START'],'STOP':dico_gene_kmer[gene_name]['STOP'], 'length': dico_gene_kmer[gene_name]['length'],'nb_kmers': nb_kmers,'kmer': dico_gene_kmer[gene_name]['kmer']} 
                                for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= nb_kmers_filter}
        
        
        #cleaning memory
        del bed_file
        del bed
        del dico_count_kmer_in_gene
        del dico_gene_kmer

        print(f"Genes number containing >= {nb_kmers_filter} kmers and a mapq > {mapq} : {len(dico_gene.keys())}\n")    
        
        print(f"\n\n*******************\n\nSTEP: \t\t MERGING GENES STAT\n\n*******************")
        
        gene_artefact=[]
        for gene_name, infos in tqdm(dico_gene.items()):
            kmers=infos['kmer']
            pvalue=[]
            if args.baseChr and args.baseUN:
                if baseUN in infos['region']:
                    infos['Classe']= 'PAN'
                else:
                    if baseChr in infos['region']:
                        infos['Classe'] = infos['region']
                    else:
                        infos['Classe'] = 'PAN'
            dico_gene[gene_name].pop('kmer')
            if col_log2FC==True:
                log2FC=[]
            for kmer in kmers:
                if kmer in dico_kmer_pval_log:  
                    pvalue.append(dico_kmer_pval_log[kmer]['pvalue'])
                    if col_log2FC==True:
                        log2FC.append(dico_kmer_pval_log[kmer]['log2FC'])
            if len(pvalue)<nb_kmers_filter:
                gene_artefact.append(gene_name)
            else:
                mediane = stat.median(pvalue)
                minimum = min(pvalue)
                infos['medianne_Pval']=mediane
                infos['min_Pval']=minimum
            if col_log2FC==True:
                if len(log2FC) < nb_kmers_filter:
                    gene_artefact.append(gene_name)
                else:
                    moyenne = stat.mean(log2FC)
                    ecart_type = stat.stdev(log2FC)
                    infos['mean_Log2FC']=moyenne
                    infos['sd_Log2FC']=ecart_type

        for gene_name in list(set(gene_artefact)):
            del dico_gene[gene_name]
        data=pd.DataFrame.from_dict(dico_gene, orient='index')

        print(f"Genes number containing >= {nb_kmers_filter} kmers and a mapq > {mapq} after merging data : {len(dico_gene.keys())}\n") 
        # Liste de régions avec des duplications
        if args.baseChr and args.baseUN:
            # Crée une copie du DataFrame en filtrant les lignes
            data = data[data['Classe'] != 'PAN'].copy()
            xtick=[]
            
            ### sort par col ref puis par position
            data = data.sort_values(by = ['Classe', 'START'], ascending = [True, True], na_position = 'first')
            data['-logp']= - np.log10(data['medianne_Pval'])
            run_pos = 0
            cum_pos = []
            
            # calcule cumulative positions to give a x axe size in manhattan
            for chrom, group_df in data.groupby('Classe'):
                cum_pos.append(group_df['START'] + run_pos)
                run_pos += group_df['START'].max()
            
            # new columns
            data['cum_pos'] = pd.concat(cum_pos)
            data['Genenb'] = data.index
            
            #Manhattan plot

            g = sns.relplot(data=data, x='cum_pos', y='-logp', aspect=4, hue="Classe", palette='Set1', linewidth=0, s=8, legend= None )
            g.ax.set_xlabel('Chromosome')
            g.ax.set_ylabel('-logp')
            g.ax.set_xticks(data.groupby('Classe')['cum_pos'].median())
            g.ax.set_xticklabels(data['Classe'].unique())
            # Rotation des étiquettes de l'axe x
            g.ax.tick_params(axis='x', rotation=45)
            # Ajustement des marges pour faire de la place à l'axe des y
            plt.subplots_adjust(right=0.1)  # Augmentez la valeur si nécessaire
            g.ax.axhline(3, linestyle="--", linewidth=1)
            g.fig.suptitle("Manhattan plot")
            # Ajout de tight_layout pour ajuster les marges
            plt.tight_layout()
            plt.show()

        if args.similar:
            # Définir l'ordre souhaité des clés à l'intérieur de chaque dictionnaire
            ordre_des_cles = ['region','START','STOP','length','nb_kmers','medianne_Pval','min_Pval','mean_Log2FC','sd_Log2FC','Similar','Note']
        else:
             # Définir l'ordre souhaité des clés à l'intérieur de chaque dictionnaire
            ordre_des_cles = ['region','START','STOP','length','nb_kmers','medianne_Pval','min_Pval','mean_Log2FC','sd_Log2FC']    
        # Créer un nouvel Ordered Dictionary avec les clés dans l'ordre souhaité
        
        dico_gene_order = {gene: collections.OrderedDict((cle, dico_gene[gene][cle]) for cle in ordre_des_cles) for gene in dico_gene}
        dico_gene=dico_gene_order


        print("\n\n\n")


        print(f"\n\nTSV File {output} Structure (5 first lines):\n\n")
        count=0
        for key in dico_gene.keys():
            if count==0:
                print("GENE | ", end='')
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
            
            
        print(f"\n\n{output} was sucessly created.\n\n")

    else:
        if col_pvalue==False:
            print("\n\nColums | pvalue | not checked in, Verify your input file.\n The header name could be 'pvalue'\n\n**/ERREOR STOPPING ANALYZE\**\n\n")
            print(" (°~/°) ")
        if col_tag==False:
            print("\n\nColums | tag | not checked in, Verfy your input file.\n The header name could be 'tag'\n\n**/ERREOR STOPPING ANALYZE\**\n\n")
            print(" (°~/°) ")
else:
    if args.diff:
        print(f"\n\n\n\n Diff argument not available, Verify your path of input file or argument -diff \n\n**/ERREOR STOPPING ANALYZE\**\n\n")
    if args.bed:
        print(f"\n\nBed argument not available, Verify your path of input file or argument -b \n\n**/ERREOR STOPPING ANALYZE\**\n\n")