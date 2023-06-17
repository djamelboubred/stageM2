import subprocess
from subprocess import PIPE
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import argparse
from collections import Counter
import statistics as stat
from tqdm import tqdm

#conda activate myenv
#python v2bed.py -diff diff-counts.tsv -index tag pvalue log2FC -b kmer_60k.bed -gff gene_60k.gff -CF liste_gene_pi_analisis.csv -CF2 liste_gene_fst_analisis.csv

parser = argparse.ArgumentParser()


parser.add_argument("-diff", "--diff", type=str, help="File with pvlaue and logFoldChange after DE analysis")
parser.add_argument('-index', nargs='+', help='List of indices') #prend en argument une liste de nom de colonnes en premier celle avec les kmer en second la pvalue et en dernier le logFoldChange
parser.add_argument("-b", "--bed", type=str, help="Bed file")
parser.add_argument("-gff", "--gff", type=str, help="GFF file format with gene anotation")
parser.add_argument("-CF", "--CF", type=str, help="Comparative file with gene see by different analysis")
parser.add_argument("-CF2", "--CF2", type=str, help="Comparative file with gene see by different analysis")
parser.add_argument("-mapq", "--mapq", type=int, help="mapq filter values by defaults mapq = 15")
args = parser.parse_args()

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
    bed_file=pd.read_csv(args.bed, sep='\t')
    bed=pd.DataFrame(bed_file)
    print(f"\n\n*******************\n\nLancement de l'étape: \t\t NB KMERS TO GENE\n\nLe fichier {args.bed} est chargé\n\n*******************")

    dico_count_kmer_in_gene={}
    dico_gene_kmer= {}
    gene_identique=False # booléen dans le cas ou le gene est le même gène que la ligne d'avant on a la valeur True
    nb_kmer=0 #itération qui compte le nombre de kmer total qui sont pris en compte en fonction de nos kmer

    for i in tqdm(range(len(bed))):
        if int(bed.iloc[i,4]) > mapq:
            nb_kmer+=1
            region=bed.iloc[i,0]
            gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
            kmer=bed.iloc[i, 3].strip().split('_')[2]
            if gene_name in dico_gene_kmer:
                # La clé existe, ajouter le k-mer à la liste existante
                dico_gene_kmer[gene_name]['kmer'].append(kmer)
                
            else:
                # La clé n'existe pas, créer une nouvelle entrée avec la clé et le k-mer
                gene_start=int(bed.iloc[i, 15])
                gene_stop=int(bed.iloc[i, 16])
                length=gene_stop-gene_start
                if region[0]!='C' and region[1]!='h' and region[2]!='r':
                    region='Pan'           
                             
                dico_gene_kmer[gene_name] = {'kmer': [kmer],'length':length, 'region':region}


            #
            #Vérification
            #
            
            if gene_name in dico_gene_kmer:
                dico_count_kmer_in_gene[gene_name] = len(dico_gene_kmer[gene_name]['kmer'])
            else:
                print(f"La clé '{gene_name}' n'existe pas dans le dictionnaire 'dico_gene_kmer'. Sa valeur de mapQ est = '{int(bed.iloc[i,4])}")


    premiere_cle = next(iter(dico_gene_kmer))
    print(dico_gene_kmer[premiere_cle])

    #Visualisation de la première clé du dico_count_kmer_in_gene
    premiere_cle = next(iter(dico_count_kmer_in_gene))
    print(dico_count_kmer_in_gene[premiere_cle])

    nb_kmer_tot=len(bed)
    print(f"Le Nombre de gène dans le fichier {args.bed} est {len(dico_gene_kmer.keys())}")
    print(f"Le nombre de kmer ayant une valeur de mapq > {mapq} est {nb_kmer}/{nb_kmer_tot}")

    #Suppression du fichier bed en mémoire
    del bed_file

    #Suppression de la variables bed
    del bed


    dico_gene = {gene_name: {'nb_kmers': nb_kmers,'region':dico_gene_kmer[gene_name]['region'], 'length': dico_gene_kmer[gene_name]['length'], 'kmer': dico_gene_kmer[gene_name]['kmer']} 
                        for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= 10}
    
    #Suppression des dictionnaire intermédiaire
    del dico_count_kmer_in_gene
    del dico_gene_kmer

    premiere_cle = next(iter(dico_gene))
    print(dico_gene[premiere_cle])


    print(f"Le nombre de gene avec un nombre de kmer >= 10 et une valeur de mapq > {mapq} est {len(dico_gene)}")

    #supprimer le dictionnaire dico_count_kmer_in_gene

    #del dico_count_kmer_in_gene
    
    if args.diff:
        print(f"\n\n*******************\n\nLancement de l'étape: \t\t MERGING GENES STAT\n\n*******************")
        for gene_name, infos in tqdm(dico_gene.items()):
            kmers=infos['kmer']
            pvalue=[]
            logFC=[]
            for kmer in kmers:
                if kmer in dico_kmer_pval_log:
                    pvalue.append(dico_kmer_pval_log[kmer]['pvalue'])
                    logFC.append(dico_kmer_pval_log[kmer]['logFC'])
                else:
                    print(f"ERROR KEY.VALUES {kmer}")
            moyenne = stat.mean(logFC)
            ecart_type = stat.stdev(logFC)
            mediane = stat.median(pvalue)
            minimum = min(pvalue)
            infos['mean_LogFC']=moyenne
            infos['sd_LogFC']=ecart_type
            infos['medianne']=mediane
            infos['min_pvalue']=minimum

        del dico_gene['kmer']

        premiere_cle = next(iter(dico_gene_kmer))
        print(dico_gene_kmer[premiere_cle])

    if args.CF:

    file_comp=pd.read_csv(args.CF, sep=',')
    file_comp = pd.DataFrame(file_comp)

    

else:
    print("\n\n*******************\n\nLancement de l'étape: \t\t NB KMERS TO GENE\n\n*******************")
    print("\n\n******************\n\nFichier BED ABSENT\n\nVeuillez relancer l'analyse en vérifiant que l'arguement -bed est bien renseigner\n\n*******************")


    

if args.CF2:
    
    file_comp2=pd.read_csv(args.CF2, sep=',')
    file_comp2 = pd.DataFrame(file_comp2)


if args.gff:

    gff=pd.read_csv(args.gff, sep='\t')
    gff=pd.DataFrame(gff)