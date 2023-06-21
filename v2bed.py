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
import random

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
parser.add_argument("-s", "--similar", action='store_true',help="Create col Similar if Mentioned in bed file similar code gene")
parser.add_argument("-o", "--output", type=str, help="Name of output file of gene list")
parser.add_argument("-chunk", "--chunk", type=int, help="Nb of kmer DE are selected by group by default 100 000")

args = parser.parse_args()

dico_kmer_pval_log={}

if args.mapq:
    mapq=args.mapq
else:
    mapq=15    

print(f"Le seuil de mapping est {mapq}")

"""
 


"""

if args.diff:
    diff_file=pd.read_csv(args.diff, sep='\t')
    diff=pd.DataFrame(diff_file)
    # Supprimer le fichier de la mémoire
    del diff_file
    data_diff={}
    if args.chunk:
        # Vérifier que le nombre de lignes dans le fichier est suffisant pour la sélection aléatoire
        if len(len(diff)) < args.chunk:
            print('Le fichier ne contient pas suffisamment de lignes pour la sélection aléatoire.')
        else:
            # Sélectionner args.chunk valeurs d'indices uniques aléatoires
            indices_aleatoires = random.sample(range(len(diff)), args.chunk)

            # Sélectionner les lignes correspondant aux indices aléatoires
            lignes_aleatoires = [diff.iloc[i,] for i in indices_aleatoires]

            # Afficher les lignes sélectionnées
            print(f'{args.chunk} lignes sélectionnées :')

        data_diff=lignes_aleatoires   
    else:
        # Sélectionner args.chunk valeurs d'indices uniques aléatoires
        indices_aleatoires = random.sample(range(len(diff)), 100000)
        # Sélectionner les lignes correspondant aux indices aléatoires
        for i in indices_aleatoires:
            data_diff[i]={'tags':diff.iloc[i,0],'pvalue':diff.iloc[i,1], 'log2FC':diff.iloc[i,4], 'Comptage':[diff.iloc[i,5:]]}
        # Afficher les lignes sélectionnées
        print('100 000 lignes sélectionnées :')

    # Ajouter une colonne "Expression" basée sur le seuil de p-value et le logFC
    for i, infos in tqdm(data_diff.items()):
        x=infos['log2FC']
        y=infos['pvalue']
        if x > 0 and y < 0.05:
            infos['Expression'] = 'Up-regulated'
        elif x < 0 and y < 0.05:
            infos['Expression'] = 'Down-regulated'
        else:
            infos['Expression'] = 'Unchanged'
    
    log2FC_values = []
    pvalue_values = []
    expression_values = []

    for key, values in data_diff.items():
        log2FC_values.append(values['log2FC'])
        pvalue_values.append(values['pvalue'])
        expression_values.append(values['Expression'])
        #définis les couleurs pour le volcano plot    
    colors = {'Up-regulated': 'red', 'Down-regulated': 'blue', 'Unchanged': 'grey'}
    premiere_cle = next(iter(data_diff))
    print(data_diff[premiere_cle])    
    # Créer un graphique de type Volcano Plot
    plt.scatter(log2FC_values, -1 * np.log10(pvalue_values), c=[colors[expression] for expression in expression_values])

    # Ajouter les labels des axes et la légende
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 p-value')
    plt.legend(colors)

    # Afficher le scatter plot
    #plt.show()
    # Sauvegader le graphique
    plt.savefig("volcano_diff_counts.pdf") 
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
                if args.similar:
                    if bed.iloc[i,20].strip().split(';')[2].split('=')[0]=="Note":
                        gene_similar = bed.iloc[i,20].strip().split(';')[2].split('=')[1].split(':')[0].split(' ')[2]            
                    else:
                        gene_similar = bed.iloc[i,20].strip().split(';')[3].split('=')[1].split(':')[0].split(' ')[2]
                    dico_gene_kmer[gene_name] = {'kmer': [kmer],'length':length, 'region':region,'START': gene_start, 'STOP': gene_stop, 'Similar':gene_similar}
                else:
                    dico_gene_kmer[gene_name] = {'kmer': [kmer],'length':length, 'region':region, 'START': gene_start, 'STOP': gene_stop}


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

    if args.similar:
        dico_gene = {gene_name: {'Similar':dico_gene_kmer[gene_name]['Similar'],'nb_kmers': nb_kmers,'region':dico_gene_kmer[gene_name]['region'],'START': dico_gene_kmer[gene_name]['START'],'STOP':dico_gene_kmer[gene_name]['STOP'], 'length': dico_gene_kmer[gene_name]['length'], 'kmer': dico_gene_kmer[gene_name]['kmer']} 
                            for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= 10} 
    else:
        dico_gene = {gene_name: {'nb_kmers': nb_kmers,'region':dico_gene_kmer[gene_name]['region'],'START': dico_gene_kmer[gene_name]['START'],'STOP':dico_gene_kmer[gene_name]['STOP'], 'length': dico_gene_kmer[gene_name]['length'], 'kmer': dico_gene_kmer[gene_name]['kmer']} 
                            for gene_name, nb_kmers in dico_count_kmer_in_gene.items() if nb_kmers >= 10}
    
    #Suppression des dictionnaire intermédiaire
    del dico_count_kmer_in_gene
    del dico_gene_kmer

    premiere_cle = next(iter(dico_gene))
    print(dico_gene[premiere_cle])


    print(f"Le nombre de gene avec un nombre de kmer >= 10 et une valeur de mapq > {mapq} est {len(dico_gene)}")
    
    if args.diff:
        print(f"\n\n*******************\n\nLancement de l'étape: \t\t MERGING GENES STAT\n\n*******************")
        for gene_name, infos in tqdm(dico_gene.items()):
            kmers=infos['kmer']
            pvalue=[]
            logFC=[]
            dico_gene[gene_name].pop('kmer')
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
            infos['medianne_Pval']=mediane
            infos['min_Pval']=minimum

        premiere_cle = next(iter(dico_gene))
        print(dico_gene[premiere_cle])

    if args.CF and args.CF2:
        print(f"\n\n*******************\n\nLancement de l'étape: \t\t TAGS GENES WITH CF FILE\n\n*******************")
        comp_file=pd.read_csv(args.CF, sep=',')
        comp = pd.DataFrame(comp_file)

        comp2_file=pd.read_csv(args.CF2, sep=',')
        comp2 = pd.DataFrame(comp2_file)

        nb_match=0
        nb_match1_2=0
        nb_match2=0
        #Suppression de comp_file

        del comp_file
        del comp2_file

        for gene_name, infos in tqdm(dico_gene.items()):
            name=False # Booléen pour true si trouvé avec le nom du gène ou le nom similaire sinon on passe 
            # à la recherche par start et stop du gène en fonction de la région étudié
            for i in range(len(comp)):
                gene_comp=comp.iloc[i,6]
            #for gene_comp in comp.iloc[:,6].values:
                if gene_comp.lower() in gene_name.lower():# == gene_comp:
                    # Ajout du tags 'pi'
                    infos['tags'] = 'Pi'
                    nb_match +=1
                    print("MATCH1") 
                    name=True
                if gene_comp.lower() in infos['Similar'].lower(): # == gene_comp:
                    # Ajout du tags 'pi'
                    infos['tags'] = 'Pi'
                    nb_match +=1
                    print("MATCH2")
                    name=True
                if name==False:
                    #print(f"{infos['region'].lower()} ==> {comp.iloc[i,0].lower()}")
                    if infos['region'].lower() == comp.iloc[i,0].lower():
                        # [-------------------] gene in bed
                        #       ---------       gene in comp file
                        #print("test")
                        if int(infos['START']) < int(comp.iloc[i,4]) and int(infos['STOP']) > int(comp.iloc[i,5]):
                        #    print("OUI I")
                            nb_match +=1
                        # [-------------------]        gene in bed
                        #                   ---------- gene in comp file
                        if int(infos['START']) < int(comp.iloc[i,4]) and int(infos['STOP']) > int(comp.iloc[i,5]) and int(infos['STOP']) > int(comp.iloc[i:4]):
                        #    print('OUI II')
                            nb_match +=1
                        #   [-------------------] gene in bed
                        # --------                gene in file comp
                        if int(infos['START']) > int(comp.iloc[i,4]) and int(infos['STOP']) > int(comp.iloc[i,5]) and int(infos['START']) < int(comp.iloc[i,5]):
                        #    print('OUI III')
                            nb_match +=1
                        #   [-------------------]       gene in bed
                        # -------------------------     gene in comp file
                        if int(infos['START']) < int(comp.iloc[i,4]) and int(infos['STOP']) < int(comp.iloc[i,5]):
                        #    print('OUI IV')
                            nb_match +=1
                    else:
                        # Sinon, initialiser le champ "tags" à une valeur par défaut
                        infos['tags'] = '0'        
            
            name=False # Booléen pour true si trouvé avec le nom du gène ou le nom similaire sinon on passe 
            # à la recherche par start et stop du gène en fonction de la région étudié
            for i in range(len(comp2)):
                gene_comp=comp2.iloc[i,8]
            #for gene_comp in comp.iloc[:,6].values:
                if gene_comp.lower() in gene_name.lower():# == gene_comp:
                    if tags_value == 'Pi':
                        # Ajout du tags 'fst'
                        infos['tags'] = 'Pi,fst'
                        nb_match1_2+=1
                    else:
                        # Ajout du tags 'fst'
                        infos['tags'] = 'fst'
                        nb_match2+=1
                    name=True
                if gene_comp.lower() in infos['Similar'].lower(): # == gene_comp:
                    if infos['tags'] == 'Pi':
                        # Ajout du tags 'fst'
                        infos['tags'] = 'Pi,fst'
                        nb_match1_2+=1
                    else:
                        # Ajout du tags 'fst'
                        infos['tags'] = 'fst'
                        nb_match2+=1
                    name=True
                if name==False:
                    if infos['region'].lower() == comp2.iloc[i,0].lower():
                        # [-------------------] gene in bed
                        #       ---------       gene in comp file
                        
                        if int(infos['START']) < int(comp2.iloc[i,6]) and int(infos['STOP']) > int(comp2.iloc[i,7]):
                            if infos['tags'] == 'Pi':
                                # Ajout du tags 'fst'
                                infos['tags'] = 'Pi,fst'
                                nb_match1_2+=1
                            else:
                                # Ajout du tags 'fst'
                                infos['tags'] = 'fst'
                                nb_match2+=1
                        
                        # [-------------------]        gene in bed
                        #                   ---------- gene in comp file
                        if int(infos['START']) < int(comp2.iloc[i,6]) and int(infos['STOP']) > int(comp2.iloc[i,7]) and int(infos['STOP']) > int(comp2.iloc[i,6]):
                            if infos['tags'] == 'Pi':
                                # Ajout du tags 'fst'
                                infos['tags'] = 'Pi,fst'
                                nb_match1_2+=1
                            else:
                                # Ajout du tags 'fst'
                                infos['tags'] = 'fst'
                                nb_match2+=1
                        #   [-------------------] gene in bed
                        # --------                gene in file comp
                        if int(infos['START']) > int(comp2.iloc[i,6]) and int(infos['STOP']) > int(comp2.iloc[i,7]) and int(infos['START']) < int(comp2.iloc[i,7]):
                            if infos['tags'] == 'Pi':
                                # Ajout du tags 'fst'
                                infos['tags'] = 'Pi,fst'
                                nb_match1_2+=1
                            else:
                                # Ajout du tags 'fst'
                                infos['tags'] = 'fst'
                                nb_match2+=1
                        #   [-------------------]       gene in bed
                        # -------------------------     gene in comp file
                        if int(infos['START']) < int(comp2.iloc[i,6]) and int(infos['STOP']) < int(comp2.iloc[i,7]):
                            if infos['tags'] == 'Pi':
                                # Ajout du tags 'fst'
                                infos['tags'] = 'Pi,fst'
                                nb_match1_2+=1
                            else:
                                # Ajout du tags 'fst'
                                infos['tags'] = 'fst'
                                nb_match2+=1
                    else:
                        # Sinon, initialiser le champ "tags" à une valeur par défaut
                        infos['tags'] = '0'

        print(f"Le nombre de gène qui ont matché avec le fichier{args.CF} est {nb_match}")
        print(f"Le nombre de gène qui ont matché avec le fichier {args.CF2} est {nb_match2}")
        print(f"Le nombre de gène qui ont matché avec les fichiers {args.CF} et {args.CF2} est {nb_match1_2}")
        print("\n\nVoici la structure du dictionnaire final :\n\n")
        premiere_cle = next(iter(dico_gene))
        print(dico_gene[premiere_cle])
    if args.output:
        output=args.output+".tsv"
    else:
        output="output.tsv"
    with open(output, 'w', newline='') as f:
    # Créer un objet writer pour écrire dans le fichier TSV
        writer= csv.writer(f, delimiter='\t')
        for key, value in dico_gene.items():
            writer.writerow([key, value])
else:
    print("\n\n*******************\n\nLancement de l'étape: \t\t NB KMERS TO GENE\n\n*******************")
    print("\n\n******************\n\nFichier BED ABSENT\n\nVeuillez relancer l'analyse en vérifiant que l'arguement -bed est bien renseigner\n\n*******************")


#if args.CF2:
#    
#    file_comp2=pd.read_csv(args.CF2, sep=',')
#    file_comp2 = pd.DataFrame(file_comp2)


if args.gff:

    gff=pd.read_csv(args.gff, sep='\t')
    gff=pd.DataFrame(gff)