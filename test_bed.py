import subprocess
from subprocess import PIPE
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import argparse
from collections import Counter
import statistics




parser = argparse.ArgumentParser()


parser.add_argument("-diff", "--diff", type=str, help="File with pvlaue and logFoldChange after DE analysis")
parser.add_argument('-index', nargs='+', help='List of indices') #prend en argument une liste de nom de colonnes en premier celle avec les kmer en second la pvalue et en dernier le logFoldChange
parser.add_argument("-b", "--bed", type=str, help="Bed file")
parser.add_argument("-gff", "--gff", type=str, help="GFF file format with gene anotation")
parser.add_argument("-CF", "--CF", type=str, help="Comparative file with gene see by different analysis")
parser.add_argument("-CF2", "--CF2", type=str, help="Comparative file with gene see by different analysis")
args = parser.parse_args()

diff=pd.read_csv(args.diff, sep='\t')
diff=pd.DataFrame(diff)

file_comp=pd.read_csv(args.CF, sep=',')
file_comp = pd.DataFrame(file_comp)

file_comp2=pd.read_csv(args.CF2, sep=',')
file_comp2 = pd.DataFrame(file_comp2)

bed=pd.read_csv(args.bed, sep='\t')
bed=pd.DataFrame(bed)

gff=pd.read_csv(args.gff, sep='\t')
gff=pd.DataFrame(gff)

dico_kmer_pval_log={}

if args.index:
    # Affichage du premier élément de la liste
    print("Premier élément de la liste :", args.index[0])
    print("Deuxième élément de la liste :", args.index[1])
    print("Troisième élément de la liste :", args.index[2])
    for i in range(len(diff)):
        dico_kmer_pval_log[diff.loc[i, args.index[0]]]=[diff.loc[i, args.index[1]],diff.loc[i, args.index[2]]]
        # Récupérer la valeur de la p-value (5ème colonne) et la pvalue (2ème colonne
else:
    print("Aucune liste d'indices fournie.")
    for i in range(len(diff)):
        dico_kmer_pval_log[diff.loc[i, 'tag']]=[diff.loc[i, 'pvalue'],diff.loc[i, 'log2FC']]
        # Récupérer la valeur de la p-value (5ème colonne) et la pvalue (2ème colonne

premiere_cle = next(iter(dico_kmer_pval_log))
# Supprimer le premier élément
print(dico_kmer_pval_log[premiere_cle])

dico_read_length_gene={}
dico_coverage={}
dico_count_gene={}
gene_name='first'
thresholds = [5, 10, 25, 50, 75, 90]
coverage_threshold = [0] * len(thresholds)
dico_gene_kmer= {}
gene_identique=False # booléen dans le cas ou le gene est le même gène que la ligne d'avant on a la valeur True
for i in range(len(bed)):
    if int(bed.iloc[i,4]) > 15:
       #initie le premier gène du fichier ne créant mes variable de position pour le kmer et le gène        
        if gene_identique==False:
            gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
            kmer_start=int(bed.iloc[i, 1])
            kmer_stop=int(bed.iloc[i, 2])
            gene_start=int(bed.iloc[i, 15])
            gene_stop=int(bed.iloc[i, 16])
            correction=0
            gene_identique=True
       
        #vérifie qu'on regarde le même gène (on regarde un kmer du gène)
        # Ajouter le premier k-mer
        if gene_name in dico_gene_kmer:
            # La clé existe, ajouter le k-mer à la liste existante
            dico_gene_kmer[gene_name]['kmer'].append(bed.iloc[i, 3].strip().split('_')[2])
        else:
            print("ABSENT")
            # La clé n'existe pas, créer une nouvelle entrée avec la clé et le k-mer
            dico_gene_kmer[gene_name] = {'kmer': [bed.iloc[i, 3].strip().split('_')[2]]}
        if gene_name==bed.iloc[i, 20].strip().split(';')[1].split('=')[1]:
            if kmer_start > int(bed.iloc[i, 1]):
                kmer_start=int(bed.iloc[i, 1])
            if kmer_stop < int(bed.iloc[i, 2]):
                kmer_stop=int(bed.iloc[i, 2])
        
        #condition dans le cas notre kmer appartient à un gène différent, on calcul le coverage pour le gène d'avant et on initie les valeurs pour le gène actuel
        if gene_name != bed.iloc[i, 20].strip().split(';')[1].split('=')[1]:

            # Vérification et mise à jour du dictionnaire
            if gene_name in dico_read_length_gene:
                #vérifie que le kmer_start est plus petit que celui déjà enregistrer et donc on à une augmentation du reads par la droite
                if kmer_start < dico_read_length_gene[gene_name]['kmer_start']:
                    dico_read_length_gene[gene_name]['kmer_start'] = kmer_start
                #vérifie que le kmer_start est plus petit que celui déjà enregistrer et donc on à une augmentation du reads par la gauche
                if kmer_stop > dico_read_length_gene[gene_name]['kmer_stop']:
                    dico_read_length_gene[gene_name]['kmer_stop'] = kmer_stop
                if kmer_start > dico_read_length_gene[gene_name]['kmer_stop']:
                    correction= kmer_start - dico_read_length_gene[gene_name]['kmer-stop']
                    if correction < dico_read_length_gene[gene_name]['correction']:
                        dico_read_length_gene[gene_name]['correction']=correction
            else:
                dico_read_length_gene[gene_name]={'kmer_start':kmer_start, 'kmer_stop':kmer_stop,'gene_start':gene_start, 'gene_stop':gene_stop, 'correction':correction}
                if gene_name == "maker-MH_428188-snap-gene-0.0":
                    print("Le gene name maker-MH_428188-snap-gene-0.0 est trouvé")
            gene_name=bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
            kmer_start=int(bed.iloc[i, 1])
            kmer_stop=int(bed.iloc[i, 2])
            gene_start=int(bed.iloc[i, 15])
            gene_stop=int(bed.iloc[i, 16])
            correction=0
        #print(dico_gene_kmer)
        ####
####  L'erreur est à ce niveau j'ai mon dico_gene_kmer qui n'a pas les clé gene_name alors qu'elle ont une valeur que je dois concidérer
        ####
        
        if gene_name in dico_gene_kmer:
            dico_count_gene[gene_name] = {'nb_kmers': len(dico_gene_kmer[gene_name])}
        else:
            print(f"La clé '{gene_name}' n'existe pas dans le dictionnaire 'dico_gene_kmer'. Sa valeur de mapQ est = '{int(bed.iloc[i,4])}")

        #print(dico_gene_kmer['Oglab_000002'])
        #dico_count_gene[gene_name] = {'nb_kmers':len(dico_gene_kmer[gene_name])}

print(dico_count_gene)
for gene_name in dico_read_length_gene:
    # Accéder aux informations spécifiques du gène
    gene_info = dico_read_length_gene[gene_name]
    
    # Accéder à chaque valeur du gène
    kmer_start = gene_info['kmer_start']
    kmer_stop = gene_info['kmer_stop']
    gene_start = gene_info['gene_start']
    gene_stop = gene_info['gene_stop']
    correction=gene_info['correction']
    #   |---------------|   gene
    #       |_______|   kmer
    if kmer_start >= gene_start and kmer_stop <= gene_stop:
        reads=kmer_stop-kmer_start
    #   |---------------|   gene
    #              |_______|    kmer
    if kmer_start >= gene_start and kmer_stop >= gene_stop:
        reads=gene_stop-kmer_start
    #   |---------------|   gene
    # |_______| kmer
    if kmer_start <= gene_start and kmer_stop <= gene_stop:
        reads=kmer_stop-gene_start
    #           |-----| gene
    #       |____________|  kmer
    if kmer_start <= gene_start and kmer_stop >= gene_stop:
        reads=gene_stop-gene_start


    length_gene=gene_stop-gene_start
    reads=reads-correction
    coverage=(reads/length_gene)*100
    #print(coverage, type(coverage))
    for i in range(len(thresholds)):
        #print(i)
        if coverage > thresholds[i]:
            coverage_threshold[i]+=1

    dico_coverage[gene_name] = {'coverage': coverage, 'length': length_gene}  # Ajouter les données au dictionnaire

print("Le sueil de couverture est =",thresholds)
print("Le nombre de gène ayant avec une couverture =", coverage_threshold)
#print(dico_gene_kmer["maker-MH_428188-snap-gene-0.0"])
#print(dico_read_length_gene["maker-MH_428188-snap-gene-0.0"])
#print(dico_coverage["maker-MH_428188-snap-gene-0.0"])
dico_fusion = {gene: {'nb_kmers': nb_kmers, 'coverage': dico_coverage[gene]['coverage'], 'length': dico_coverage[gene]['length']} 
                        for gene, nb_kmers in dico_count_gene.items() if nb_kmers >= 10}



dico_gene_log2FC={}

for gene_name in dico_gene_kmer:
    liste_pvalue_gene=[]
    
    for kmer in dico_gene_kmer[gene_name]:
        liste_pvalue_gene.append(dico_kmer_pval_log[kmer][1])
        
    #print(liste_pvalue_gene)
    print(len(liste_pvalue_gene))
    if len(liste_pvalue_gene) >= 10:
        moyenne = statistics.mean(liste_pvalue_gene)
        ecart_type = statistics.stdev(liste_pvalue_gene)
        dico_fusion[gene_name]["MeanLogFC"] = moyenne
        dico_fusion[gene_name]["SDLogFC"] = ecart_type


# Parcours du dictionnaire dico_fusion
nb_match=0
nb_match2=0
nb_match1_2=0
for gene_name, infos in dico_fusion.items():
    # Vérification si gene_name est présent dans la 5ème colonne du dataframe
    gene_name_without_loc = gene_name.replace('LOC_', '')
    if gene_name_without_loc in file_comp.iloc[:, 6].values:
        # Ajout du tags 'pi'
        infos['tags'] = 'Pi'
        nb_match +=1
    if gene_name_without_loc in file_comp2.iloc[:, 8].values:
        # Récupération de la valeur de la clé 'tags' avec la méthode get()
        tags_value = infos.get('tags')
        # Vérification si la valeur de 'tags' est 'pi'
        nb_match2 +=1
        if tags_value == 'Pi':
            # Ajout du tags 'fst'
            infos['tags'] = 'Pi,fst'
            nb_match1_2+=1
        else:
            # Ajout du tags 'fst'
            infos['tags'] = 'fst'     
    else:
        # Sinon, initialiser le champ "tags" à une valeur par défaut
        infos['tags'] = '0'

#print(dico_fusion)
outfile='gene_DE_mapQ15_kmer10.tsv'

print("la taille du dico_coverage =", len(dico_coverage))
print("La taille du dico_counts_gene=", len(dico_count_gene))

print("La taille du dico_fusion (seuil map15 et kmer > 10)=", len(dico_fusion))

print("Le nombre de gène qui concorde avec l'analyse Pi est de ", nb_match)
print("La taille du fichier", args.CF ,"est",len(file_comp))

print("Le nombre de gène qui concorde avec",args.CF2 ,"est de ", nb_match2)
print("La taille du fichier",args.CF2, "est",len(file_comp2))

print("Le nombre de gène qui concorde avec",args.CF,"et", args.CF2,"est de ", nb_match1_2)
print("la taille du dico_coverage =", len(dico_fusion))


# Écriture du dictionnaire dans un fichier TSV
with open(outfile, 'w', newline='') as fichier:
    writer = csv.writer(fichier, delimiter='\t')
    writer.writerow(['Gene', 'Nb_kmers', 'Coverage', 'Length', 'Tags'])  # Écriture de l'en-tête
    for gene, infos in dico_fusion.items():
        nb_kmers = infos.get('nb_kmers', 0)  # Récupération de 'nb_kmers' avec une valeur par défaut de 0
        coverage = infos.get('coverage', 0)  # Récupération de 'coverage' avec une valeur par défaut de 0
        length = infos.get('length', 0)  # Récupération de 'length' avec une valeur par défaut de 0
        tags = infos.get('tags',0) #Récupération de 'tags' avec une valeurs par défault de 0
        writer.writerow([gene, nb_kmers, coverage, length, tags])


def nb_gene_in_gff(gff):
    nb_gene_chr=[]
    name=[]
    for i in range(1,14):
        if i < 10:
            chr="Chr0"+str(i)
        if i > 9:
            chr="Chr"+str(i)
        if i == 13:
            chr="Chr"
            name.append("Pangenome")
        if chr != "Chr":
            name.append(chr) 
        liste_gene_chr=[]    
        print(chr)
        for i in range(len(gff)):
            if chr!="Chr":
                if gff.iloc[i,0] == chr:
                    liste_gene_chr.append(gff.iloc[i, 8].strip().split(';')[1].split('=')[1]) 
            if chr=="Chr":
                if gff.iloc[i,0][0] != chr[0] and gff.iloc[i,0][1] != chr[1] and gff.iloc[i,0][2] != chr[2]:
                    liste_gene_chr.append(gff.iloc[i, 8].strip().split(';')[1].split('=')[1])  
        nb_gene_chr.append(len(set(liste_gene_chr)))

    print("Le nombre de gène contenu dans gff=",nb_gene_chr)

    width=-0.4
    x = list(range(len(name)-1))
    y = nb_gene_chr[:-1]
    last_x = [len(name)-1]
    last_y = nb_gene_chr[-1:]

    # Définir la couleur de l'histogramme
    colors = ['orange'] * (len(name)-1) + ['r']
    # Créer un histogramme avec deux couleurs différentes
    plt.bar(name, nb_gene_chr, width=width,color=colors, label="GFF",align='edge')

def info_gene_kmer_in_genome(bed):
    nb_gene_chr=[]
    nb_kmer_chr=[]
    liste_gene_chr_with_10_kmer_min=[]
    liste_gene_chr_with_50_kmer_min=[]
    liste_gene_chr_with_75_kmer_min=[]
    name=[]
    for i in range(1,14):
        if i < 10:
            chr="Chr0"+str(i)
        if i > 9:
            chr="Chr"+str(i)
        if i == 13:
            chr="Chr"
            name.append("Pangenome")
        if chr != "Chr":
            name.append(chr) 
        liste_gene_chr=[]
        print(chr)
        for i in range(len(bed)):
            #compteur=0
            if chr!="Chr":
                if bed.iloc[i,0] == chr:
                    liste_gene_chr.append(bed.iloc[i, 20].strip().split(';')[1].split('=')[1]) 
            if chr=="Chr":
                if bed.iloc[i,0][0] != chr[0] and bed.iloc[i,0][1] != chr[1] and bed.iloc[i,0][2] != chr[2]:
                    liste_gene_chr.append(bed.iloc[i, 20].strip().split(';')[1].split('=')[1])  
        nb_gene_chr.append(len(set(liste_gene_chr)))
        
        nb_kmer_chr.append(len(liste_gene_chr))

    print("Le nombre de gene est =",nb_gene_chr)
    print("Le noùbre de kmer est =",nb_kmer_chr)


    plt.figure()
    # Définition de la largeur des barres
    width = 0.4
    x = list(range(len(name)-1))
    y = nb_gene_chr[:-1]
    last_x = [len(name)-1]
    last_y = nb_gene_chr[-1:]
    nb_gene_in_gff(gff)
    # Définir la couleur de l'histogramme
    colors = ['b'] * (len(name)-1) + ['green']
    
    print(name)

    # Créer un histogramme avec deux couleurs différentes
    plt.bar(name, nb_gene_chr,width=width ,color=colors, label="BED",align='edge')
    plt.xlabel("Chromosomes")
    plt.ylabel("Nb Gene")
    title='Distribution du nombre de gene pour chaque chromoses et dans le pangenome'
    plt.title(title)
    plt.legend(loc="upper right")
    #plt.savefig("nb_gene_chr.pdf")
    plt.show()
    total=0
    for i in range(12):
        total += nb_kmer_chr[i]
    
    plt.figure()
    x = list(range(len(name)-1))
    y = nb_gene_chr[:-1]
    last_x = [len(name)-1]
    last_y = nb_gene_chr[-1:]

    # Définir la couleur de l'histogramme
    colors = ['b'] * (len(name)-1) + ['r']

    # Créer un histogramme avec deux couleurs différentes
    plt.bar(name, nb_kmer_chr, color=colors)
    plt.xlabel("Localisation")
    plt.ylabel("Nb kmer")
    title='Distribution du nombre de kmer pour chaque chromoses:'+str(total)+',pangenome:'+str(nb_kmer_chr[len(nb_kmer_chr)-1])
    plt.title(title)
    plt.savefig("nb_kmer_chr.pdf")
    plt.show()


    fig, ax1 = plt.subplots()

    # Création du barplot sur l'axe 1
    ax1.bar(name,nb_gene_chr, color='blue', alpha=0.5)
    ax1.set_xlabel('Localisation')
    ax1.set_ylabel('Nb gene', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    # Création du second plot sur l'axe 2
    ax2 = ax1.twinx()  # Partage de l'axe x
    ax2.plot(name, nb_kmer_chr, marker='o', color='red')
    ax2.set_ylabel('Nb kmer', color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    # Affichage de la figure
    plt.show()

print("Lancement de info_gene_kmer_in_genome")
info_gene_kmer_in_genome(bed)
