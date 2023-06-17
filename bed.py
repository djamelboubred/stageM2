import subprocess
from subprocess import PIPE
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
#
parser.add_argument("-b", "--bed", type=str, help="Fichier au format bed")
parser.add_argument("-o", "--out", type=str, help="Prefix output")
parser.add_argument("-ab", "--abond", type=str, help="Table d'abondance")
parser.add_argument("-fq", "--fastq", type=str, help="Fichier fatsq contenant les kmers")
parser.add_argument("-gff", "--gff", type=str, help="Fichier gff contenant les annotations des genes")
parser.add_argument("-FC", "--FC", type= str, help="Fichier Comparative qui à une colonne nommé gene va être comparé au fichier contenant les nom de gène trouvé suite à l'analyse DE" )
#
args = parser.parse_args() 
file_Comparative=pd.read_excel(args.FC)
file_Comparative = file_Comparative.rename(columns={'GeneID': 'gene'})
bed=pd.read_csv(args.bed, sep='\t')
bed=pd.DataFrame(bed)
gff=pd.read_csv(args.gff, sep='\t')
gff=pd.DataFrame(gff)
liste_gene={}
liste_gene_kmer={}
#prend en argument un DataFrame
def nb_kmer_in_gene(bed):
    occurence_gene={}
    for i in range(len(bed)):
        gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
        ###Visualisation de la longueuer de chaque gène
        if gene_name in occurence_gene:
                occurence_gene[gene_name] += 1
        else:
            occurence_gene[gene_name] = 1

    print("Le nombre de gène contenue dans le fichier bed est = ",len(occurence_gene))   #donne le nombre de gène 

    C5000 = 0
    C1000 = 0
    C500 = 0
    C200 = 0
    C75 = 0
    C10 =0
    C1=0
    C50=0
    for i in occurence_gene:
        if occurence_gene[i] > 5000:
            C5000+=1
        if occurence_gene[i] > 1000:
            C1000+=1
        if occurence_gene[i] > 500:
            C500+=1
        if occurence_gene[i] > 200:
            C200+=1        
        if occurence_gene[i] > 75:
            C75+=1    
        if occurence_gene[i] > 50:
            C50+=1    
        if occurence_gene[i] > 10:
            #liste_gene_with_nb_kmer_10.append(occurence_gene[i].values())
            C10+=1    
        if occurence_gene[i] > 1:
            C1+=1        

    nb_gene_filtrer=[C1,C10,C50,C75,C200,C500,C1000,C5000] #liste qui contient le nombre de gène qui on un nombre de kmer supérieur au filtre déterminer 
    filtre=["1","10",'50','75','200','500','1000','5000']
    print(nb_gene_filtrer)
    #
    ##visualisation
    #
    gene = list(occurence_gene.keys()) # occurence gene
    nb_kmer = list(occurence_gene.values()) # occurence du nombre de kmer par gene
    plt.figure()
    plt.hist(nb_kmer,bins=100, color='blue')
    plt.xlabel("Nb kmer")
    plt.ylabel("Nb Gene")
    plt.title('distribution du nombre de kmer par gène')
    plt.savefig("nb_kmer_in_gene.pdf")
    plt.show()
    #plot le nombre de gène ayant un nombre de kmer supérieur au filtre choisis
    #
    plt.figure()
    plt.bar(filtre, nb_gene_filtrer, color='blue')
    plt.xlabel("Nb kmer")
    plt.ylabel("Nb gène")
    plt.title('Distribution du nombre de gène ayant un nombre de kmer supérieur au filtre:')
    plt.savefig("nb_kmer_in_gene_filtre.pdf")
    plt.show()


    #regarde le nombre de kmer associé au gène est sort une table du nombre de gène ayant le nombre de kmer minimum choisis
    return occurence_gene

#récupere le nombre de gène pour chaque chromosomes dans le fichier gff
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

    print(nb_gene_chr)

    #plt.figure()
    #plt.subplot(1, 2, 1)
    width=-0.4
    x = list(range(len(name)-1))
    y = nb_gene_chr[:-1]
    last_x = [len(name)-1]
    last_y = nb_gene_chr[-1:]

    # Définir la couleur de l'histogramme
    colors = ['orange'] * (len(name)-1) + ['r']
    # Créer un histogramme avec deux couleurs différentes
    plt.bar(name, nb_gene_chr, width=width,color=colors, label="GFF",align='edge')
    #plt.xlabel("Chromosomes")
    #plt.ylabel("Nb Gene")
    #title='Distribution du nombre de gene pour chaque chromosomes et dans le pangenome (gff)'
    #plt.title(title)
    #plt.savefig("nb_gene_chr_in_gff.pdf")
    #plt.show()

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
            # Compter le nombre d'occurrences de chaque gène
        #counts = {}
        #for gene in liste_gene_chr:
        #    counts[gene] = counts.get(gene, 0) + 1
#
        ## Compter le nombre de gènes qui apparaissent au moins 10 fois
        #num_genes_with_min10_count = 0
        #num_genes_with_min50_count = 0
        #num_genes_with_min75_count = 0
        #for count in counts.values():
        #    if count >= 10:
        #        num_genes_with_min10_count += 1
        #    if count >= 50:
        #        num_genes_with_min50_count += 1
        #    if count >= 75:
        #        num_genes_with_min75_count += 1
        #liste_gene_chr_with_10_kmer_min.append(num_genes_with_min10_count)
        #liste_gene_chr_with_50_kmer_min.append(num_genes_with_min50_count)
        #liste_gene_chr_with_75_kmer_min.append(num_genes_with_min75_count)
        nb_gene_chr.append(len(set(liste_gene_chr)))
        
        nb_kmer_chr.append(len(liste_gene_chr))

    print(nb_gene_chr)
    print(nb_kmer_chr)
    print(liste_gene_chr_with_10_kmer_min)
    print(liste_gene_chr_with_50_kmer_min)
    print(liste_gene_chr_with_75_kmer_min)

    plt.figure()

    x=nb_gene_chr
    y=nb_kmer_chr
    corr_matrix = np.corrcoef(x, y)
    corr = corr_matrix[0, 1]

    # afficher la corrélation
    print("La corrélation entre les deux listes est :", corr)
    
    # créer un graphique pour montrer la corrélation
    plt.scatter(x, y)
    plt.title("Corrélation entre deux listes")
    plt.xlabel("Liste X")
    plt.ylabel("Liste Y")
    plt.show()


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


def liste_to_tsv(output,liste_gene):
    output=output+".tsv"
    with open(output, 'w', newline='') as f:
    # Créer un objet writer pour écrire dans le fichier TSV
        writer= csv.writer(f, delimiter='\t')

        # Écrire chaque clé et valeur du dictionnaire en tant que ligne dans le fichier TSV
        for key, value in liste_gene.items():
            writer.writerow([key, value])
    file=pd.read_csv(output, sep='\t', names=['gene', 'nb_kmer'])
    
    file=check_gene_presence(file,file_Comparative)
    file.to_csv(output, sep='\t', index=True)
def check_gene_presence(df1, df2):
    # Créer une colonne Tag dans le dataframe 1
    df1['Tag'] = df1['gene'].apply(lambda x: 'P' if x in df2['gene'].tolist() else 'A')
    return df1    
    #Appliquer la fonction à la colonne 'Gene' du DataFrame 2
    #file_Comparative['Tag'] = file_Comparative.apply(lambda row: intersect_file_Comparative_gene_DE_gene_with_kmer(row['gene'], file), axis=1, args=(file,))

    #file_Comparative['Tag'] = file_Comparative['gene'].apply(intersect_file_Comparative_gene_DE_gene_with_kmer(file))

#
## Plus utile car kmer.fastq modifié
#
   
##def convert_kmer(bed,kmer_fastq): #prend un fichier bed et un fichier fastq contenant les kmer associé au numéro (fichier fastq créer de  manière manuelle)
##    df = pd.DataFrame(columns=['gene', 'k-mer']) #création d'un dataframe qui va contenir les kmer contenue dans chacun des gènes
##    cor=False
##    with open(kmer_fastq, 'r') as f:
##        lines = f.readlines()
##        for i in range(len(bed)):
##            gene_name = bed.iloc[i, 20].strip().split(';')[1].split('=')[1]  
##            for kmer in lines:
##                if kmer[0] == ">":
##                    kmer= kmer.strip().split('>')[1]
##                    if int(bed.iloc[i,3]) == int(kmer):                
##                        #print(i)
##                        cor=True
##                        #print(kmer)
##                if cor==True:
##                    print("test")
##                    df.loc[i] = [gene_name, kmer]
##                    print(df)
##                    cor=False
##                    break
##    print("Voici la liste de gen_kmer",df)                
##    return df    


def liste_gene_pvalue_tsv(output,bed,abondance):   #prend en entrée une table d'abondance qui contient les pvalue de chaque kmer et un data frame qui contient 2 colonnes [gene, kemr]        
    #output=output+"_pvalue.tsv"
    genetokmer=pd.DataFrame(columns=['gene', 'k-mer'])
    for i in range(len(bed)):
        gene_name=bed.iloc[i, 20].strip().split(';')[1].split('=')[1]
        kmer=bed.iloc[i, 3].strip().split('_')[2]
        genetokmer.loc[i]=[gene_name,kmer]
    print(genetokmer)
    f = pd.read_csv(abondance, sep='\t')
    f.rename(columns={'tag': 'k-mer'}, inplace=True)
    resultat = pd.merge(genetokmer, f, on='k-mer', how='left')
    print(resultat)
    df=resultat[['gene','k-mer', 'pvalue']]
    print(df)
    df.to_csv('geneDe_pvalue.tsv', sep='\t', index=False)   # a supprimer mettre un appel de la fonction à la place
    return df

#def intersect_file_Comparative_gene_DE_gene_with_kmer(gene,file):
#    # Définir la fonction qui ajoute le tag si le gène est présent dans le DataFrame 1
#    if gene in file['gene'].tolist():
#        return 'P'
#    else:
#        return 'A'



print("Lancement de nb_kmer_in_gene")
liste_gene=nb_kmer_in_gene(bed)

#print("Lancement de nb_gene_in_gff")
#nb_gene_in_gff(gff)

print("Lancement de info_gene_kmer_in_genome")
info_gene_kmer_in_genome(bed)

print("Lancement de liste_gene_tsv")
liste_to_tsv(args.out,liste_gene)

print("add_tags")
#print("Lancement de convert_kmer")
#liste_gene_kmer=convert_kmer(bed,args.fastq)

print("Lancement de liste_gene_pvalue")
liste_gene_pvalue=liste_gene_pvalue_tsv(args.out,bed,args.abond)

# Appliquer la fonction à la colonne 'Gene' du DataFrame 2
#file_Comparative['Tag'] = file_Comparative['gene'].apply(intersect_file_Comparative_gene_DE_gene_with_kmer(liste_gene_pvalue))





