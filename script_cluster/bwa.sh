    #!/bin/bash

    ##Configuration SLURM ########
    ##on définit le nom du job 
    #SBATCH --job-name=BWA_ALN_kmer

    ##définir la partition utilisée
    #SBATCH -p highmem

    ##définir le nombre de  à utiliser 
    #SBATCH -c 12

    cd /scratch
    mkdir -p ${USER}_${SLURM_JOB_ID}_bwa
    cd ${USER}_${SLURM_JOB_ID}_bwa
    mkdir -p output
    cd output
    ref=$1
    fasta_kmer=$2

    # nas3:/data3/projects/integration/OglaRS2.ADWL02-allCtgsIRIGIN_TOG5681.dedup8095-NR.fasta 
    # nas3:/data3/projects/integration/rnaseq_riz/output/clean/testing/mapping_DE/kmer_DE.fasta*
    # nas3:/data3/projects/integration/gene_60k* 
    echo "COPYING DATA" 

    scp ${ref} .
    scp ${fasta_kmer} .
    
##Recupére le nom des fichiers
    fasta_kmer = $(basename ${fasta_kmer})
    ref = $(basename ${ref})
    path_output = $(realpath ${fasta_kmer})
    path_output = $(dirname ${path_output})
    # GERER LES FICHIER AVEC EXTENSION GZIP
    if (file ${fasta_kmer} | grep -q gzip ) ; then
        echo "UNCOMPRESSED ${fasta_kmer}"
        gzip -d ${fasta_kmer}
        fasta_kmer = $(basename ${fasta_kmer} .gz)
    fi
    if (file ${ref} | grep -q gzip ) ; then
        echo "UNCOMPRESSED ${ref}"
        gzip -d ${ref}
        ref = $(basename ${ref} .gz) 
    fi

    ##Charger les logiciel à utiliser 
    module load bioinfo/samtools/1.9
    module load bioinfo/bwamem2/2.2.1 

    #INDEXATION DE A REF
    bwa index -p index ${ref}

    echo "LANCEMENT BWA ALN"
    #LANCEMENT DE BWA ALN
    extension=$(echo ${fasta_kmer#*.})
    name_kmer = $(basename ${fasta_kmer} ${extension})
    extension=$(echo ${ref#*.})
    name_ref = $(basename ${ref} ${extension})
    output = "${name_ref}_VS_${name_kmer}"

    # LANCEMENT DE BWA ALN
    # -l taille du reads (31)
    # -n pourentage de taux d'erreur accepté lors du mapping
    bwa aln -t 12 -l 31 -n 0.03 index ${fasta_kmer} > ${output}.sai

    #CONVERSION SAI TO SAM
    bwa samse index ${output}.sai ${fasta_kmer} -f ${output}.sam

    #Etape : TRIE BAM
    # -l niveau zero de compréssion gzip , @ niveau de threads
    samtools sort -l O -@ 2 -o ${output}_sorted.bam ${output}.bam

    # Etape: Recupération des kmers mappés (-F élimine le kmer avec tags choisis (4 unmapped))
    samtools view -h -F 4 -b ${output}_sorted.bam > ${output}_sorted_mapped.bam

    #Etape: Filtre avec mapq > 1 (map= 0, correspond aux kmers non mappés )
    samtools view -h -q 1 -F 4 -b ${output}_sorted.bam > ${output}_sorted_mapped_mapQ.bam

    #ETAPE: d'indexation du bam qu'on va garder
    samtools index kmer_sorted_mapped.bam

    #ETAPE: FLAGSTATS récupere le nombre de kmers mappés pour chacune des étapes précedentes
    samtools flagstat ${output}_sorted.bam > ${output}_flagst.txt
    samtools flagstat ${output}_sorted_mapped.bam > ${output}_mapped_flagst.txt
    samtools flagstat ${output}_sorted_mapped_mapQ.bam > ${output}_mapped_mapQ_flagst.txt

    tree .
    ls -l 
    rm ${output}
    rm index*
    rm ${output}.sai
    rm ${output}*.sam
    rm ${output}_sorted.bam
    rm ${output}_sorted_mapped.bam

    #TRANSFERT DES RESULTATS
    scp -r output ${path_output}

