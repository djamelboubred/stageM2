import subprocess
from subprocess import PIPE
import argparse


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", type=str, help="Input file with in first colums the tags (kmer), without header. Obtain with DEKUPL Pipeine (raw_pvalues_counts.tsv)")
parser.add_argument("-o", "--output", type=str, help="Output name file, name of fasta file")
args = parser.parse_args()

file_in= args.input
file_out = args.output

command = ["wc","-l",f"{file_in}"]
shape = subprocess.run(command, stdout=PIPE, stderr=PIPE)
shape = int(str(shape.stdout).strip().split(' ')[0].split('\'')[1])

with open(file_in, 'r') as b, open(file_out, 'w') as f:

    for i in range(0, int(shape)):

        contig = b.readline().rstrip()

        contig = contig.split('\t')[0]
        if contig != 'contig' and contig !='tag':
            f.write(f">{i}_tags_{contig}\n{contig}\n")