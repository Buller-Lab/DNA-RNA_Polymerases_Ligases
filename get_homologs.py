import os
import pandas as pd
import subprocess
import sys
import argparse

parser = argparse.ArgumentParser(description='Run MMseqs2 easy-search + linclust pipeline')
parser.add_argument('-i', '--input', required=True, help='Input seeds FASTA file')
parser.add_argument('-d', '--database', required=True, help='Path to MMseqs2 database')
parser.add_argument('-o', '--output-dir', required=True, help='Output directory / project name')
args = parser.parse_args()

input_fasta = args.input
database = args.database
output_dir = args.output_dir

if not os.path.exists(input_fasta):
    print(f"Error: Input file not found: {input_fasta}")
    sys.exit(1)

path = output_dir
if not os.path.exists(path):
    os.makedirs(path)
    print("Created project directory")
else:
    print("Output directory already exists — will overwrite files")

subprocess.run(
    f"mmseqs easy-search {input_fasta} {database} {output_dir}/alnRes.m8 tmp "
    f"--format-output 'target,tseq,pident,tlen,taxid,taxname,taxlineage' "
    f"--min-seq-id 0.25 --max-seqs 100000000000 --split-memory-limit 90G",
    shell=True, check=True
)

hits = pd.read_csv(
    f'{output_dir}/alnRes.m8',
    sep='\t',
    header=None,
    comment='#'
)
hits.columns = ['ID', 'Sequence', 'SeqID_Seed', 'Length', 'TaxID', 'Organism', 'TaxInfo']

with open(f'{output_dir}/hits.fasta', 'w') as hits_fasta:
    for _, row in hits.iterrows():
        hits_fasta.write(f">{row['ID']}\n{row['Sequence']}\n")

subprocess.run(
    f"mmseqs easy-linclust {output_dir}/hits.fasta {database} {output_dir}/clusterRes_60 tmp "
    f"--min-seq-id 0.6 -c 0.8 --cov-mode 0",
    shell=True, check=True
)

cluster_60 = pd.read_csv(
    f'{output_dir}/clusterRes_60_cluster.tsv',
    sep='\t',
    header=None,
    names=['representative', 'member'],
    comment='#'
)

hits.to_csv(
    f'{output_dir}/alnRes.m8',
    sep='\t',
    index=False,
    header=['target', 'tseq', 'pident', 'tlen', 'taxid', 'taxname', 'taxlineage']
)

cluster_60.to_csv(
    f'{output_dir}/clusterRes_60_cluster.tsv',
    sep='\t',
    index=False,
    header=['representative', 'member']
)

print("Done. Final files with headers:")
print(f"  - {output_dir}/alnRes.m8          (with header)")
print(f"  - {output_dir}/clusterRes_60_cluster.tsv  (with header)")