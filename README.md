## This is the GitHub Repository for the publication: 

#  Recent Advances in the Evolution of Polymerases and Ligases for the Synthesis of Modified Nucleic Acids

Sumire Honda Malca,[a] Peter Stockinger,[b] Miquel Estévez-Gay,[b] and Rebecca Buller*[a],[b]


[a]	Competence Center for Biocatalysis, Institute of Chemistry and Biotechnology, Zurich University of Applied Sciences, Einsiedlerstrasse 31, 8820 Wädenswil, Switzerland


[b]	Department of Chemistry, Biochemistry and Pharmacy, University of Bern, Freiestrasse 3, 3012 Bern, Switzerland



Keywords: X · Y · Z 



In this repository we provide a link to a webapp with all related enzyme families, as well as methods, data and code that were used to build the respective databases.

## Databases & Webapp

For each enzyme family described in this review, thousands of enzyme sequences were assembled in a searchable database, linking their NCBI accessions, organisms of origin, and SeqID to (experimentally) confirmed polymerases/ligses. To provide an additional filtering opportunity, we have added extremophilic annotations (thermophilic, psychrophilic, halophilic, acidophilic, and alkaliphilic), which can support mining for enzymes with improved process properties, such as thermal- or solvent stability.

Databases can be found as .m8 files (can be opened via Excel), Cytoscape Sessions, as well as in an interactive webapp:

https://buller-lab.github.io/DNA-RNA_Polymerases_Ligases/

## Methods and Code Availability
Databases of each fold archetype (superfamily)  were build using MMseqs2 (https://github.com/soedinglab/MMseqs2) easy-search (against NCBI NR database) and easy-cluster (60% SeqID cutoff) functionalities.
Seed sequences are provided as fasta files (seeds.fasta). Extremophilic annotations were collected from BacDive database (https://bacdive.dsmz.de/) and were linked via TaxIDs (see ).
The python script that was used to collect and cluster the sequences is also provided (see get_homologs.py). 

## Citation
Please cite the following review if you found this ressource helpful:
XYZ
