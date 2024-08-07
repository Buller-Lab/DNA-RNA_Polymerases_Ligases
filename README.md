# CyclAmidImid

Code and databases for the publication: 
## Nature’s Toolbox for the Hydrolysis of Lactams and Cyclic Imides
Peter Stockinger, [a] and Rebecca Buller*[a]

[a] Competence Center for Biocatalysis, Institute of Chemistry and Biotechnology, Zurich University of Applied Sciences, Einsiedlerstrasse 31, 8820 Wädenswil, Switzerland

Keywords: Amidase · Imidase · Hydrolase · Amidohydrolase ·Lactamase ·Lactams ·Imides ·Cyclic ·Amides ·Biocatalysis · Hydantoinases · Dihydropyrimidinases 

![image](https://github.com/user-attachments/assets/33f5a252-cd52-4b03-b85d-a35f022a7b8f)


## Background

The historical development of various naming and classification schemes to describe enzymes with hydrolytic activity toward lactams and cyclic imides has contributed to a fragmented landscape that poses challenges in the search for novel biocatalysts, particularly for organic synthesis in the pharmaceutical and fine chemical industry. To identify (novel) hydrolytic enzymes suitable for application, we thus suggest to not only rely on functional classifiers (such as EC numbers: https://brenda-enzymes.org/) but also use structural folds as a search criterium (CATH IDs: https://www.cathdb.info/). In this context, we have generated a data repository of promising hydrolase superfamilies to facilitate the identification of cyclic imide or lactam hydrolyzing enzymes.

## Databases & Webapp

For each structural fold described in this review, thousands of enzyme sequences were assembled in a searchable database, linking their NCBI accessions, organisms of origin, and SeqID to (experimentally) confirmed lactamases/imidases. To provide an additional filtering opportunity, we have added extremophilic annotations (thermophilic, psychrophilic, and halophilic), which can support mining for enzymes with improved process properties, such as thermal- or solvent stability.

Databases can be found as .m8 files (can be opened via Excel), Cytoscape Sessions, as well as in an interactive webapp:
https://ccbiozhaw.github.io/CyclAmidImid/

## Methods and Code Availability
Databases of each fold archetype (superfamily)  were build using MMseqs2 (https://github.com/soedinglab/MMseqs2) easy-search (against NCBI NR database) and easy-cluster (60% SeqID cutoff) functionalities.
Seed sequences are provided as fasta files with the ending _seeds.fasta. Extremophilic annotations (thermophilic, psychrophilic & halophilic) were collected from BacDive database (https://bacdive.dsmz.de/) and were linked via TaxIDs.
Python scripts that were used to collect and cluster the sequences are also provided.

## Citation
Please cite the following review if you found this ressource helpful:
DOI t.b.d.
