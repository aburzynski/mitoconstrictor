# mitoconstrictor

Set of python scripts for mtDNA annotation. Run as:
python gff2gb.py fasta_file
Prerequisites:
python 2
Biopython
(several other python libraries easily installed with pip)
wise2
tmhmm
infernal
arwen
glimmer3
critica
vmatch
blast
Databases:
custom blast formatted database of mitochondrial genomes
hmm v2 profiles of mitochondrial proteins (included, can be build with hmmer tools)
cm profiles of rRNAs and tRNAs (included, can be customized with locarna)
pfam database (non-essential)
Produces annotated GenBank file as well as several pdf and text files describing the features of the mitogenome.
