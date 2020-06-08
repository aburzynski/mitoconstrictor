# mitoconstrictor

Set of python scripts for mtDNA annotation. Run as:</br>
python gff2gb.py fasta_file</br>
Prerequisites:</br>
python (currently version 2.7 is supported)</br>
Biopython </br>
(several other python libraries easily installed with pip)
wise2 </br>
tmhmm </br>
infernal </br>
arwen </br>
glimmer3 </br>
critica </br>
vmatch </br>
blast </br/>
Databases:</br>
custom blast formatted database of mitochondrial genomes </br>
hmm v2 profiles of mitochondrial proteins (included, can be build with hmmer tools) </br>
cm profiles of rRNAs and tRNAs (included, can be customized with locarna) </br>
pfam database (non-essential) </br>
Produces annotated GenBank file as well as several pdf and text files describing the features of the mitogenome.
