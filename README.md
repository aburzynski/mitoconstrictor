# mitoconstrictor

Set of python scripts for mtDNA annotation. Run as:</br>
python gff2gb.py fasta_file</br>
Prerequisites:</br>
-python (currently version 2.7 is supported)</br>
-Biopython
(and several other python libraries easily installed with pip)</br>
-wise2 https://www.ebi.ac.uk/~birney/wise2/</br>
-tmhmm https://services.healthtech.dtu.dk/software.php</br>
-infernal http://eddylab.org/infernal/</br>
-arwen </br>
-glimmer3 (icm profile for glimmer included)</br>
-critica </br>
-vmatch </br>
-blast+ </br>
-EMBOSS (getorf)</br>
-Databases:</br>
+database of mitochondrial genomes, blast-formatted, with corresponding lineage specification (for mgcode and profile customization, included, can be customized at will)</br>
+individual blast-formatted sub-databases of major clades (for critica, examples included)</br>
+hmmer v2 profiles of mitochondrial proteins (included, these can be customized with hmmer)</br>
+cm profiles of rRNAs and tRNAs (included, these can be customized with locarna)</br>
+pfam database (part A, hmm file pressed with hmmer, non-essential, not included) </br>
Produces annotated GenBank file as well as several pdf and text files describing the features of the mitogenome.</br>
</br>
Several dependencies run on linux only but Win10 linux subsystem is also adequate. Everything is commandline-only.</br>

Publication with some examples of usage is freely available at:</br>
https://www.nature.com/articles/s41598-020-67976-6
