import sys
import os
from Bio import SeqUtils
from Bio.Seq import Seq
import RNA

verbose = True
verbose = False
experimental=True
experimental=False

arwen_path = sys.argv[1]
ttable_id = sys.argv[2]
seqid = arwen_path.split(".",1)[0]

#print arwen_path,ttable_id

fasta_file = seqid + '.fa'
with open(fasta_file, 'r') as f:
 contents = f.read()
 sequence = ''.join(contents.split('\n')[1:])
 sequence_len = len(sequence)
 sequence+=sequence

file_path = seqid + '.ar'
if not (os.path.isfile(file_path)):
   command = 'arwen -br -w -gc' + ttable_id +' '+ fasta_file + ' > ' + file_path
#   command = 'arwen -br -w -ps90 -gc' + ttable_id +' '+ fasta_file + ' > ' + file_path

#   print command
   os.system(command)
else:
   print "Re-using ", file_path   

with open(file_path, 'r') as f:
	arwen_lista_linii = f.readlines()
data_lines = arwen_lista_linii[2::3]
sequence_lines = arwen_lista_linii[3::3]
structure_lines = arwen_lista_linii[4::3]
#print len(arwen_lista_linii),arwen_lista_linii[1]

source = "arwen"
ftype = "tRNA"
phase = "."
with open(seqid + ".gff", 'a') as gff:
 for i in range(0,len(data_lines)):
        linia = data_lines[i]
	elementy_linii = linia.split()
	kolumna4i5 = elementy_linii[2].split(',')
	start_str = kolumna4i5[0]
	end_str = kolumna4i5[1][:-1]
	if start_str[0] == "c":
		start_str = start_str[2:]
		strand = "-"
	else: 
		strand = "+"
                start_str = start_str[1:]
        start = int(start_str)
        end = int(end_str)
        if start > end:
         end = end + sequence_len
        f_length = end-start+1
        f_struct = '.'*f_length
        for j in range(0,f_length):
         try:
          if structure_lines[i][j] in ['(',')']:
           f_struct = f_struct[:j]+structure_lines[i][j]+f_struct[j+1:]
          elif structure_lines[i][j] =='A':
           f_struct = f_struct[:j]+'_'+f_struct[j+1:]
         except:
          pass

#        md = RNA.md()
        RNA.cvar.temperature = 20
## temperature set a bit lower than default for now... Consider more fancy adjustment.
## free energy of predicted structure to be used as score. 
## Sequence output by arwen truncated to structure... 
## Not sure why this is needed, but the extra nt are there and need to be removed...
## the sequence may contain dots. They substitute IUPAC ambiguities.
## Consider getting the sequence from fasta instead - still, the warning will remain

        if verbose: print '\n',linia
        asequence = sequence_lines[i][:f_length].replace('.','n')
        bsequence=sequence[start-1:end]
        if strand == "-":
         bsequence = str(Seq(bsequence).reverse_complement())
        mfe = RNA.energy_of_structure(bsequence,f_struct,0)
        if verbose: 
          print f_struct, mfe
          print asequence
          print bsequence
        score = str(-1*mfe)
        
#	score = elementy_linii[3] #this was actually wrong - it is relative anticodon position!

	gname = elementy_linii[1].split("-")[1]
        label = SeqUtils.seq1(gname)
        anticodon_position = int(elementy_linii[3])-1
#        anticodon = elementy_linii[4]
        if True:###gname in ['Ser','Leu','Met']:
#          print elementy_linii[4]
          if gname in ['Ser','Arg','Gly'] and elementy_linii[4][2:4] == 'ct': label = label +'2'
          if ttable_id == '4' and gname == 'Arg' and elementy_linii[4][2:4] == 'cg': label = label + '1'
          if ttable_id == '13' and gname == 'Gly' and elementy_linii[4][2:4] == 'cc': label = label + '1'
          if gname == 'Ser' and elementy_linii[4][2:4] == 'ga': label = label +'1'
#          if gname == 'Leu' and elementy_linii[4][3] == 'g': label = label + '1'
#          if gname == 'Leu' and elementy_linii[4][3] == 'a': label = label + '2'
          if gname == 'Leu' and elementy_linii[4][3] == 'g': label = label + '2'
          if gname == 'Leu' and elementy_linii[4][3] == 'a': label = label + '1'
          if gname == 'Met' and elementy_linii[4][1] == 't': label = label + '2'
	aux = "product=tRNA-" + gname  + ";anticodon=" + elementy_linii[4] + ';anticodon_position='+ str(anticodon_position)+";label="+label+';structure='+f_struct
#        print aux
	gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
	gff.write(gotowa_linia)
##Structure saved in "forward" context, needs to be reversed for '-' features

## Variation in anticodon and aa specificity by not relying solely on arwen judgement implemented in the following
## TEST it! Use "experimental" switch for easy testing; seems to make things worse, not better,
##          so it is not used by default in production.
## checking +/- one position of anticodon... Exactly two alternatives are checked.
## to be considered must meet the following conditions:
#  -enough space around so at least one more dot before hitting stem
#  -preceding T and trailing R
#  The same score and location but different gname (label), anticodon and anticodon position (and modified structure)
## starts with anticodon, then check the conditions, then translate (revcomplemented) it into gname, 
## the rest should be easy (copied from the main loop).

        if experimental:
         if verbose: print 'considering early alternative anticodon'
         alternative_early_anticodon_position = anticodon_position - 2 #python counting vs gff counting...
         if bsequence[alternative_early_anticodon_position].upper()== 'T' and (bsequence[alternative_early_anticodon_position+5].upper()=='A' or bsequence[alternative_early_anticodon_position+4].upper()=='G'):
          if verbose: print 'sequence feasible', bsequence[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
          if f_struct[alternative_early_anticodon_position]<>'(':
           if verbose: print 'structure feasible', f_struct[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
           anticodon=bsequence[alternative_early_anticodon_position+1:alternative_early_anticodon_position+4]
           if verbose: print anticodon,'passed sructural constrains'
           aminoacid=Seq(anticodon).reverse_complement().translate(table=int(ttable_id))
           if verbose: print aminoacid, 'specificity'
           gname= SeqUtils.seq3(aminoacid)
           label=str(aminoacid)
           if gname in ['Ser','Arg','Gly'] and anticodon[1:3] == 'ct': label = label +'2'
           if ttable_id == '4' and gname == 'Arg' and anticodon[1:3] == 'cg': label = label + '1'
           if ttable_id == '13' and gname == 'Gly' and anticodon[1:3] == 'cc': label = label + '1'
           if gname == 'Ser' and anticodon[1:3] == 'ga': label = label +'1'
           if gname == 'Leu' and anticodon[2] == 'g': label = label + '2'
           if gname == 'Leu' and anticodon[2] == 'a': label = label + '1'
           if gname == 'Met' and anticodon[0] == 't': label = label + '2'
           e_structure=f_struct[:alternative_early_anticodon_position+1]+'___.'+f_struct[alternative_early_anticodon_position+5:]
	   aux = "product=tRNA-" + gname  + ";anticodon=(" + str(anticodon) + ');anticodon_position='+ str(anticodon_position-1)+";label="+label+';structure='+e_structure
           if verbose: 
             print e_structure
             print bsequence
	   gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
	   gff.write(gotowa_linia)
          else:
           if verbose: print 'bad structure',f_struct[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
         else:
           if verbose: print 'wrong context', bsequence[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
         if verbose: print 'considering late alternative anticodon'
         alternative_early_anticodon_position = anticodon_position  #python counting vs gff counting...
         if bsequence[alternative_early_anticodon_position].upper()== 'T' and (bsequence[alternative_early_anticodon_position+5].upper()=='A' or bsequence[alternative_early_anticodon_position+4].upper()=='G'):
          if verbose: print 'sequence feasible', bsequence[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
          if f_struct[alternative_early_anticodon_position+4]<>')':
           if verbose: print 'structure feasible', f_struct[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
           anticodon=bsequence[alternative_early_anticodon_position+1:alternative_early_anticodon_position+4]
           if verbose: print anticodon,'passed sructural constrains'
           aminoacid=Seq(anticodon).reverse_complement().translate(table=int(ttable_id))
           if verbose: print aminoacid, 'specificity possible'
           gname= SeqUtils.seq3(aminoacid)
           label=str(aminoacid)
           if gname in ['Ser','Arg','Gly'] and anticodon[1:3] == 'ct': label = label +'2'
           if ttable_id == '4' and gname == 'Arg' and anticodon[1:3] == 'cg': label = label + '1'
           if ttable_id == '13' and gname == 'Gly' and anticodon[1:3] == 'cc': label = label + '1'
           if gname == 'Ser' and anticodon[1:3] == 'ga': label = label +'1'
           if gname == 'Leu' and anticodon[2] == 'g': label = label + '2'
           if gname == 'Leu' and anticodon[2] == 'a': label = label + '1'
           if gname == 'Met' and anticodon[0] == 't': label = label + '2'
           e_structure=f_struct[:alternative_early_anticodon_position]+'.___'+f_struct[alternative_early_anticodon_position+4:]
	   aux = "product=tRNA-" + gname  + ";anticodon=(" + str(anticodon) + ');anticodon_position='+ str(anticodon_position-1)+";label="+label+';structure='+e_structure
           if verbose: 
            print e_structure
            print bsequence
	   gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
	   gff.write(gotowa_linia)
          else:
           if verbose: print 'bad structure',f_struct[alternative_early_anticodon_position:alternative_early_anticodon_position+5]
         else:
           if verbose: print 'wrong context', bsequence[alternative_early_anticodon_position:alternative_early_anticodon_position+5]

