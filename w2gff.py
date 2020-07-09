import sys
import re
import os
import os.path

#first run genewise on provided fasta file against all the usual hmms,
#in sequence, immediately parsing the output into gff file
#for empty outputs repeat analysis with lower cutoff

verbose=True
verbose=False
fasta_file = sys.argv[1]
fasta_name = fasta_file.split('.', 1)
mitochondria =['ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND5','ND6','ND4L','ATP8']

ttable_id = sys.argv[2]

#circular genomes assumed => work on concatenated file.
#read original fasta and write it twice to temporary file, second time without the first line
#consider re-writing this dirty hack using Biopython

temp_file = fasta_name[0] + '.tmp'

with open(fasta_file, 'r') as f:
 contents = f.read()
 linie = contents.split('\n')
 seq_length = len(''.join(linie[1:]))

with open(temp_file, 'w') as f:
 f.write(contents)
 f.write('\n'.join(linie[1:]))

for gene in mitochondria:
 cutoff = 50.0
 redo = False
 print 'seeking ', gene
 file_path = fasta_name[0] + '.' + gene
 while True: #iteration of cutoff if the output is empty
#  if True:
  if not (os.path.isfile(file_path)) or redo:
   command = 'genewisedb '+ gene +'.hmm ' + temp_file + ' -hmmer -codon codon.table.'+ttable_id+' -alg 333 -aalg 333L -cut ' + str(cutoff) +' -silent -kbyte 100000000 -nohis -sum -gff > ' + file_path
#   print command
   os.system(command)
   print str(cutoff)
   with open(file_path, 'a') as f:
    f.write(str(cutoff))
  else:
    with open(file_path, 'r') as f:
      linie = f.readlines()
    cutoff = float(linie[-1:][0])
    print 'Re-using existing file, cutoff: ',cutoff,
  with open(file_path, 'r') as file_input:
	wise_data = file_input.read()
	zakres_linii = wise_data.split("//")
  good = False
  while len(zakres_linii) > 1: #check for empty output
                            #have to consider low scoring output: if finding score was higher than filtering score due to algorithm pecularities..
                            #need another control variable analogous to redo. Can't re-use redo...
   linie = zakres_linii[1].split('\n')

   names = file_path.split(".",1)
   fname = names[0]
   fext = names[1]
#   score = '.'
   with open(fname+".gff", 'a') as gff:
    for linia in linie[1:-1]:
        kol = linia.split('\t')
        ftype = kol[2]        
        zapis = ftype == 'cds'
#        seqid = kol[0]
        seqid = fname
        source = kol[1]
	start = int(kol[3])
	end = int(kol[4])
        if not zapis:
          if verbose: print 'supplementary'
          score = kol[5]
        phase = kol[7]
	strand = kol[6]
        if strand == "-":
		start = int(kol[4])
		end = int(kol[3])
        if start > seq_length:    #location outside, no need to add
          zapis = 0
          if verbose: print 'outside'
        if float(score) < (cutoff-1): #additional filtering, not done by genewisedb, a bit relaxed, just in case
          zapis = 0
          if verbose: print 'low score'
	attributes = "product=" + fext 
        gotowa_linia = seqid + '\t' + source + '\t' + ftype.upper() + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes +'\n'
        if zapis:
         if verbose: print 'accepted',seqid,start,end,fext,score,cutoff
         gff.write(gotowa_linia)
         good = True
        else:
         if verbose: print 'rejected', seqid,start,end,fext,score,cutoff
   print '.. OK'
   if good: 
    zakres_linii = []
   else:
    cutoff = cutoff-1
    if cutoff < 10:
      zakres_linii = []
  if not good:
   cutoff = cutoff/2
   redo = True
   if cutoff < 10:
    print 'not found.'
    break
  else:
    break
os.remove(temp_file)
