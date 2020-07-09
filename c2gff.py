import sys
import os
import os.path
from glob import glob
import math

verbose = False
#verbose = True

file_path = sys.argv[1]
seqid = file_path.split(".",1)[0]

ttable_id = sys.argv[2]

temp_file = 'tmp.fa'

fasta_file = seqid + '.fa'
with open(fasta_file, 'r') as f:
 contents = f.read()
 linie = contents.split('\n')
 sequence_len = len(''.join(linie[1:]))


with open(temp_file, 'w') as f:
 f.write('>' + seqid + '\n')
 f.write('\n'.join(linie[1:]))
 f.write('\n'.join(linie[1:]))

file_path = seqid + '3.cds'
if not (os.path.isfile(file_path)):
   command = 'critica-script tmp '+ttable_id
   if verbose: print command
   os.system(command)
   os.rename('tmp3.cds',file_path)
   for n in glob('tmp*'):
     os.remove(n)
else:
  print 'Re-using ',file_path

with open(file_path, 'r') as critica_input:
	critica_lista_linii = critica_input.readlines()

source = "critica"
ftype = "CDS"
phase = "0"
ln = 0
with open(seqid + ".gff", 'a') as gff:
 for linia in critica_lista_linii[0:-1]:
	kol = linia.split()
        start_str = kol[1]
        end_str = kol[2]
#        score = str(1 - float(kol[3]))
#        score = str(1 / float(kol[3]))
        score = str(int(-1*math.log(float(kol[3]))))
        start = int(start_str)
        end = int(end_str)
        strand = '+'
        if start > end:
          strand = '-'
          start = end
          end = int(start_str)
        if verbose: print start, end,sequence_len
        if (start < sequence_len):
          if verbose: print 'accepted' 
          ln = ln + 1
          aux = "product=orf" + str(ln)
          gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
	  gff.write(gotowa_linia)
	
#more meaningfull score is needed. If indeed kol[3] represents probability-like stat
# then more approprate than 1-x would be 1/x. But then this would lead to very high scores...
# Moreover, there is a risk of division by zero error...
# Both could be solved by some sort of normalization 
# but this would require storing everything in memory and doing second pass just for scores.
# Not worth the effort. Instead 1/x could use (more appropriate) -log(x). 
# Same problem with zero but this should never happen anyway.