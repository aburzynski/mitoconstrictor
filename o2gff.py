import sys
import os
from Bio.Data import CodonTable
from math import log
#first run getorf on provided fasta file,
#then parse the output into gff file

# lets try more sophisticated length discrimination here.
# bootstrap the sequence, get orfs of the bootstrapped set, then evaluate  
# probability of obtaining such long orf by chance.
# count all orfs and longer orfs in the bootstrapped set, calculate the proportion
# Reversed proportion of longer orfs should have values comparable to bitscores!

# orf is a string containing nucleotide sequence (with stop) of the orf to be evaluated.
def prob(orf):
 orl =len(orf)
 n = max(2,2000 / orl)
 fh = open('tmp.fas','w')
 fh.write('>\n')
 fh.write(n*orf) # tune this for speed and numeric stability
 fh.write('\n')
 fh.close()
 command = 'fseqboot -reps 3 -noprogress -outfile tmp.phy tmp.fas > /dev/null 2>&1'
 os.system(command)
 command = 'getorf -table '+ ttable_id+' -find 3 -minsize 2 -circular T -noreverse -outseq tmp.bootorf tmp.phy -snucleotide > /dev/null 2>&1'
 os.system(command)
 fh = open('tmp.bootorf', 'r')
 n_orf=0
 n_long=0
 for linia in fh:
  if linia[0]=='>':
    n_orf+=1
    kolumna=linia.split()
    start = int(kolumna[1][1:])
    end = int(kolumna[3][:-1])
    if end-start > orl:
     n_long+=1
 fh.close()
 if n_long >0:
  score=str(n_orf/float(n_long))
 else:
  score=str(2*n_orf)
 return score
 # fraction of long(er) orfs is a good proxy of probability for obtaining such a long orf by chance.
 # its reverse should be comparable to bitscore. 
 # Or should it be other function than reverse? Log perhaps?
 # yeah, it should be proportional to -log2(p), not to 1/p
 # so the 'corrected' score could be log2(score)
 # unfortunately, this is way too flat. Dirty hack: use base=1.03, yields comparable numbers

# Analyticlal approach to the same. Much faster (10000x).
def aprob(orf):
# precision = 5
# precision = 10**(-1*(precision-1))
 span=len(orf)
 cspan = span / 3
 ttable=CodonTable.unambiguous_dna_by_id[int(ttable_id)]
# startlist=ttable.start_codons
 stoplist=ttable.stop_codons
 denom=float(span)
 f={}
 f['A']=orf.count('A')/denom
 f['C']=orf.count('C')/denom
 f['T']=orf.count('T')/denom
 f['G']=orf.count('G')/denom
# pstart=0
# for codon in startlist:
#  pstart=pstart+f[codon[0]]*f[codon[1]]*f[codon[2]]
 pstop=0
 for codon in stoplist:
  pstop=pstop+f[codon[0]]*f[codon[1]]*f[codon[2]]
# p=0
# for i in range(4,cspan-1):
#  contribution = pstart*(1-pstop)**i
#  p+=contribution
# total = p
# while contribution > precision:
#  i+=1
#  contribution = pstart*(1-pstop)**i
#  total+=contribution
# if p<total:
#  p=1-p/total
# else:
#  p=1/denom #any better idea?
 p=(1-pstop)**cspan
# return 1/p
 return -1*log(p,1.03)

fasta_file = sys.argv[1]
ttable_id = sys.argv[2]
fasta_name = fasta_file.split('.', 1)
orf_file_path = fasta_name[0] + '.orf'
gff_file_path = fasta_name[0] + '.gff'

# real data here, with length evaluation.

with open(fasta_file, 'r') as f:
 contents = f.read()
 linie = contents.split('\n')
 sequence_len = len(''.join(linie[1:]))

extended_seq = ''.join(linie[1:])
extended_seq = extended_seq + extended_seq 

if not (os.path.isfile(orf_file_path)):
 command = 'getorf -table '+ttable_id+' -find 3 -minsize 70 -circular T -outseq ' + orf_file_path + ' ' + fasta_file + ' -snucleotide1 2> /dev/null'
 os.system(command)
else:
 print 'Re-using ', orf_file_path
phase = '0'
seqid = fasta_name[0]

f = open(orf_file_path, 'r')
g = open(gff_file_path, 'a')
for linia in f:
 if linia[0]=='>':
  kolumna = linia.split()
  attributes = 'product=' + kolumna[0][1:]
  source = 'getorf'
  ftype = 'CDS'
  if len(kolumna) > 4 and kolumna[4] == '(REVERSE':
    strand = '-'
    start = str(int(kolumna[3][:-1])-3)
    end = kolumna[1][1:]
  else:
    strand = '+'
    start = kolumna[1][1:]
    end = str(int(kolumna[3][:-1])+3)
  score='0'
  try:
   score = str(int(aprob(extended_seq[int(start)-1:int(end)+3].upper())))
  except:
   pass
  g.write(seqid + '\t' + source + '\t' + ftype + '\t' + start + '\t' + end + '\t' + score + '\t' + strand + '\t' + phase + '\t' + attributes +'\n')
