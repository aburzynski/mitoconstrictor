import sys
import os
import os.path
from glob import glob

verbose = False

file_path = sys.argv[1]
seqid = file_path.split(".",1)[0]

fasta_file = seqid + '.fa'
with open(fasta_file, 'r') as f:
 contents = f.read()
 linie = contents.split('\n')
 sequence_len = len(''.join(linie[1:]))

idx_file = fasta_file+'.tis'

if not (os.path.isfile(idx_file)):
   command = 'mkvtree -pl -allout -dna -db '+ fasta_file
   if verbose: print 'running', command
   os.system(command)
else:
   if verbose: print 're-using indexes'

ln = 0

for option in ['tandem', 'supermax','p']:
 results_f = fasta_file+'.'+ option

 if not (os.path.isfile(results_f)):
   command = 'vmatch -'+option+' -l 14 -pp chain global '+ fasta_file + ' > '+ results_f
#   command = 'vmatch -e 1 -p -d -'+option+' -l 14 '+ fasta_file + ' > '+ results_f
   if verbose: print 'running', command
   os.system(command)
 else:
   if verbose: print 're-using',results_f

 with open(results_f, 'r') as v_input:
   lista_linii = v_input.readlines()

   source = "vmatch"
   ftype = "repeat_region"
   phase = "0"
   strand ='.'

   with open(seqid + ".gff", 'a') as gff:
     for linia in lista_linii:
       if linia[0] <>'#':    
	kol = linia.split()
        start1 = int(kol[2])+1
        size1 =  int(kol[0])
        score = kol[9]
        end1 = start1+size1
        start2 = int(kol[6])+1
        size2 = int(kol[4])
        end2 = start2+size2
        if (start1 < sequence_len):
          if start2 == (end1 + 1):
           start2 = start2 + 1
           end1 = end1 - 1
          ln = ln + 1
          aux = "ID=repeat" + str(ln)+';repeat_type='+option
          gotowa_linia1 = seqid + '\t' + source + '\t' + ftype + '\t' + str(start1) + '\t' + str(end1) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
          gotowa_linia2 = seqid + '\t' + source + '\t' + ftype + '\t' + str(start2) + '\t' + str(end2) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
          bstart=start1
          bend=end2
          if end2-start1 > sequence_len /2:
            bstart=start2
            bend=end1+sequence_len
            gff.write(gotowa_linia2)
            gff.write(gotowa_linia1)
          else:
            gff.write(gotowa_linia1)
            gff.write(gotowa_linia2)
          bogus_line = '\t'.join([seqid,'bogus','bogus',str(bstart),str(bend),score,strand,phase,aux])+'\n'
          gff.write(bogus_line)

#vmatch uses zero-based indices!! Need to add one to starts!
# make sure start2>end1+1, this rises error in tbl2asn.
# should that happen make more space between repeats
# - add one to start2 AND subtract one from end1, for symmetry.
