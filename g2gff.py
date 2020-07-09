import sys
import os

glimmer_path = sys.argv[1]
seqid = glimmer_path.split(".",1)[0]

ttable_id = sys.argv[2]

fasta_file = seqid + '.fa'
with open(fasta_file, 'r') as f:
 contents = f.read()
 sequence_len = len(''.join(contents.split('\n')[1:]))

file_path = seqid + '.predict'
if not (os.path.isfile(file_path)):
   command = 'glimmer3 -z'+ttable_id+' -g150 ' + fasta_file + ' prostomia.icm ' + seqid + ' 2> /dev/null'
   os.system(command)
else:
   print "Re-using ", file_path   

with open(file_path, 'r') as f:
	glimmer_lista_linii = f.readlines()

source = 'glimmer'
ftype = 'CDS'
phase = '0'
with open(seqid + '.gff', 'a') as gff:
 for linia in glimmer_lista_linii[1:-1]:
        kol = linia.split()
        gname = kol[0]
	score = kol[4]
	strand = kol[3][:1]
        start_str = kol[1]
	end_str = kol[2]
        if strand == '-':
         start_str = kol[2]
         end_str = kol[1]
	start = int(start_str)
        end = int(end_str)
        if start > end:
          end = end + sequence_len
        aux = "product=" + gname
        gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + str(start) + '\t' + str(end) + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux + '\n'
	gff.write(gotowa_linia)
