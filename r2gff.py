import sys
import re
import os.path
#import time
import os
from Bio.SeqUtils import nt_search, seq3
from Bio.Seq import reverse_complement

def flip_struct(x):
        br_dict={'(':')',')':'(','>':'<','<':'>','[':']',']':'[','{':'}','}':'{','.':'.','_':'_',':':':',',':',','-':'-','~':'~'}
        span = len(x)
        r = x[::-1]
        for j in range(0,span):
           r = r[:j]+br_dict[r[j]]+r[j+1:]
        return r ## reversed structure with adjusted brackets.

## structure sanitization procedure:
## pass zero: replace underscores with dots - not really needed. 
## Could be left for the main program...If really needed- hairpin conservation may actually be important...
##  -do not touch *[ nn] characters in the first pass, only spaces
##  -replace each *[ nn]* with '.'*nn string (by slice-in) in structure only (second pass)
##  -global bracket adjustment still needed
##  After that, refolding should not happen

def sanitize_structure(x,y): #x structure, y sequence
   opening =  ['(','[','{','<']
   closing =  [')',']','}','>']
   legit = ['AT','TA','GC','CG','GT','TG']
   last = len(x)
   i = last
   x=x.replace('_','.') # Figures look better with this replace.
   if verbose: print "\nsanitizing"
   if verbose: print "first pass\n",x,'\n',y
   while i > 0:
        i = i - 1
        if y[i] =='-':
          if x[i] in opening:
           bcount = 1
           for j in range(i+1,len(x)):
            if x[j] in closing:
             bcount = bcount - 1
            elif x[j] in opening:
             bcount = bcount + 1
            if bcount == 0: #found a closing bracket. Slice out the annotation string (immutable!)
             old_annotation = x
             new_annotation = old_annotation[:j]+'.'+old_annotation[j+1:]
             x = new_annotation
             break
          elif x[i] in closing:
           bcount = -1
           for j in range(i-1,0,-1):
            if x[j] in closing:
             bcount = bcount - 1
            elif x[j] in opening:
             bcount = bcount + 1
            if bcount == 0: #found an opening bracket. Slice out the annotation string (immutable!)
             old_annotation = x
             new_annotation = old_annotation[:j]+'.'+old_annotation[j+1:]
             x = new_annotation
             break
          y = y[:i]+y[i+1:]
          x = x[:i]+x[i+1:]

   # second pass 
   #  -identify large inserts in sequence (y)
   #  -replace each respective fragment of structure 
   #   (equivalent to *[ nn]* in sequence) with '.'*nn string (by slice-in)
   #   use 'N' to keep the sequence in-sync
   if verbose: print '\n',x,'\n',y,'\nsecond pass: expand large inserts'
   xy=[]
   i=0
   maxi = len(y)
   while i < maxi:
    if deepverbose: print '\n',y[i:i+2]
    if y[i:i+2]=='*[':
      for j in range(i+1,len(y)):
        if y[j] =='*':
          xy.append((i,j+1))
          i = j
          break
    i+=1
#   print y, len(xy), 'inserts'
   xy = sorted(xy, reverse=True)
   for i,j in xy:
     nn = int(y[i+2:j-2])
#     print y[i:j], nn
#     print x
     x = x[:i]+'.'*nn+x[j:]
     y = y[:i]+'N'*nn+y[j:] ##Will need to replace the sequence later
#     print x
#     print y 
   if verbose: print '\n',x,'\n',y,'\n after second pass'

# need one more thing: the ends may be truncated, as indicated by <[ n]* at 5' and *[ n]> at 3'
# if this is detected, approprate number of characters has to be truncated off the ends of structure.
# how many? If n is 0, then just the length of the string, from * to the end.
# at this point there should be no other * in the structure so may use that shortcut.
# if n is not 0, not sure what happens but should still be OK.

# this is indicative of a possibility to do more precise annotation
# of truncated forms. This would require access to the model part of the alignment
# to see how much of it is missing from each side.
# doable but not here - in the body of thos script, closer to the beginning.

   if verbose: print 'end polishing'

   if y[0] == '<':
    i=y.find('*')
    x = x[i+1:]
    y = y[i+1:]
   if y[-1] == '>':
    i=y.find('*')
    x =x[:i]
    y =y[:i]

   if verbose: print '\n',x,'\n',y,'\n','after end polishing'

   # last pass, verify pairs, if not in legit substitute both with dots. 
   # It is still needed
   if deepverbose: print 'final pass'
   last = len(x)-1
   i=-1
   while i < last:
      i = i + 1
      bcount = 0
      if x[i] in closing:
        if deepverbose: print x[i],y[i],'closing at', i
        bcount -=1
        for j in range(i-1,0,-1):
         if x[j] in opening:
          bcount+=1
         elif x[j] in closing:
          bcount-=1
         if bcount == 0:
          if deepverbose: print 'match at',j,x[j],y[j]
          if not (y[i]+y[j] in legit):
           if deepverbose: print 'not goot, dotting'
           x = x[:i]+'.'+x[i+1:]
           x = x[:j]+'.'+x[j+1:]
           if deepverbose:
            print x
            print y
          elif deepverbose: print 'OK'
          break
        if bcount == 1:
          if deepverbose: print 'no opeing bracket found!'
          if deepverbose:  print 'not good, dotting'
          x = x[:i]+'.'+x[i+1:]
          x = x[:j]+'.'+x[j+1:]
          if deepverbose:
            print x
            print y
      elif x[i] in opening:
        if deepverbose:print x[i],y[i],'opening at', i
        bcount +=1
        for j in range(i+1,len(x)):
#         print bcount,j
         if x[j] in opening:
          bcount+=1
         elif x[j] in closing:
          bcount-=1
         if bcount == 0:
          if deepverbose:print 'match at',j,x[j],y[j]
          if not (y[i]+y[j] in legit):
           if deepverbose:  print 'not good, dotting'
           x = x[:i]+'.'+x[i+1:]
           x = x[:j]+'.'+x[j+1:]
           if deepverbose:
            print x
            print y
          elif deepverbose: print 'OK'
          break
        if bcount == 1:
          if deepverbose: print 'no closing bracket found!'
          if deepverbose:  print 'not good, dotting'
          x = x[:i]+'.'+x[i+1:]
          x = x[:j]+'.'+x[j+1:]
          if deepverbose:
            print x
            print y
   return x,y

##MAIN
#t0=time.time()
deepverbose = False
verbose = False
#verbose = True

file_path = sys.argv[1]
seqid = file_path.split(".",1)[0]

temp_file = 'tmp.fa'

fasta_file = seqid + '.fa'
with open(fasta_file, 'r') as f:
 contents = f.read()
 linie = contents.split('\n')
 sequence_len = len(''.join(linie[1:]))

extended_seq = ''.join(linie[1:])
extended_seq = extended_seq + extended_seq 

with open(temp_file, 'w') as f:
 f.write('>' + seqid + '\n')
 f.write('\n'.join(linie[1:]))
 f.write('\n'.join(linie[1:]))

# option --anytrunc is needed rarely, only when no hit is returned for one or two genes...
# consider conditioning its use on the results of the initial run without it.
# it would be much easier to evaluate this after gff parsing but then it will be too late... (gff already there)


file_path = seqid + '.rRNA'
if not (os.path.isfile(file_path)):
   print 'running cmscan with rRNA profiles...'
   command = 'cmscan --cpu 2 --notextw -o ' + seqid + '.rRNA rrna.cm ' + temp_file 
   os.system(command)
else:
   print 'Re-using ', file_path

###parsing structure info, not tbl file.

with open(file_path, 'r') as file_input:
	infernal_data = file_input.read()
	data = infernal_data.split("Hit alignments:",1)
        genes = data[1].split("\n>> ")[1:]

# still, genes may be a good proxy of what is there.
# if there are less than three genes, rerun with the option.
# Since the full duplicate is used we actually expect four genes...
# ... not good either: does not trigger the switch if there is one good and one poor hit of 12S...
# Have to inspect the first line to see if they are different...
models=set()
for gene in genes:
 model=gene.split('\n')[0]
 models.add(model)
models=list(models)
# TO DO: prepare lineage specific cm files for problematic cases (instead)

if len(models) < 2:
   print 'too few complete rRNA genes, scanning for truncated genes...'
   command = 'cmscan --anytrunc --cpu 3 --notextw -o ' + seqid + '.rRNA rrna.cm ' + temp_file 
   os.system(command)

with open(file_path, 'r') as file_input:
	infernal_data = file_input.read()
	data = infernal_data.split("Hit alignments:",1)
        genes = data[1].split("\n>> ")[1:]

ftype = 'rRNA'
source = "infernal"
phase = "."
relaxed=['!']

with open(seqid + ".gff", 'a') as gff:

 for gene in genes:
  linia = gene.split('\n')
  column = linia[3].split()
  if column[1] in relaxed: 
## OK, so linia[0] has only name, linia[1] headers, linia[2] is a separator, linia[3] has all numbers
## linia[4]is empty, linia[5] has marks of unclear purpose, linia[6] holds the structure
## linia[7] is model aln, linia[8] shows consensus, linia[9] is target aln, 
## linia[10] marks confidence(?), maybe one more empty line is present.
## data in linia[3] are parsed first. Interesting addition is column[14]
## it may indicate possible truncation, if either 5' or 3' are in it.
## should that happen, closer inspection of linia[7], just its ends
## would allow better annotation of the truncated feature.

## this should be done before sanitization but after isolating target structure.

## cmscan stores reversed start/end coordinates for '-' strand
## better to correct this immediately, make sure it is consistent
   start = column[9]
   end = column[10]
   score = column[3]
   product = linia[0].strip()
## now parse structure and sequence (aln)
   structure = linia[6].rstrip(' CS')
   aln = linia[9].upper().replace('U','T')
   short_structure = structure.lstrip()
   offset = len(structure) - len(short_structure)
   structure = short_structure
   aln = aln[offset:len(structure)+offset]
   aux=''
# so about here potential truncation should be checked.
   if '5' in column[14] or '3' in column[14]: #this triggers model structure check
    model = linia[7][offset:len(structure)+offset]
    missing5 = 0
    missing3 = 0
#   possible 5' truncation, extract the number... This should be available also in column[6]-1...
    if model[:2]=='<[':
      for j in range(1,len(model)):
        if model[j] =='*':
           missing5 = int(model[2:j-1])
           break
#  possible 3' truncation, get this number too... This one is nowhere else to find
    if model[-2:]==']>':
      for j in range(len(model)-2,0,-1):
        if model[j:j+2]=='*[':
           missing3 = int(model[j+2:-2])
           break
# now check which one is bigger enough and leave the annotation in aux
    if missing5 > 5*missing3:
      aux="note=_3';"
    if missing3 > 5*missing5:
      aux="note=_5';"
# less acute cases are left unannotated

## need conversion to common coordinates/lengths
   if deepverbose: print '\nstructure:',structure,'\n      aln:',aln
   (better_structure, better_aln) = sanitize_structure(structure,aln) 
   structure = better_structure
   aln = better_aln
   strand = column[11]
   if strand == '-':
    start = column[10]
    end = column[9]
    structure=flip_struct(structure)
   aln_len=int(end)-int(start)+1
   if len(aln) == aln_len == len(structure):
     if deepverbose: print len(aln),len(structure)
   else: # should never happen, parsing failed
     print "can't adopt cmscan structure"
#     print structure
#     print aln
     structure='.'*aln_len
   if verbose:
    print '' 
    print product
    print structure,start,end
    print len(structure), aln_len
   aux += 'product=' + product+';label='+product+';structure='+structure
#when to do flipping of structure on reverce strand? 
# For tRNA it is in the main program but this is convenient since anticodon is decided there.
# but for rRNA it will be better to store final structures in gff. Potential later sanitization may require re-flipping!
   gotowa_linia = seqid + '\t' + source + '\t' + ftype + '\t' + start + '\t' + end + '\t' + score + '\t' + strand + '\t' + phase + '\t' + aux +'\n'
   if int(start) <= sequence_len:
         gff.write(gotowa_linia)
#t1=time.time()
#print t1-t0
