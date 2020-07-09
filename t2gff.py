### TO DO: implement all codon tables (default is 5, major alternatives are: 2,4,9,13, minor alternatives are 14,24,33)
##  Need pref_code update for each ttable - check if this would be enough.

import sys
import re
import os
from Bio.SeqUtils import nt_search, seq3
import RNA
from Bio.Seq import reverse_complement

RNA.cvar.temperature = 20
#do_filter = False
do_filter = True
deepverbose = True
deepverbose = False
verbose = False
#verbose = True


## structure sanitization procedure:
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
#   if verbose: print "\nsanitizing"
   if deepverbose: print "first pass\n",x,'\n',y
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
   if deepverbose: print '\n',x,'\n',y,'\nsecond pass'
   xy=[]
   i=0
   maxi = len(y)
   while i < maxi:
    if deepverbose: print '\n',y[i:i+2]
    if y[i:i+2]=='*[':
      for j in range(i+1,len(y)):
        if y[j] =='*':
          xy.append((i,j+1))
          i = j+1
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
   if deepverbose: print '\n',x,'\n',y,'\n mid second pass'

   # last pass, verify pairs, if not in legit substitute both with dots. 
   # It is still needed
   if deepverbose: print '\n',x,'\n',y,'\n','final pass sanitization'
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

def acodon_loop(s): #find position of anticodon loop in a structure

# instead of perfect match (find) need partial match... implement through substring count
# better would be to count pairs... or at least verify them
# add a requirement for one base margin in the loop. Experimental!

   good3= False
   good5= False
   good = False
   for scope in range(26,8,-1): # Find left arm of anticodon stem: arbitrary start from 26 down to 8 
    substructure = s[scope:scope+5]
    cleft = substructure.count('<') #This is based on infernal structures so "<" is OK.
#    print cleft,substructure
    if cleft >= 2: # 3 misses some ultra-short stems...
     good3 = True
     break
#   left=scope+5 with this possibility of missing loop border exists...
                ## better identify it by finding non-bracket position...
   i=scope+cleft
   for i in range(scope+cleft,len(s)-5):
    if s[i] not in ['<','>','(',')']: ## may expand the list...should include ","? No.
     break
# left is first non-bracket position so the beginning of the loop
   left = i
# now find the end of the loop
   for scope in range(left+5,min(len(s)-5,left+15)): 
#tuning needed? So the minimum size of the loop is 5? This does not seem correct - should be 3+2+2=7...
#unless some very tight acodon loops actually do exist... need to investigate...
#And maximum is 15? 
    substructure = s[scope:scope+5]
    cright = substructure.count('>')
#    print cright,substructure
    if cright >=cleft:
     good5 = True
     break
   for i in range(scope+1,left,-1):
    if s[i] not in ['<','>','(',')']:
     break
   right = i+1
#right is the last position in the loop (or at least it should be). Or the first after (in pythonic counting).
   if (good3 and good5 and (right - left) > 5):
#in any case 3 is too liberal condition for a minimal loop size. 
#Moreover with the above such short loop wouldn't have been found. Remove or tighten.
     good  = True
     if verbose:
#      print s
      print s[left-5:right+5]
#      print s[left:right]
#      print left,right
   else:
     if verbose:
      print s
      print 'anticodon loop not identified'
## optional verification of pairing to be implemented here
#   if verbose: print s[left-1:right+1],left,right
   return left+1,right-1,good #the extra margin around stem. May need to be removed if some aa is not found.

###MAIN

if len(sys.argv) < 2:
 print 'usage: python t2gff.py fasta_filename ttable_id' 
else:

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
  f.write('\n'.join(linie[1])) ## dirty hack to cover circular cases
                               ## true sequence is in linie[1:], true length is also correct
                               ## for simplicity we store the "extended" version for slicing
 extended_seq = ''.join(linie[1:])+linie[1]
# print extended_seq

 file_path = seqid + '.cms'
 if not (os.path.isfile(file_path)):
   command = 'cmscan trn_a.cm ' + temp_file + '> ' +  file_path 
   os.system(command)
   os.remove(temp_file)
 else:
   print 'Re-using ', file_path

 relaxed = ['!','?']
 pref_code5  = {'F':'TGAAR','C':'TGCAR','Y':'TGTAR','H':'TGTGR','D':'TGTCR','I':'TGATR','N':'TGTTR','L1':'YTAAR','S1':'TTGAR','W':'TTCAR','P':'TTGGR','R':'TTCGR','Q':'TTTGR','V':'TTACR','A':'TTGCR','G':'TTCCR','E':'TTTCR','T':'TTGTR','K':'TTTTR','M':'YCATR','M2':'YTATR','L2':'YTAGR','S2':'TKCTR'}
 pref_code2  = {'F':'TGAAR','C':'TGCAR','Y':'TGTAR','H':'TGTGR','D':'TGTCR','I':'TGATR','N':'TGTTR','L1':'YTAAR','S1':'TTGAR','W':'TTCAR','P':'TTGGR','R':'TTCGR','Q':'TTTGR','V':'TTACR','A':'TTGCR','G':'TTCCR','E':'TTTCR','T':'TTGTR','K':'TTTTR','M':'YCATR','L2':'YTAGR','S2':'TRCTR','M2':'YTATR'}
 pref_code4  = {'F':'TGAAR','C':'TGCAR','Y':'TGTAR','H':'TGTGR','D':'TGTCR','I':'TDATR','N':'TGTTR','L1':'YTAAR','S1':'TTGAR','W':'TTCAR','P':'TTGGR','R1':'TTCGR','Q':'TTTGR','V':'TTACR','A':'TTGCR','G':'TTCCR','E':'TTTCR','T':'TTGTR','K':'TTTTR','M':'YCATR','R2':'YYCTR','L2':'YTAGR','S2':'TRCTR'}
 pref_code9  = {'F':'TGAAR','C':'TGCAR','Y':'TGTAR','H':'TGTGR','D':'TGTCR','I':'YDATR','N':'TDTTR','L1':'YTAAR','S1':'TTGAR','W':'TTCAR','P':'TTGGR','R':'TTCGR','Q':'TTTGR','V':'TTACR','A':'TTGCR','G':'TTCCR','E':'TYTCR','T':'TTGTR','K':'TCTTR','M':'YCATR','L2':'YTAGR','S2':'TKCTR'}
 pref_code13 = {'F':'TGAAR','C':'TGCAR','Y':'TGTAR','H':'TGTGR','D':'TGTCR','I':'TGATR','N':'TGTTR','L1':'YTAAR','S1':'TTGAR','W':'TTCAR','P':'TTGGR','R':'TTCGR','Q':'TTTGR','V':'TTACR','A':'TTGCR','G1':'TTCCR','E':'TTTCR','T':'TTGTR','K':'TTTTR','M':'YCATR','M2':'YTATR','L2':'YTAGR','S2':'TRCTR','G2':'YYCTR'}

 if ttable_id == '2':
   pref_code = pref_code2
 elif ttable_id == '4':
   pref_code = pref_code4
 elif ttable_id == '9':
   pref_code = pref_code9
 elif ttable_id == '13':
   pref_code = pref_code13
 else:
   pref_code = pref_code5
 
 with open(file_path, 'r') as file_input:
	infernal_data = file_input.read()
	data = infernal_data.split("Hit alignments:",1)
        genes = data[1].split("\n>> ")

 list_of_records = []

# parse the main output from cms file into memory, store all important data in a list of records
# each record is also a list: 
# product(string)[0], 
# possible_aa(list of strings)[1], 
# start(string)[2], 
# end(string)[3], 
# score(string)[4], 
# possible_anticodon(list of strings)[5], 
# retain(boolean)[6], 
# strand(string)[7]
  
 for gene in genes[1:]:
  linia = gene.split('\n')
  column = linia[3].split()
  if column[1] in relaxed: 
## cmscan stores reversed start/end coordinates for '-' strand
## better to correct this immediately, make sure it is consistent
   start = column[9]
   end = column[10]
   score = column[3]
   strand = column[11]
   if strand == '-':
    start = column[10]
    end = column[9]
   product = linia[0].strip()
## now parse structure and sequence (aln)
   structure = linia[6].rstrip(' CS')
   aln = linia[9].upper().replace('U','T')
   short_structure = structure.lstrip()
   offset = len(structure) - len(short_structure)
   structure = short_structure
   aln = aln[offset:len(structure)+offset]
## need conversion to common coordinates/lengths
   if deepverbose: print '\nstructure:',structure,'\n      aln:',aln
   (better_structure, better_aln) = sanitize_structure(structure,aln) 
   structure = better_structure
   aln = better_aln
   if 'N' in aln.upper():
     #indicative of excessive sanitization, replace with true sequence
     true_aln = extended_seq[int(start)-1:int(end)]
     if strand == '-':
       true_aln = reverse_complement(true_aln)
     aln = true_aln.upper()
   aln_len = int(end) - int(start) +1 
   if len(aln) == aln_len == len(structure):
     if deepverbose: print len(aln),len(structure)
     mfe = RNA.energy_of_structure(aln.replace('T','U'),structure.replace('<','(').replace('>',')'),0)
   elif int(start) <= sequence_len: #final sanity check, should never happen
     if  verbose: print "can't adopt cmscan structure. Re-folding at", start, end
     (structure, mfe) = RNA.fold(aln)
   score = -1*mfe
#   if mfe > 0: 
#     print 'adjusting score'
#     score = 0
   if verbose:
    print '' 
    print product
    print structure,start,end
    print aln, score

#now guess anticodon and possible_aa based on pref code table (dictionary)
   left,right,not_bad = acodon_loop(structure)
   possible_aa = []
   possible_anticodon = []
   possible_offset = []
   possible_structure = []
   no_left_arm= 'S' in product
   if not_bad: 
# bad case happens if anticodon loop is not found, record will not be stored
# try to minimize those by sophisticating detection
    if verbose: print 'not bad',
    for aa in pref_code:
     hits = nt_search(aln[left:right], pref_code[aa])
     if len(hits) >1:
      possible_aa.append(aa)
      offset = int(hits[1]) + left
#      offset = int(hits[1]) + left + 1
      ##this is not used internally so must be in non-pythonic coordinates! Add "1" for pythonic god!
      possible_offset.append(offset)
      possible_anticodon.append(aln[offset+1:offset+4])
      if verbose: print aa,aln[offset+1:offset+4],

      #ended up with 3 lists in the same order: 
      #one with possible aa [0]
      #the second with respective (real) anticodons [4], 
      #third with their relative positions [6]
      #make separate record for each valid aa+anticodon combination
      #avoiding the need for sanitization of product
      # the original product is entirely ignored... 

   if len(possible_aa) > 0:
    for i in range(0, len(possible_aa)):
      dane_do_zapisu = [possible_aa[i], start, end, score, possible_anticodon[i], strand, possible_offset[i], structure]
      list_of_records.append(dane_do_zapisu)
      if verbose: print '\naccepted', possible_aa[i],
   else:
     if verbose: print 'no aa...', aln,aln[left:right],pref_code

#filter out nearly identical hits: similar location, same aa and anticodon
# there is something wrong with this procedure... It leaves true identical records... Why?
# same record in checked and accepted gets duplicated somehow?
# to prevent this from happening the results should be stored in a dictionary structure, not a list
# keyed by the defining elements: 0, 1,2,3(?),4,5... possibly 6 (same as in filtering step)

 if verbose: print '\n',len(list_of_records),'records before filter'
 if do_filter:
  delta=10 ##This needs tuning...here orelse during clustering this should approach half of the usual trn length... 
  checked = []
  accepted = []
  while len(list_of_records) > 0:
   rec1=list_of_records.pop()
   best =True
   for rec2 in list_of_records+checked:
    ovlap =  ((abs(int(rec1[1]) - int(rec2[1])) < delta) and (abs(int(rec1[2]) - int(rec2[2])) < delta))
    same = (rec1[0]==rec2[0]) and (rec1[4] == rec2[4])
    if ovlap and same and float(rec1[3]) < float(rec2[3]):## which location is retained?...
      if deepverbose: print 'duplicate',rec1[0:5]
      best = False
      break
   if best and not (rec1 in accepted):
      accepted.append(rec1)
   checked.append(rec1)
 list_of_records = accepted

 if verbose: print len(list_of_records), 'in the middle'
# accepted = [list_of_records.pop()]
 accepted = []
 while len(list_of_records) > 0:
  r = list_of_records.pop(0)
  duplicate = False
  for q in accepted:
   if r[0:6]==q[0:6]:
     duplicate =True
     break
  if not duplicate:
     accepted.append(r)

 list_of_records = sorted(accepted, key = lambda x: float(x[3]), reverse = True)
 if verbose: print len(list_of_records),'after filter'
 with open(seqid + ".gff", 'a') as gff:
  for rec in list_of_records:
     strand = rec[5]
     start = rec[1]
     end = rec[2]
     aux = 'product=tRNA-' + seq3(rec[0][0]) + ';anticodon=(' + rec[4].lower() + ');anticodon_position=' + str(rec[6]) +';label=' + rec[0] + ';structure='+rec[7]
     linia_output = seqid + '\t' + 'infernal\ttRNA' + '\t' + start + '\t' + end + '\t' + str(rec[3]) + '\t' + strand + '\t.\t' + aux + '\n'
     if int(start) <= sequence_len:     
         gff.write(linia_output)

