import sys
import re
import os
import os.path
from glob import glob
from Bio import SeqIO, SeqUtils, pairwise2, AlignIO
from Bio.Alphabet import IUPAC, generic_rna
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Graphics import GenomeDiagram 
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from Bio.Seq import Seq
from itertools import cycle
from reportlab.lib import colors
import more_itertools
import copy
import math
from colorama import init
#from collections import OrderedDict
#from taxon_names_resolver import Resolver, TaxDict
#from ete3 import NCBITaxa
#import requests
#import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
#import matplotlib.patches as patches
import warnings
from scipy import stats
from numpy import median,std
import RNA
from math import copysign
from collections import Counter
from datetime import date

warnings.filterwarnings("ignore")

init()
print ''

deepverbose = True
deepverbose = False
verbose = False
#verbose = True

#aa_list = ['C','T','F','Y','H','D','I','N','L1','L2','S1','S2','W','P','R','Q','V','A','G','E','K','M','M2']
mitochondria =['ATP6','COX1','COX2','COX3','CYTB','ND1','ND2','ND3','ND4','ND5','ND6','ND4L','ATP8']
OXPHOS = {'I':['ND1','ND2','ND3','ND4','ND5','ND6','ND4L'],
          'III':['CYTB'],
          'IV':['COX1','COX2','COX3'],
          'V':['ATP6','ATP8']}

aa_prop = {'C':'pol','T':'pol','F':'aro','Y':'aro','H':'bas','D':'aci','I':'ali','N':'pol','L':'ali','S':'pol','W':'aro','P':'pro','R':'bas','Q':'pol','V':'ali','A':'ali','G':'gly','E':'aci','K':'bas','M':'ali','*':'non','X':'any','B':'any','Z':'any','J':'ali'}
length_limit ={'ATP6':600,'COX1':1300,'COX2':500,'COX3':700,'CYTB':1000,'ND1':800,'ND2':1000,'ND3':300,'ND4':1000,'ND5':1500,'ND6':400,'ND4L':250,'ATP8':200,'ORF':130,'orf':85,'small rRNA subunit (12S)':800, '12S':800,'large rRNA subunit (16S)':1000,'16S':1000}

## ORF limit should be more flexible, could be used for interactive tuning. 
#  Since this is in dictionary - could be modified accordingly. This is the easiest way to get rid of bogus ORFs overlapping true genes

#TO DO: tune per-gene length_limit based on real data!!! Needed in containment filtering, may be used in Oyala - alternatively. 
#       Determines acceptable limit of tRNA-non_tRNA overlap 
#       Used also in CDS filtering limiting the acceptable length of overlap

overlapl = lambda f1, f2:  any(x in f1 for x in [f2.start, f2.end + (-1)]) or any(x in f2 for x in [f1.start, f1.end + (-1)]) 
# overlap of simple locations

total_c = lambda f1, f2: all(x in f2 for x in f1)
# containment of locations (f1 in f2)

def safe_divide(x,y):
 try:
  z=float(x)/y
 except:
  z=0
 return z


def overlap (f1, f2): # overlap of features, even with both compound locations
 for el1 in f1.location.parts:
  for el2 in f2.location.parts:
   c = overlapl(el1,el2)
   if c: break    
  if c: break
 return c

def contain(f1,f2):
 return total_c(f1.location,f2.location) 

# virtualization: combine 2-part compound location at ends into one in virtual coordinates
# for evaluation of offsetted location containment/overlap

def virtualize(loc): 
  p0=loc.parts[0]
  p1=loc.parts[-1]
  direct = loc.strand
  vs = p0.start
  ve = p1.end + sequence_len
  if direct == -1:
   vs = p1.start
   ve = p0.end + sequence_len
  return FeatureLocation(vs,ve,direct)

#need separate case for de-virtualization of simple location (shift)
#could be checked before calling? If both ends are > sequence_len...
# could use offset_loc with delta = 0 for such cases...

def de_virtualize(loc):
 direct = loc.strand
 true_end = loc.end - sequence_len
 seg1 = FeatureLocation(loc.start, sequence_len, strand=direct)
 seg2 = FeatureLocation(0, true_end, strand=direct)
 if direct == -1:
          location = CompoundLocation([seg2, seg1])
 else:
          location = seg1 + seg2
 return location

#switch two locations to commmon virtual space

def pair_virt(loc_f1,loc_f2):
 if len(loc_f1.parts) > 1:
   loc_f1 = virtualize(loc_f1)
   if loc_f2.start + sequence_len <= loc_f1.end:
    loc_f2 = loc_f2 + sequence_len 
 if len(loc_f2.parts) > 1:
   loc_f2 = virtualize(loc_f2)
   if loc_f1.start + sequence_len <= loc_f2.end:
    loc_f1 = loc_f1 + sequence_len
 if loc_f1.start > sequence_len and loc_f2.start > sequence_len:
  loc_f1 = loc_f1 + sequence_len*(-1)
  loc_f2 = loc_f2 + sequence_len*(-1)
 return loc_f1,loc_f2 

def offset_loc(loc,delta):
 ## equality condition in de_virt criterion was problematic
 ## changed allowing ends at sequence_len
 ## pos=sequence_len is NOT in, and so is 'end'
 x = loc
 if len(x.parts)>1: ## WRONG condition!! breaks truly compound locations!!! Need alternative procedure for such cases!
  x = virtualize(x)
 x = x + delta
 if x.start < 0:
  x = x + sequence_len
 if x.start >= sequence_len:
  x = x + sequence_len*(-1)
 if x.end > sequence_len:
  x = de_virtualize(x)
 return x

def offset_loc_compound(loc,delta):
 loc = loc + delta
 x=[]
 for y in loc.parts:
  if y.start < 0:
   y = y + sequence_len
  if y.end > sequence_len and y.start < sequence_len:
   y =  de_virtualize(y)
  if y.start >= sequence_len:
   y = y + sequence_len*(-1)
  for i in y.parts:
    x.append(i)
  if len(x) >1:
   z = CompoundLocation(x)
  else:
   z = x[0]
 return z

def offset_universal(loc, delta):
 size = len(loc.parts)
 strand = loc.strand
 if size == 1:
   return offset_loc(loc,delta)
 else:
   raw_list = []
   for i in loc.parts:
    raw_list.append(i)
   new_list = []
   j=0
   for i in range(0,size-1):
       j+=1
       if j == size:
        j=0
       if raw_list[i].start == 0 and raw_list[j].end == sequence_len and strand == -1:
         intermediate = FeatureLocation(raw_list[j].start,raw_list[i].end+sequence_len,strand)
         new_list.append(intermediate)
       elif raw_list[i].end == sequence_len and raw_list[j].start == 0 and strand <> -1:
         intermediate = FeatureLocation(raw_list[i].start,raw_list[j].end+sequence_len,strand)
         new_list.append(intermediate)
       else:
         new_list.append(raw_list[i])
         new_list.append(raw_list[j])
   raw_list = new_list
   if deepverbose:
     print'before shift:',size
     for x in raw_list:
      print x
   modified_list = []
   for i in raw_list:
    shifted = i + delta
    if shifted.start >= sequence_len:
      shifted = shifted + sequence_len*(-1)
    if shifted.end > sequence_len:
      shifted = de_virtualize(shifted)
      modified_list.extend(shifted.parts)
    else:
      modified_list.append(shifted)
   if len(modified_list) > 1:
     final_loc = CompoundLocation(modified_list)
   else:
     final_loc = modified_list[0]
   if deepverbose :
     print 'after shift:'
     for x in modified_list:
      print x
   return final_loc 
    
def diminute_end(loc):
 x = loc
 if len(x.parts)>1:
  x = virtualize(x)
 xe = x.end - 1
 x = FeatureLocation(x.start, xe, x.strand)
 if x.end > sequence_len and x.start < sequence_len:
  x =  de_virtualize(x)
 if x.start >= sequence_len:
  x = x + sequence_len*(-1)
 return x

def diminute_start(loc):
 x = loc
 if len(x.parts)>1:
  x = virtualize(x)
 xs = x.start + 1
 x = FeatureLocation(xs, x.end, x.strand)
 if x.end > sequence_len and x.start < sequence_len:
  x =  de_virtualize(x)
 if x.start >= sequence_len:
  x = x + sequence_len*(-1)
 return x

def offset_contain(f1,f2,delta):
 loc = f1.location
 fx = SeqFeature(offset_loc(loc,delta),'misc',{})
 return contain(fx,f2)

def offset_overlap(f1,f2,delta):
    loc = f1.location
    fx = SeqFeature(offset_loc(loc,delta),'misc',{})
    return overlap(fx,f2)

def big_score(f1,f2):
   score = 0
   for pos in f1.location:
    if pos in f2:
     score +=1
   return score

def in_frame(f1,f2):
 c1 = (f1.type == 'CDS') and (f2.type == 'CDS')
 c2 = (f1.strand == f2.strand)
 v1,v2 = pair_virt(f1.location,f2.location)
 if f1.strand == 1:
  c3 = (v1.start % 3) == (v2.start % 3)
 if f1.strand == -1:
  c3 = (v1.end %3) == (v2.end %3)
 return (c1 and c2 and c3)

def midpoint_BAD(loc):
 x = loc
 if len(loc.parts) > 1:
  x = virtualize(loc)
 y = x.start + (x.end - x.start) / 2
 if y > sequence_len:
  y = y - sequence_len
 return y

def midpoint(loc): # for simple locations only, really
 x = loc
 if len(loc.parts) > 1:
  x = virtualize(loc)
 y = x.start + (x.end - x.start) / 2
 if y >= sequence_len:
  y = y - sequence_len
 return y

def ending(f):
  if f.strand == -1:
   x = f.location.parts[-1].start
  else:
   x = f.location.parts[0].start
  return x
 
# returns the slice of a record based on a feature.
# returns letter_annotations, not only sequence
# !!! BUGFIX!!! returned wrong part for compound locations if STRAND unaccounted for!!!
# not sure, re-check...
# will not work correctlty for truly compound locations (ea. set of features)
# for that separate slices based on each feature would have to be created.
# CAUTION! Features crossing location parts are NOT preserved! No easy work-arund exists.

def slice_rec(rec,fx):
  if fx.strand <> -1:
   sub_rec_start = fx.location.parts[0].start
   sub_rec_end = fx.location.parts[-1].end
  else:
   sub_rec_start = fx.location.parts[-1].start
   sub_rec_end = fx.location.parts[0].end
  if sub_rec_start > sub_rec_end:
   sub_rec = rec[sub_rec_start:] + rec[:sub_rec_end]
  else:
   sub_rec = rec[sub_rec_start:sub_rec_end]
  return sub_rec

def window_set(rec, wdw): #returns a list of records, being sliced in 1100 (+/-1) overlappig windows
 values = []
 stepik = max(sequence_len/1100,1) #tune
 if stepik > wdw:
     stepik=wdw-4
 delta = -1* (wdw)/2
 initial_loc = FeatureLocation(0,wdw-stepik,strand=+1)
 working_loc = FeatureLocation(0,wdw,strand=+1)
 initial_loc = offset_loc(initial_loc, delta)
 working_loc = offset_loc(working_loc, delta)
 while midpoint(working_loc) > 0:
#  print midpoint(working_loc),' ',
  initial_loc = offset_loc(initial_loc,-1)
  working_loc = offset_loc(working_loc,-1)
 f0 = SeqFeature(initial_loc,'misc',{})
 fx = SeqFeature(working_loc,'misc',{})
 while True:
  sub_rec = slice_rec(rec,fx)
  coordinate = midpoint(fx.location)
  sub_rec.id=str(coordinate)
  sub_rec.description =''
  values.append(sub_rec)
  fx.location = offset_loc(fx.location,stepik)
  if contain(f0, fx):
     sub_rec = slice_rec(rec,fx)
     coordinate = midpoint(fx.location)
     sub_rec.id=str(coordinate)
     sub_rec.description =''
     values.append(sub_rec)
     delta = sequence_len-1-coordinate
     if delta > stepik/2:
      fx.location = offset_loc(fx.location,delta)
      sub_rec = slice_rec(rec,fx)
      coordinate = midpoint(fx.location)
      sub_rec.id=str(coordinate)
      sub_rec.description =''
      values.append(sub_rec)
     else:
      del values[-1]
      fx.location = offset_loc(fx.location,delta)
      sub_rec = slice_rec(rec,fx)
      coordinate = midpoint(fx.location)
      sub_rec.id=str(coordinate)
      sub_rec.description =''
      values.append(sub_rec)
     break
 return values

#Consider switching all windowed procedures to the above common procedure generating windowed sub_records
# more memory consuming but should be OK for mitogenomes

## Caution: must use record slices to propagate letter_annotations
def window_skew(rec, i, wdw, L1, L2, stat):
 values = []
 stepik = max(sequence_len/1100,1) #tune
 initial_loc = FeatureLocation(0,wdw-stepik,strand=+1)
 delta = -1* (wdw)/2
 initial_loc = offset_loc(initial_loc, delta)
 f0 = SeqFeature(initial_loc,'misc',{})
 working_loc = FeatureLocation(0,wdw,strand=+1)
 working_loc = offset_loc(working_loc, delta)
 fx = SeqFeature(working_loc,'misc',{})
 while True:
  sub_rec = slice_rec(rec,fx)
  s = sub_rec.seq.upper()
  x = 0
  y = 0
  for pos in range(len(sub_rec)):
      if sub_rec.letter_annotations[stat][pos] in i:
          if s[pos] == L1: x+=1
          if s[pos] == L2: y+=1
  try:
   skew = (x-y) / float(x+y)
  except ZeroDivisionError:
   skew = 0.0
  coordinate = midpoint(fx.location)
  values.append((coordinate,skew))
  fx.location = offset_loc(fx.location,stepik)
  if contain(f0, fx):
   values = sorted(values)
   avg_point = (values[0][1]+values[-1][1])/2
   values.append((sequence_len-1, avg_point))
   values.append((0, avg_point))
   break
 return values

def AT_skew(rec, wdw=200):
 return window_skew(rec, [0,1,2,3], wdw, 'A', 'T', 'codon')

def scaling_factor(VL):
 m = 0
 for a in VL:
   for x,y in a:
     b = abs(y)
     if b > m:
       m = b
 return m

def normalize_values(value_list,max_value,ofset=0):
  new_list = []
  for x,y in value_list:
    y = (y-ofset) / max_value
    new_list.append((x,y))
#  avg_point = (new_list[0][1]+new_list[-1][1])/2.0
#avg_point should be calculated more accurately than simple average
# this is really sub-pixel resolution issue
#but if the gap is smaller than half pixel the second augment is not needed?
# This is a BUG! One more record should be calculated, with the center at zero!
#nope, this is calculated perfectly well. The last one (missing) should always be very close to the first one.
#  new_list = [(0, avg_point)]+new_list
#  new_list.append((sequence_len-1, avg_point))
#  new_list.append((sequence_len-1, new_list[0][1]))
  return new_list

 ### for calculation of filtered contents stats of a feature (and record)
 ### NEED letter annotations!
 ### May be used in sliding window, individual features or whole records
 ### CAN'T be used on sets of features, unless separate letter annotation is created first
   # i: list of letter annotation values
   # fx: feature
   # r: record, needed for sequence and letter annotations
   # L: list of counted letters (uppercase!)
   # stat: type of letter annotation used (codon, prop etc..)

# simplified version with counting only
def r_con(r,i,L,stat):
 s = r.seq.upper()
 XL=dict.fromkeys(L,0)
 denom = 0
 for pos in range(len(r)):
      if r.letter_annotations[stat][pos] in i:
          denom +=1
          for y in L:
           if y == s[pos]:
             XL[y]+=1
 return sum(XL.values()), denom

def r_contents(r,i,L,stat):
 numer,denom = r_con(r,i,L,stat)
 try:
   z = numer / float(denom)
 except ZeroDivisionError:
   z = 0.0
 return z 

def f_contents(f, r, i, L, stat):
 sub_rec = slice_rec(r,f)
 return r_contents(sub_rec,i,L,stat)

## the old windowed procedures are obsolete; all windowed stats are done through window_set
## review and remove/rewrite depending on needs
def window_contents(rec, i, wdw, L1, L2, stat):
 values = []
 stepik = max(sequence_len/1100,1) #tune
 initial_loc = FeatureLocation(0,wdw-stepik,strand=+1)
 delta = -1* (wdw)/2
 initial_loc = offset_loc(initial_loc, delta)
 f0 = SeqFeature(initial_loc,'misc',{})
 working_loc = FeatureLocation(0,wdw,strand=+1)
 working_loc = offset_loc(working_loc, delta)
 fx = SeqFeature(working_loc,'misc',{})
 while True:
  contents = f_contents(fx,rec,i,[L1,L2],stat)
  coordinate = midpoint(fx.location)
  values.append((coordinate,contents))
  fx.location = offset_loc(fx.location,stepik)
  if contain(f0, fx):
   values = sorted(values)
   avg_point = (values[0][1]+values[-1][1])/2
   values.append((sequence_len, avg_point))
   values.append((0, avg_point))
   break
 return values

def GC_contents(rec, wdw=200):
 return window_contents(rec,[0,1,2,3],wdw,'G','C','codon')

def GCn(rec, wdw=200):
 return window_contents(rec,[0,1],wdw,'G','C','prop')

def GC3(rec, wdw=200):
 return window_contents(rec,[3],wdw,'G','C','codon')

def loc2gb(loc):
    from StringIO import StringIO
    out_handle = StringIO()
    f = SeqFeature(loc,'misc',{})
    r = SeqRecord(rec.seq)
    r.features.append(f)
    out_handle.write(r.format('gb'))
    result = out_handle.getvalue().split('misc')[1].split('ORIGIN')[0].strip()
    if deepverbose: print 'converting', f.location, ' to ', result
    return result

def gb2loc(location_text): #expected input format: pos:gb_location
 a=location_text.split(':')
 s=1
 loc = FeatureLocation(0,sequence_len, strand=s)
 if a[0] == 'pos':
#    print a[1][:10]
    if a[1][:10] == 'complement' :
      #use negative strand, strip envelope:
#      print 'reverse'
      s = -1
      a[1] = a[1][11:-1]
    p = a[1].split('..')
    if len(p) == 1:
      p.append(p[0])
    loc = FeatureLocation(int(p[0])-1,int(p[1]),strand=s)
 return loc

def acodon_pos(g):
        from StringIO import StringIO
        out_handle = StringIO()
#        deepverbose = True
        if deepverbose: print 'estimating anticodon location for',g.qualifiers['product']
        rel_pos = int(g.qualifiers['anticodon_position'])
        while True: ## space for condition if rel_pos too low to consider valid.
                    ## for reasoably good candidates this is obsolete...
         anticodon = g.qualifiers['anticodon'][-4:-1].upper()
         g_sequence=g.extract(rec.seq)
         subseq = g_sequence[rel_pos:].upper()
         pos = subseq.find(anticodon)
         if pos == -1: ## consider simplifying since rel_pos is now reliable..
          rel_pos = rel_pos - 1
          if verbose: print 'correcting acodon position...'
          if rel_pos < 10: 
           return 'unknown',0
         else:
          if g.location.strand == 1:
           if len(g.location.parts)>1:
             ofset = virtualize(g.location).start
           else: 
             ofset = g.location.start
#           pos = subseq.find(anticodon)
           abs_pos = ofset + pos  + rel_pos
           if deepverbose: print ofset,pos,rel_pos,'ofset, pos, rel_pos: components for anticodon location (+)'
           abs_loc = FeatureLocation(abs_pos, abs_pos + 3, strand = +1)
           if deepverbose: print ' temporary anticodon location created: abs_loc:', abs_loc
           break
          if g.location.strand == -1:
           if len(g.location.parts)>1:
             ofset = virtualize(g.location).end
           else: 
             ofset = g.location.end
           pos = subseq.find(anticodon)
           abs_pos = ofset - pos - rel_pos
           if deepverbose: print ofset,pos,rel_pos,'components for anticodon location (-)'
           abs_loc = FeatureLocation(abs_pos-3, abs_pos, strand = -1)
           if deepverbose: print ' temporary anticodon location created', abs_loc
           break
        if abs_loc.start > sequence_len:
          abs_loc = abs_loc + sequence_len*(-1)
        if abs_loc.start < 0:
          abs_loc = abs_loc + sequence_len
        if abs_loc.end > sequence_len:
          abs_loc = de_virtualize(abs_loc)
        if deepverbose: print ' estimated anticodon position:', abs_loc
#        deepverbose = False
        return loc2gb(abs_loc),rel_pos+pos

def flip_struct(x):
        br_dict={'(':')',')':'(','>':'<','<':'>','[':']',']':'[','{':'}','}':'{','.':'.','_':'_',':':':',',':',','-':'-','~':'~'}
        span = len(x)
#        print x, span
        r = x[::-1]
#        print r, span
        for j in range(0,span):
#           print j,r[j],br_dict[r[j]]
           r = r[:j]+br_dict[r[j]]+r[j+1:]
        return r ## reversed structure with adjusted brackets.

def annotate_trna(g):
         new_location = g.location 
         new_qual = {}
         abs_pos,rel_pos = acodon_pos(g)
#         rel_pos = int(g.qualifiers['anticodon_position'])
# rel_pos in raw annotations is unreliable!!!
# use abs_pos to get it here!
# However, better fix this in t2gff to prevent problems
# after that may finally simplify acodon_pos
         ss = g.qualifiers['structure']
         ss = ss.replace('_','.')
         ss = ss[:rel_pos]+'___'+ss[rel_pos+3:]
         new_feature = SeqFeature(new_location, 'tRNA', new_qual)
         new_feature.qualifiers['anticodon'] = '(pos:' + abs_pos + ',aa:' + g.qualifiers['product'][-3:] + ',seq:' + g.qualifiers['anticodon'][-4:-1] + ')'
         new_feature.qualifiers['product'] = g.qualifiers['product']
         new_feature.qualifiers[g.qualifiers['source']] = g.qualifiers['score']
         new_feature.qualifiers['label'] = g.qualifiers['label']
         new_feature.qualifiers['structure'] = ss##g.qualifiers['structure']
         if g.location.strand == -1:
            new_feature.qualifiers['structure'] = flip_struct(ss)
         return new_feature

def annotate_cds(CDS_candidates):

#modify with sorting, maybe dictionarizing?
# need to have best possible name (from Genewise), best possible location (from getorf) for starters.
    if verbose: print 'combining orfs at locations'
#    f = CDS_candidates.pop(0)
    max_l = 0
#    f_score = 
    max_score = 0
    new_name = 'ORF'
#    new_location = f.location
#    new_qual['evidence']=''
#    max_l = 1
    #just a proposal for location and name. For refinement as follows.
#    if verbose: print f.location
    for g in CDS_candidates:
      if verbose: print g.location, g.qualifiers['source'], g.qualifiers['product']
      new_score = float(g.qualifiers['score'])/len(g) #THIS IS THE MOST PROBLEMATIC part
#      new_qual['evidence'] = new_qual['evidence'] + ' ' + g.qualifiers['source'] #mainly for debugging, not really needed
#adopt the name of the best scorer. Genewise should allways win. No longer true.???
# alternatively just adopt the Genewise first, then best scorer... Without the re-calculation, possibly?
# after sanitization there should be one Genewise hit at most... Unless its chimera, but these should be deatt with earlier.
# so the whole "score" business is not that important, if there is no GeneWise hit the neame could be anything.
# NO NEED TO CALCULATE OR STORE THE "score", trun off "downgrading"
      if new_score > max_score:
#        new_qual.update({g.qualifiers['source']:g.qualifiers['score']})
        max_score = new_score
#        new_name = g.qualifiers['product']
      new_qual.update({g.qualifiers['source']:g.qualifiers['score']})
      if g.qualifiers['source'] == 'GeneWise':
#        max_score = new_score
        new_name = g.qualifiers['product']
      if len(g) > max_l and g.qualifiers['source'] == 'getorf':
        max_l = len(g)
        new_location = copy.copy(g.location) #adopt the longest orf
                                  #make sure it _is_ ORF! 
                                  #genewise can give longer but bogus annotations so should be disregarded here
                                  # in fact length is not important at all, there should be only one getorf feature..
                                  ##TO BE implemented and TESTED        

#    if max_score < 0.06: 
   #simple hack for bogus names. Perhaps this should be left to a later stage, 
   # with similar accounting as with trn. Evidently one max_score is not good for all cases.
   ### less important after initial sanitization but still: possible truncation of good CDS
   ### BETTER lineage-specific HMMs => increase of the threshold! <-This is the correct solution.
   #just the name is changed, the CDS will still be there with all scores
   #add original name to notes
#      new_qual['note']='downgraded from '+ new_name
#      new_name = 'ORF'
    new_feature = SeqFeature(new_location, 'CDS', {})
    new_feature.qualifiers.update(new_qual)
    new_feature.qualifiers['product'] = new_name
    new_feature.qualifiers['score'] = max_score     
    if verbose: print new_feature.location, new_feature.qualifiers['product'],len(new_feature)
    return new_feature

def finalize(loc,f):
  if deepverbose: print 'finalizing ', loc,
  if loc.end > sequence_len:
   loc = offset_loc(loc,0)
  fx = SeqFeature(loc,f.type,{})
  fx.qualifiers.update(f.qualifiers)
  if deepverbose: print ' into ', fx.location
  return fx


def over5(f,g): #first is allways cutting, second is being cut
 if verbose: print "5'", f.qualifiers['product'], g.qualifiers['product'], ' at ', g.location
 codon_OK = 1
 need = True
 if len(f.location.parts) == len(g.location.parts) == 1: 
  r = f.location
  c = g.location
 else:
  r,c = pair_virt(f.location,g.location)
 c_end = c.end
 c_start = c.start
 delta = -2*f.strand ##allow 2bp tolerance
# c0 = c
 if verbose: print ' in common space: ', r, c
 while overlapl((r+delta),c) and need:
    while True:
### TO DO: fix this! the procedure with changing delta was wrong, did not stop after finding proper start codon!!!
###        would make CDS shorter than they need to be. Consider re-writing and simplifying.
     if f.strand == -1:
      c_end = c_end - 3
     if f.strand == 1:
      c_start = c_start + 3
     c = FeatureLocation(c_start, c_end, c.strand)
     gx = SeqFeature(c,g.type,{})
     if gx.extract(rec.seq+rec.seq)[0:3].upper() in start_codons:
      if verbose: print '.',
#      delta = 0
      codon_OK = 1
      break
     else:
      codon_OK = 0
      if verbose: print '-',
      if len(gx) < 79: # more robust criterion of failure needed? Perhaps use per gene limits already here?
        if verbose: print 'too short'
#        c = c0    
        codon_OK = 1
        need = False
        break
 if need:
#  fx = finalize(r,f)
  gx = finalize(c,g)
  if verbose:
    for delta_base in range(3):
      delta = -1*delta_base*f.strand
      if not overlapl((r+delta),c):
       break 
    print g.qualifiers['product'], ' adjusted with tolerance ', delta, gx.location
 else: 
#  fx = f
  if verbose: print 'adjustment impossible'
  gx = g
 return gx

def over3(f,g): #first should always be cutting (ea. rna), second cut (ea. cds)
 if verbose:  print "3'", f.qualifiers['product'], g.qualifiers['product']
 OK=0
 if len(f.location.parts) == len(g.location.parts) == 1: 
  r = f.location
  c = g.location
  if deepverbose: print 'simple: ', r, c
 else:
  r,c = pair_virt(f.location,g.location)
  if deepverbose:  print 'complex, in common space:', r, c
 c_end = c.end
 c_start = c.start
 r_start = r.start
 r_end = r.end
 if f.strand == -1: 
  for c_start in [r_end, r_end -1, r_end +1, r_end-2,r_end+2]:
   if (c_start > c_end) or (c_start < r_start): 
    OK = 0
    if verbose: print 'topology'
    break  
   c = FeatureLocation(c_start,c_end,c.strand)
   r = FeatureLocation(r_start,c_start,r.strand)
   phase = len(c) % 3
   last = (rec+rec)[c.start:c.start+phase].reverse_complement().seq + 'AAA'
   if last[:3].upper() in stop_codons:
     OK = 1
     if verbose: print 'OK', last[:3]
     break
 if f.strand == 1:
  for c_end in [r_start, r_start +1, r_start -1,r_start+2, r_start-2]:
   if (c_end < c_start) or (c_end > r_end): 
    OK = 0
    if verbose: print 'topology'
    break  
   c = FeatureLocation(c_start,c_end, c.strand)
   r = FeatureLocation(c_end, r_end, r.strand)
   phase = len(c) % 3
   last = (rec+rec)[c.end - phase:c.end].seq + 'AAA'
   if last[:3].upper() in stop_codons:
     OK = 1
     if verbose: print 'OK', last[:3]
     break
 if OK:
#  fx = finalize(r,f) ## is this needed?
  gx = finalize(c,g)
 else:
  if verbose: print 'failed, no suitable stop codon in the vicinity(+/-2) of cutting gene'
#  fx = f
  gx = g
##       special cases where no cds cutting is needed 
##       because the stop codon is already there
##       Really, this is only for tbl2asn compatibility
 if OK and big_score(f,g) == 2:
   trn_seq = fx.extract(rec.seq)
   cds_seq = gx.extract(rec.seq)
   if cds_seq[-1].upper() + trn_seq[:1].upper() in stop_codons:
    gx = g
    if verbose: print 'two bp overlap, existing stop, restoring CDS'
 if OK and big_score(f,g) == 1:
   trn_seq = fx.extract(rec.seq)
   cds_seq = gx.extract(rec.seq)
   if cds_seq[-2:].upper() + trn_seq[0].upper() in stop_codons:
    gx = g
    if verbose: print 'one bp overlap, existing stop, restoring CDS'
 if OK and big_score(f,g) == 3:
   trn_seq = f.extract(rec.seq)
   if trn_seq[0:2].upper() in stop_codons:
    gx = g
    if verbose: print 'three bp overlap, existing stop, restoring CDS'
 return gx

# first cutting second cut. Version with walking up the cutting point, for resolving chimeras.
# may be used for sorting of overlaps but with caution: problem of phase difference!
# better use over3 for that!
def under3(f,g): 
 if verbose:  print "3'", f.qualifiers['product'], g.qualifiers['product']
 OK=False
 if len(f.location.parts) == len(g.location.parts) == 1: 
  r = f.location
  c = g.location
  if deepverbose: print 'simple: ', r, c
 else:
  r,c = pair_virt(f.location,g.location)
  if deepverbose:  print 'complex, in common space:', r, c
 c_end = c.end
 c_start = c.start
 r_start = r.start
 r_end = r.end
 if c.strand == 1:
   c_end = r_start
 if c.strand == -1:
   c_start = r_start
 while True:
     c = FeatureLocation(c_start, c_end, c.strand)
     for phase in [1,2]:
      if f.strand == -1:
       last = (rec+rec)[c.start + phase : c.start + 3].reverse_complement().seq + 'AAA'
      else:
       last = (rec+rec)[c.end - 3 : c.end - phase].seq + 'AAA'      
      if last[:3].upper() in stop_codons:
        OK = True
        if verbose: print 'OK', last[:3], c, phase
        break
     if OK:
       if f.strand == 1:
         c_end = c_end - phase
       if f.strand == -1:
         c_start = c_start + phase
       c = FeatureLocation(c_start, c_end, c.strand)
       break
     if f.strand == 1:
       c_end = c_end - 3
     if f.strand == -1:
       c_start = c_start + 3
     if deepverbose: print '-',last[:3]
     if (c_end < c_start): 
       OK = False
       if verbose: print 'topology'
       break
 if len(c) < 79: # more robust criterion of failure needed?
        if verbose: print 'too short'
        OK = False
 if OK:
  gx = finalize(c,g)
 else:
  if verbose: print 'failed, no suitable stop codon upstream'
  gx = g
 return gx,OK

def Oyala(rna,cds):
     if verbose: print 'entering Oy'
     delta = 1
     rx = rna
     cx = cds
     while delta < len(cds):
        if offset_overlap(rna,cds,delta) and not offset_overlap(rna,cds,delta*(-1)):
             if cds.strand == 1:
                 cx = over5(rna,cds)
                 if verbose: print 'over5 on +', delta  
             else:
                 cx = over3(rna,cds)
                 if verbose: print 'over3 on -',delta
             break
        if offset_overlap(rna,cds,delta*(-1)) and not offset_overlap(rna,cds,delta):
             if cds.strand == 1:
                 cx = over3(rna,cds)
                 if verbose: print 'over3 on -',delta
             else:
                 cx = over5(rna,cds)
                 if verbose: print 'over5 on +',delta
             break
        delta = delta + 1
     return cx

# takes a list of dictionaries and returns new dictionary with sum of values
# no empty keys, no side effects through modified DL members (?)
def dsum(DL):
 def merge(d1,d2):
  return dict(d1,**d2)
 WL=copy.deepcopy(DL)
 d = reduce(merge,WL)
 for k in d:
  d[k]=0
  for x in DL:
   d[k]+=x.get(k,0)
 return d

# subtracts values in d2 from values in d1, returns new d1
# assumes that d2 keys are all present in d1
# keys with zero counts are removed

def dsub(d1,d2):
 d=copy.deepcopy(d1)
 zero=[]
 for k in d2:
  try:
   d[k]-=d2[k]
   if d[k]==0:
    zero.append(k)
  except KeyError:
    continue
 for k in zero:
  del d[k]
 return d

# for summing up per-record annotations:
# makes dictionary keyed by given annotations, values are lists 
# non-existing annoations are not represented (no empty lists)
def get_stats(RL):
 good_keys=['utable','ttable','GCn','GCnTM']
 total={}
 for r in RL:
  for k in good_keys:
   try:
    total[k].append(r.annotations[k])
   except:
    total[k]=[r.annotations[k]]
 return total

def f2list(FL):
 loc_list = []
 for f in FL:
  for loc in f.location.parts:
   loc_list.append(loc)
 return sorted(loc_list, key = lambda x: x.start) # converts list of features to a list of simple locations, sorted by location start

def amalgamate_list(loc_list): # for amalgamation of repeats,
                               # takes a (redundant) list of simple locations as an argument, returns a list of non-overlapping simple locations
                               # side effect: the locations in the returned list have None strand, regardless of the source!
                               # but the source is not affected
                               #for general-purpouse amalgamation need slightly different procedure:
                               #  -> amalgamated location must have the end calculated as max of the two ends
   loc_list_final = []
   loc_list = sorted(loc_list, key = lambda x: x.start, reverse = False)
   loc = loc_list.pop(0)
   loc= FeatureLocation(loc.start,loc.end, None)
   
   while len(loc_list) >0:
    loc1 = loc_list.pop(0)
    loc1 = FeatureLocation(loc1.start,loc1.end, None)
    if loc1.start <= loc.end:
      loc = FeatureLocation(loc.start,loc1.end, None)
    else:
      loc_list_final.append(loc)
      loc= FeatureLocation(loc1.start,loc1.end, None)
      
   if not (loc in loc_list_final):
       loc_list_final.append(loc)
   return loc_list_final

def amalgamate_zero(FL): # for amalgamation of features, should be much faster
                         # version with  new features inheriting strandness
   loc_list = f2list(FL)
   return amalgamate_variable_list(loc_list)

def amalgamate_variable_list(loc_list):   
  loc_list_final = []
  loc_list = sorted(loc_list, key = lambda x: x.start, reverse = False)
  if len(loc_list) >0: 
   loc = loc_list.pop(0)
   loc= FeatureLocation(loc.start,loc.end, loc.strand)
   
   while len(loc_list) >0:
    loc1 = loc_list.pop(0)
    loc1 = FeatureLocation(loc1.start,loc1.end, loc1.strand)
    if overlapl(loc,loc1):
      loc = FeatureLocation(loc.start,max(loc1.end,loc.end), loc.strand)
    else:
      loc_list_final.append(loc)
      loc= FeatureLocation(loc1.start,loc1.end, loc1.strand)
      
   if not (loc in loc_list_final):
       loc_list_final.append(loc)
  return loc_list_final

def circularize(loc_list): # de-virtualization without compounding, gets list, returns list but element(s) extending past sequence_len are split.
                           # strand is copied, order of elements is not checked
                           #important update: de-compound compound locations first (list2list)!
    final_list = []
    loc_list_copy = list2list(loc_list)
    while len(loc_list_copy) > 0:
         x = loc_list_copy.pop(0)
         if x.end > sequence_len:
            x1 = FeatureLocation(x.start,sequence_len, x.strand)
            x2 = FeatureLocation(0,x.end - sequence_len,x.strand)
            final_list.extend([x1,x2])
         else:
            final_list.append(x)
    return final_list


def amalgamate(FL, tolerance=100): ##FL is a list of features, [f1,f2,..fn]
 def combine(x,y):
  z = FeatureLocation(x.start,y.end,x.strand)
  return z
 def distance(x,y):
  dist = y.start - x.end
  if dist < 0:
    dist = dist + sequence_len
  return dist 
 list_parts = amalgamate_zero(FL)
 if len(list_parts) > 1:  
  while True:
   list_parts = sorted(list_parts, key = lambda x: x.start)
   part_no = len(list_parts)
   parts_left = part_no
   if deepverbose: print 'starting with ',part_no 
   if part_no > 1:
    #now do the combining
    for i in range(part_no-1, 1, -1):
      if deepverbose: print i-1, list_parts[i-1], i, list_parts[i]
      if distance(list_parts[i-1],list_parts[i]) < tolerance+1:
        loc = combine(list_parts[i-1], list_parts[i])
        del list_parts[i]
        del list_parts[i-1]
        list_parts = list_parts + [loc]
        parts_left = len(list_parts)
        if deepverbose: print parts_left, ' parts left'
        if i < parts_left: break
   if parts_left == part_no: break
   #separate check for the amalgamation of first and last
   #with propper de-virtualization 
  if distance(list_parts[-1], list_parts[0]) < tolerance+1:
   loc = combine(list_parts[-1], list_parts[0] + sequence_len)
   del list_parts[-1]
   del list_parts[0]
   loc = de_virtualize(loc)
   list_parts.append(loc)
 if deepverbose: print 'ending with ', len(list_parts)
 return list_parts   ##returns a list of locations, gaps are included with arbitrary tolerance, the last one may be compound

## produces a list of locations present in f1 but not in f2.
## iterates over the range covered  by both to decide.
## how to prevent structure? Maybe use record slicing?

def subtract(f1,f2): ##f1 - f2 -> list of locations, some may be compound. Potentally very slow. Avoid.
 XXspacer = lambda x: x == 'XX'
 f1_loc,f2_loc = pair_virt(f1.location,f2.location)
 start = min(f1_loc.start, f2_loc.start)
 stop = max(f1_loc.end, f2_loc.end)
 join_pos = []
 for pos in range(start,stop):
  if (pos in f1) and (pos not in f2):
   join_pos = join_pos + [pos]
  else:
   join_pos = join_pos + ['XX']
 locations = list(more_itertools.split_at(join_pos, XXspacer))
 list_parts = []
 if f1.strand == f2.strand:
   strand = f1.strand
 else:
   strand = None
 for i in locations:
  if len(i) > 0:
   s_start = i[0]
   s_end = i[-1]+1
   loc = FeatureLocation(s_start, s_end, strand)
   if loc.start > sequence_len:
     loc = loc + sequence_len*(-1)
   if loc.end > sequence_len:
     loc = de_virtualize(loc)
   list_parts.append(loc)
 return list_parts

def raw_ncr_locations(rec):
 xxx = lambda x: x == 'XX'
 ncpos =[]
 for pos in range(0,len(rec)):
  if not any(pos in f for f in rec.features):
   ncpos = ncpos + [pos]
  else:
   ncpos = ncpos + ['XX']
 ncrs_locations = list(more_itertools.split_at(ncpos, xxx))
 list_ncrs = []
 for i in ncrs_locations:
  if len(i) > 0:
   s_start = i[0]
   s_end = i[-1]+1
   ncr_loc = FeatureLocation(s_start, s_end, None) #ncrs don't have a strand
   list_ncrs.append(ncr_loc)
 list_ncrs = sorted(list_ncrs, key = lambda x: x.start)
 return list_ncrs #returns a sorted (by start) list of simple locations not assigned to any feature

def collect_ncrs(rec):
 list_ncrs = raw_ncr_locations(rec)
 first_ncr = list_ncrs[0]
 last_ncr = list_ncrs[-1]
 if first_ncr.start == 0 and last_ncr.end == sequence_len and len(list_ncrs) > 1:
   seg1 = last_ncr
   seg2 = first_ncr
   if verbose: print 'removing', last_ncr, first_ncr
   list_ncrs.remove(last_ncr)
   list_ncrs.remove(first_ncr)
   ncr_loc = seg1 + seg2
   list_ncrs.append(ncr_loc)
 list_features = []
 for loc in list_ncrs:
  ncr = SeqFeature(loc,'misc',{})
  list_features.append(ncr)
 if deepverbose:
    print "list of ncr FEATUREs' locations:"
    for f in list_features:
     print f.location
 return list_features #returns list of features, not locations. The last feature may have CompoudLocation

def list2list(LL): # gets a list of locations, some may be compound
  loc_list = []
  for x in LL:
   for loc in x.parts:
    loc_list.append(loc)
  return loc_list # returns list of simple locations 

def local_dominant(f,fp,fm): #compares overlap of f with two references - for local strandness
      dominant = False
      current_strand = None
      op = overlap(f,fp)
      if deepverbose: print op, fp.strand
      om = overlap(f,fm)
      if deepverbose: print om, fm.strand
      if op and not om:
        dominant = True
        current_strand = fp.strand
      if om and not op:
        dominant = True
        current_strand = fm.strand
      return dominant, current_strand

def parse_gff(gff_name,rec):
 length = len(rec)
 with open(gff_name, 'r') as gff_input:
  gff_file_contents = gff_input.read()
  gff_lista_linii = gff_file_contents.split("\n")
  for raw_gff in gff_lista_linii:
   gff = raw_gff.split()
   if len(gff) > 8:
    start = int(gff[3])
    end = int(gff[4])
    f_type = gff[2]
    direct = None
    if gff[6] == '+': 
      direct = 1
    elif gff[6] == '-':
      direct = -1
    qual = {}
    key_val_list = gff[8].split(';')
    for key_val in key_val_list:
     key_list = key_val.split('=',1)
     qual.update({key_list[0]:key_list[1]})
    qual.update({'source':gff[1], 'score':gff[5]})
    location = FeatureLocation(start-1, end, strand=direct)
    if location.start < 0: #this is unusual... but should be supported(?). Alternative: stop with error.
      location = location + length
    if location.end > length: 
      location = de_virtualize(location) #this means that already at gff parsing compund locations are possible.
    new_feature = SeqFeature(location, type=f_type, qualifiers=qual)
    rec.features.append(new_feature)
 return rec 

def parse_gff3(gff_name,rec): #used in prot_rec annotation
# length = len(rec.seq)
 with open(gff_name, 'r') as gff_input:
  gff_file_contents = gff_input.read()
  gff_lista_linii = gff_file_contents.split("\n")
  for raw_gff in gff_lista_linii:
   gff = raw_gff.split('\t')
   if len(gff) > 8:
    start = int(gff[3])
    end = int(gff[4])
    f_type = gff[2]
    direct = None
    if gff[6] == '+': 
      direct = 1
    elif gff[6] == '-':
      direct = -1
    qual = {}
    key_val_list = gff[8].split(';')
    for key_val in key_val_list:
     key_list = key_val.split('=',1)
     if len(key_list)>1:
      qual.update({key_list[0]:key_list[1]})
    qual.update({'source':gff[1], 'score':gff[5]})
#    if f_type == 'ORF':
#       s_nt = start - 1
    if f_type == 'protein_match':
#     if qual['Name'] == 'TRANSMEMBRANE':
#       x_nt = (start-1)*3 + s_nt
#       y_nt = end*3 + s_nt
       loc = FeatureLocation(start-1,end,strand=direct)
#       if loc.end > length: 
#         loc = de_virtualize(loc)
       clean_feature = SeqFeature(loc, f_type, {})
       clean_feature.qualifiers = qual
#       clean_feature.qualifiers['source'] = qual['source']
#       clean_feature.qualifiers['note'] = qual['Name'].lower()
       rec.features.append(clean_feature)
 return rec 

def parse_bed(bed_name,rec): #for optional(?) annotaions in bed format
 with open(bed_name, 'r') as bed_input:
   bed_contents = bed_input.read()
   bed_lines = bed_contents.split('\n')
   for raw_bed in bed_lines:
     bed= raw_bed.split('\t')
     if len(bed) > 3: # only meaningful entries, please
       start=int(bed[1])
       end=int(bed[2])
#       f_type=' '.join(bed[3].split('_')[0:2]) # tuned for polyA from ContextMapper
#       f_type='polyA site'
       f_type=bed[3].split('_')[0] #not really iomportant, will be overwritten
       score=int(bed[4])
       direct = None
       if bed[5] == '+':
         direct=1
       elif bed[5] == '-':
         direct=-1
       if start < end: #basic syntax requirement, if not - report faulty bed file
         loc=FeatureLocation(start,end,strand=direct)
         if loc.end > len(rec): # for cases of bed from concatamers - to account for circular maps
          loc=de_virtualize(loc)
         clean_feature=SeqFeature(loc,f_type,{})
         clean_feature.qualifiers['score']=score #store as int to allow easy subsequent filtering
         clean_feature.qualifiers['source']='bed'
         rec.features.append(clean_feature)
       else:
         if verbose: print 'bed file parsing problem'
 return rec

def parse_csv(csv_file):
 v_list=[]
 with open(csv_file,'r') as fh:
  for line in fh:
   value=max(float(line.split('\t')[1]),1.0) #dirty hack for zero-coverage, needed for log scale
   v_list.append(value)
 return v_list

def codon_usage(r): # build codon usage table, using record r, based on ALL CDS features
                    # for each gene/group (GOI) separately 
                    # will need to be called with record annoatated accordingly (ea. with GOI only)
 t ={}
# if verbose: print 'cux start', len(r)
 for f in r.features:
  if f.type == 'CDS':
#   if verbose: print len(f), len(f.qualifiers)
   dna = f.extract(r.seq)
   maxlen = 3*(len(dna)/3) #this should round the length down to full codons.
   for i in range(0,maxlen,3):
    codon = str(dna[i:i+3]).upper()
    if all(x in ['A','C','T','G'] for x in list(codon)): #unambiguity condition
     try:
      t[codon]+=1
     except KeyError:
      t[codon] =1
 if deepverbose:  print " codons counted {0:5d} record length {1:5d} codon types {2:5d} features:{3:5d}".format(sum(t.values()),len(r),len(t),len(r.features))
 return t # return a dictionary with codons as keys and counts as values
          # this table is genetic code neutral and does not have zeroes, only represented codons are counted!
          # can be sorted by codon (key x[0]) or number (x[1]) but this is usually not needed
          # augmenting with ttable dependent aa is trivial through codon translation

def ECNa(t): #effective number of codons - alternative implementation
             # aat: dictionary of dictionaries, keyed by aa+first letter of codon
             # 
             # acount: total number of aa (==codons) scored in t
             # the stat should be between len(aat) and len(t),
             # to normalize (20 - 61)
             # use simple aritmetic linear transformation 
 an = 0
 if True: #len(t) > 20:
  aat ={}
  for x in t:
#x is key so codon; t[x] is codon count. x[0] is first nt of a codon
#aa is derived by translation of a codon but aa families are distinguished by addition of first codon nt
   dna = Seq(x,IUPAC.ambiguous_dna)
   aa = str(dna.translate(table=ttable))+x[0]
   try:
    aat[aa].update({x:t[x]})
   except KeyError:
    aat[aa] = {x:t[x]}

# builds a dictionary containing aa-family counts
# after update the value should 
# store all counts, keyed by codon
# efficiently it is just clustering of of counts

  for aa in aat:
   n = 0
   ccount = sum(aat[aa].values()) # total codon count, per aa-family
   for cdn in aat[aa]:
    nn = aat[aa][cdn] # individual codon count
    n += nn*(nn) # need only sum of squares
   try:
    an += (float(ccount)*ccount) / n
   except:
    an += 1 #only case with n=1 should trigger this
  try:
   minn = len(aat) #number of aa (families) present
   maxn = len(t) #number of codons present
   nn = (an-minn)/(maxn-minn)
   an=(41*nn)+20
  except:
   an=20
 return an

def ECN(t): #calculate effective number of codons based on the dictionary produced by codon_usage
 aacount = {}
 aap2 = {}
 n = 0.0
 utable = CodonTable.unambiguous_dna_by_id[ttable_id]
 ### xc[0]: codon; xcs: string representation of codon; xcn: count; xc[1]: amino acid 
 for xc in utable.forward_table.iteritems():
  xcs = str(xc[0])
  xcn = t.get(xcs,0)
  try:
   aacount[xc[1]] += xcn
  except KeyError:
   aacount[xc[1]] = xcn
#  print xc[1], aacount[xc[1]],'count'
 for xc in utable.forward_table.iteritems():
  xcs = str(xc[0])
  xcn = t.get(xcs,0)
  ## make sure no division by zero occurs!
  p = 0
  if aacount[xc[1]]>0:
   p = float(xcn) / aacount[xc[1]]
  p = p*p
  if p > 0:
   try:
    aap2[xc[1]] += p
   except KeyError:
    aap2[xc[1]] = p
 for aa in aap2:
  n += (1 / aap2[aa])
##NORMALIZATION: experimental
 maxn=len(utable.forward_table)
 minn=len(aap2)
 if minn <20:
   maxn = 0
   for xc in utable.forward_table.iteritems():
    try:
     p = aap2[xc[1]]
     maxn += 1
    except:
     continue
# print n,minn,maxn
# try:
#  nn = (n-minn)/(maxn-minn)
# print nn
#  n=(41*nn)+20
# except:
#  n=0
##ALTERNATIVE, "harmonic" normalization TESTING needed
 minn = float(minn)
 maxn = float(maxn)
 try:
  n = 1/(0.05-(0.0336/((1/minn-1/maxn))*(1/minn -1/n)))
 except:
  n = 0.0
 return n            

def ECNw(t): #effective number of codons - closer to Wrigtht 
 aacount = {}
 aap2 = {}
 n = 0.0
 utable = CodonTable.unambiguous_dna_by_id[ttable_id]
 for xc in utable.forward_table.iteritems():
  xcs = str(xc[0])
  xcn = t.get(xcs,0)
  try:
   aacount[xc[1]] += xcn
  except KeyError:
   aacount[xc[1]] = xcn
 for xc in utable.forward_table.iteritems():
  xcs = str(xc[0])
  xcn = t.get(xcs,0)
  ## make sure no division by zero occurs!
  p = 0
  if aacount[xc[1]]==1:
   p = 1
  elif aacount[xc[1]]>1:
   p = float(xcn*(xcn-1)) / float(aacount[xc[1]]*(aacount[xc[1]]-1))
  if p > 0:
   try:
    aap2[xc[1]] += p
   except KeyError:
    aap2[xc[1]] = p
 for aa in aap2:
  n += (1 / aap2[aa])
##NORMALIZATION: experimental
 maxn=len(utable.forward_table)
 minn=len(aap2)
 if minn <20:
   maxn = 0
   for xc in utable.forward_table.iteritems():
    try:
     p = aap2[xc[1]]
     maxn += 1
    except:
     continue
# print n,minn,maxn
# try:
#  nn = (n-minn)/(maxn-minn)
# print nn
#  n=(41*nn)+20
# except:
#  n=0
 minn = float(minn)
 maxn = float(maxn)
 try:
  n = 1/(0.05-(0.0336/((1/minn-1/maxn))*(1/minn -1/n)))
 except:
  n = 0.0
 return n            

def ECG(t): #total (not average per aa) effective number of codons  
 ccount = 0
 n = 0
 for x in t:
   ccount += t[x]
   n += t[x]*(t[x]-1)
 try:
  gn = float(ccount-1)*ccount/n
 except:
  gn = 0.0
 return gn             

def EAA(t): #effective number of encoded amino acids from codon table
 acount =0
 n = 0
 aat ={}
 for x in t:
   dna = Seq(x,IUPAC.ambiguous_dna)
   aa = dna.translate(table=ttable)
   try:
    aat[aa] +=1
   except KeyError:
    aat[aa] = 1
 for aa in aat:
  acount +=aat[aa]
  n += aat[aa]*(aat[aa]-1)
 try:
  an = float(acount-1)*acount/n
 except:
  an = 0.0
 return an

def shannon(orf):
#Shannon entropy of a sequence
 orf=orf.upper()
 span=len(orf)
 denom=float(span)
 f={}
 shannon=0.0
 f['A']=orf.count('A')/denom
 f['C']=orf.count('C')/denom
 f['T']=orf.count('T')/denom
 f['G']=orf.count('G')/denom
 for p in f:
  try:
   shannon -=f[p]*math.log(f[p]) 
  except:
   pass
 return shannon

def CDSr(rec): ## extract all coding sequences of rec to a single multi-fasta file
 fh = open(name[0]+'.CDS.fasta','w')
 for f in rec.features:
  if f.type in ['CDS','rRNA']:
   dna=f.extract(rec.seq)
   r = SeqRecord(dna)
   r.id = f.qualifiers['product']
   r.description =''
   SeqIO.write(r,fh,'fasta')
 fh.close()
 return

def trnas(rec): ## extract all trna sequences of rec to a single multi-fasta file
 fh = open(name[0]+'.trn.fasta','w')
 for f in rec.features:
  if f.type in ['tRNA']:
   dna=f.extract(rec.seq)
   r = SeqRecord(dna)
   r.id = f.qualifiers['label']
   r.description = ''
   try:
    r.description = f.qualifiers['structure']
   except:
    pass
   SeqIO.write(r,fh,'fasta')
 fh.close()
 return

def onscreen_cds(r):
 cdss = []
 cdss_len = 0
 for f in r.features:
  if f.type == 'CDS':
   cdss.append(f)
   cdss_len += len(f)
 gene_count = len(cdss)

 missing=list(mitochondria)
 for x in cdss:
  if x.qualifiers['product'] in missing:
   missing.remove(x.qualifiers['product'])

 cdss = sorted(cdss, key= lambda x: midpoint(x.location))
 if len(missing)>0: print '\033[91m',','.join(missing),
 else: print '\033[92m',
 print gene_count, 
 text_map= '\033[0m CDS [ '
 for x in cdss:
  if x.qualifiers['product'] in mitochondria: text_map += '\033[92m'
  text_map += x.qualifiers['product']+' \033[0m'
 print text_map+']',cdss_len/3, 'codons'
 return

def onscreen(r): #nice onscreen reporting. Does not need to return any value
 onscreen_cds(r)
 aa_in = []
 trns={}
 rrnl = 0
 for f in r.features:
    if f.type == 'rRNA':
         rrnl += len(f)
    if f.type == 'tRNA':
         aminoacid = f.qualifiers['label']
         aa = aminoacid + '-' + f.qualifiers['anticodon'][-4:-1]
         aa_in.append(aminoacid)
         trns.update({aa:f})
 no_genes = len(trns)

 aa_desired = aa_list
 for x in aa_in+['M2']:
    if x in aa_desired:
         aa_desired.remove(x)

 if len(aa_desired) > 0: print '\033[91m',','.join(aa_desired),
 else: print '\033[92m',
 print no_genes, '\033[0m tRNA types'
 if rrnl < 1500:  print '\033[91m',
 else: print '\033[92m',
 print rrnl, '\033[0m bp in rRNA genes'
 return

def sanitize_structure(x,y): #x structure, y sequence
   opening =  ['(','[','{','<']
   closing =  [')',']','}','>']
   legit = ['AT','TA','GC','CG','GT','TG']
   last = len(x)-1
   i=-1
   while i < last:
      i = i + 1
      bcount = 0
      if x[i] in closing:
        bcount -=1
        for j in range(i-1,0,-1):
         if x[j] in opening:
          bcount+=1
         elif x[j] in closing:
          bcount-=1
         if bcount == 0:
          if not (y[i]+y[j] in legit):
           x = x[:i]+'.'+x[i+1:]
           x = x[:j]+'.'+x[j+1:]
          break
        if bcount == 1:
          x = x[:i]+'.'+x[i+1:]
          x = x[:j]+'.'+x[j+1:]
      elif x[i] in opening:
        bcount +=1
        for j in range(i+1,len(x)):
         if x[j] in opening:
          bcount+=1
         elif x[j] in closing:
          bcount-=1
         if bcount == 0:
          if not (y[i]+y[j] in legit):
           x = x[:i]+'.'+x[i+1:]
           x = x[:j]+'.'+x[j+1:]
          break
        if bcount == 1:
          x = x[:i]+'.'+x[i+1:]
          x = x[:j]+'.'+x[j+1:]
   return x


###
## Main program starts here
###

if os.name == 'nt':
    redir ='NUL'
    copy_str = 'copy'
else:
    redir ='/dev/null'
    copy_str = 'cp'

argument = sys.argv[1]
name = argument.split(".",1)
fasta_name = name[0]+'.fa'
output_name = name[0] + ".gb"
gff_name = name[0] + ".gff"
raw_name = name[0] + "_raw.gb"
#clus_name = name[0] + "_clus.gb"
tbl_name = name[0] + ".tbl"
fsa_name = name[0] + '.fsa'
val_name = name[0]+'.val'
gbf_name = name[0]+'.gbf'
gbk_name = name[0]+'.gbk'
sqn_name = name[0]+'.sqn'
pdf_name = name[0]+'.pdf'
prot_name = name[0] +'_prot.gb'
clean_pdf_name = name[0]+'_clean.pdf'
taxid_name = name[0]+'.taxid'

rec = SeqIO.read(fasta_name, "fasta")
rec.seq.alphabet = IUPAC.ambiguous_dna
sequence_len = len(rec.seq)
rec.features = []
a = rec.seq.count('A') + rec.seq.count('a')
t = rec.seq.count('T') + rec.seq.count('t')
c = rec.seq.count('C') + rec.seq.count('c')
g = rec.seq.count('G') + rec.seq.count('g')
at_c = float(a+t)
skewt = (a-t)/at_c
gc_c = float(g+c)
skewc = (g-c)/gc_c
rec.annotations.update({'topology':'circular'})
# assuming circular topology, gb annotation is unimportant but all analytical tasks assume circularity.
# there is no way to turn it off globally - each algo would have to be taken care of separately,
# which is impractical.

### TEST circular offsetting ### 
# remove after done#
#loc=FeatureLocation(0,3,1)
#while True:
# loc=offset_loc(loc,1)
# print midpoint(loc),' ',
#quit()
###

rec.annotations['data_file_division']='INV'
rec.annotations['date']=date.today().strftime("%d-%b-%Y").upper()

print '\033[92m'+ name[0]+':\033[0m', sequence_len, 'bp,', "{:.1f}".format(100*at_c/sequence_len)+'% AT, AT-skew:', "{:.3f}".format(skewt),'GC-skew:',"{:.3f}".format(skewc)

# single level identification by blasting local nucleotide database with the whole sequence
# lineage and binominal latin name of the best match are recovered from local sources
# latin name of the first match is used as provisional taxid, this is needed to choose lineage-specific  hmms 
# this used to work well but is apparently throwing wrong identifiers on some Arthropods and flatworms
# need to investigate.
# possible cause: highly biased sequences without known close relatives
# Need fix!!! Manual supply of taxid file is the only option for now.

if not (os.path.isfile(taxid_name)):
#  command = 'blastn -task dc-megablast -db '+os.environ['DATABASES']+'/Metazoa -query '+fasta_name+' -evalue 0.01 -perc_identity 80 -culling_limit 1 -word_size 12 -max_target_seqs 5 -outfmt "6 sseqid" > ' + taxid_name 
  command = 'blastn -db '+os.environ['DATABASES']+'/Metazoa -query '+fasta_name+' -evalue 0.01 -max_target_seqs 5 -outfmt "6 sseqid" > ' + taxid_name 
  if verbose: print command
  os.system(command)
try:
 lineage =[]
 bh = open(taxid_name, 'r')
 blast_result = bh.read()
 bh.close()
 blines = blast_result.split("\n")
 blast_answer = blines[0].split('_')
 taxid = ' '.join(blast_answer[:-1]) #binominal name
 bh = open('Metazoa.lineages','r')  #list of lineages, first column is latin name (or gb equivalent of it)
 #now scan the file for match of the first column, when found retrive the rest as lineage
 for current_line in bh:
  tax_map = current_line.strip('\n').split(';',1)
  if tax_map[0] == taxid:
   lineage = tax_map[1].split(',')
   break
 if len(lineage) <3:
   raise ValueError()
except: #what could go wrong? Assuming that taxid file is a problem... So empty taxid should also trigger this...
 print "unidentifiable sequence, caution advised"
 lineage=['unknown']

if verbose: print 'preparing reference for critica..'
try: 
 specific_db = lineage[2]
except:
 specific_db = 'Metazoa'

# it is probably better to remove the existing symlinks first, then create a new set in a loop
#for n in glob('critica.*'):
#     os.remove(n)
# dont't do that if db is not specific
#print specific_db
if specific_db== 'Metazoa':
 if verbose: print "problem identifying specific lineage"
else:
 command = 'blastn -task megablast -db '+os.environ['DATABASES']+'/'+ specific_db + ' -query '+fasta_name+' -max_target_seqs 5 -outfmt 6 > tmp.blast'
 os.system(command)
 with open('tmp.blast', 'r') as bh:
  line = bh.readline().strip()
  col =line.split()
# the following condition is provisional. Consider more sophistication and do make cd-hit cleanup occasionally.
  if not (int(float(col[2])) == 100 and int(col[3])>0.8*sequence_len):
   if verbose: print 'extending reference',specific_db, col
   dbh = open(os.environ['DATABASES']+'/'+specific_db+'.fa','a')
   dbh.write(rec.format("fasta"))
   dbh.close()
   command = 'makeblastdb -in '+ os.environ['DATABASES']+'/'+specific_db+'.fa -dbtype nucl -blastdb_version 4 -out '+os.environ['DATABASES']+'/'+specific_db
   os.system(command)
  else:
   if verbose: print 'database',specific_db,'unchanged'
 bh.close()

# symlinking stuff here
# need to re-create symlinks each time (?)
# for critica to work with current mod: it expects blast database in the current directory, named "critica"
# so we need symlinks pointing at all parts of the "specific_db" and named "critica"+appropriate extension.

for ext in ['.nhr','.nin','.nsq']:#,'.ndb','.not','.ntf','.nto']: 
     os.remove('critica'+ext)

for ext in ['.nhr','.nin','.nsq']:#,'.ndb','.not','.ntf','.nto']: 
# blast db version 4 uses three indexes and works, symlinking does not work for version 5
 try:
  os.symlink(os.environ['DATABASES']+'/'+ specific_db + ext, 'critica'+ext)
 except:
  continue
if len(lineage)>3:
 print 'Similar to',taxid, 'from: \033[33m', lineage[-3], '\033[0m'

aa_list =  ['C','T','F','Y','H','D','I','N','L1','L2','S1','W','P','R','Q','S2','V','A','G','E','K','M','M2']
ttable_id = 5
if 'Vertebrata' in lineage:
    ttable_id = 2
    rec.annotations['data_file_division']='VRT'
if any(x in lineage for x in ['Echinodermata','Platyhelminthes','Hemichordata']):
    ttable_id = 9
    aa_list =  ['C','T','F','Y','H','D','I','N','L1','L2','S1','W','P','R','Q','S2','V','A','G','E','K','M']
if any(x in lineage for x in ['Cnidaria', 'Ctenophora','Placozoa','Porifera']): ##!##
    ttable_id = 4
    aa_list =  ['C','T','F','Y','H','D','I','N','L1','L2','S1','W','P','R1','Q','S2','V','A','G','E','K','M','R2']
if any(x in lineage for x in ['Tunicata']):
    ttable_id = 13
    aa_list =  ['C','T','F','Y','H','D','I','N','L1','L2','S1','W','P','R','Q','S2','V','A','G1','E','K','M','G2','M2']
if any(x in lineage for x in ['Catenulida']):#,'Rhabditophora']):
    aa_list =  ['C','T','F','Y','H','D','I','N','L1','L2','S1','W','P','R','Q','S2','V','A','G','E','K','M','M2']
    ttable_id = 5

#if any(x in lineage for x in ['Echinostomatidae','Trichobilharzia']):
#    ttable = 14
#if any(x in lineage for x in ['Echinostoma','Alaria']):
#    ttable = 21
#if 'Rhabdopleuridae' in lineage:
#    ttable = 24
## Adjustments needed, testing !!!
if 'unknown' in lineage:
    ttable_id = 1
print 'Using genetic code table #', ttable_id

ttable = CodonTable.ambiguous_dna_by_id[ttable_id]
stop_codons = ttable.stop_codons
start_codons = ttable.start_codons
ttable_str = str(ttable_id)

if skewc < 0 and not ('Vertebrata' in lineage):
 print 'Negative GC skew! Plotting reverse-complement could be more appropriate.'
 rv = SeqRecord(rec.seq.reverse_complement())
 SeqIO.write(rv,name[0]+'_rev.fa', "fasta")
 print 'reversed sequence saved in',name[0]+'_rev.fa, consider re-running with this input.' 

if not (os.path.isfile(gff_name)):
### Special or not-so-special hmm still need testing for:
### ATP6 for Tridacna - Cardioidea
### COX2 for Mya - Euheterodonta
### ND6 for Dreissena - Euheterodonta
### TO DO: set of hmms for exremelly AT-rich Aphid-like insects
###        more special cases for diverse groups like flatworms
   special_rna = False
   special_ATP8 = False
   special_ATP6 = False
   special_ND2 = False
   special_ND6 = False
   special_COX2 = False
   special_ND4L = False


   if 'Vertebrata' in lineage:
    print 'will be using conservative ATP8 profile for Vertebrata..'
    os.remove('ATP8.hmm')
    command = copy_str+' ATP8_VRT.hmm ATP8.hmm'
    os.system(command)
    special_ATP8 = True

   if 'Euheterodonta' in lineage:
    print 'will be using specific COX2 and ND6 profiles for Euheterodonta..'
    os.remove('COX2.hmm')
    os.remove('ND6.hmm')
    command = copy_str+' COX2_EUH.hmm COX2.hmm'
    os.system(command)    
    command = copy_str+' ND6_EUH.hmm ND6.hmm'
    os.system(command)
    special_COX2 = True
    special_ND6 = True

   if 'Cardioidea' in lineage:
    print 'will be using specific ATP6 profile for Cardioidea..'
    os.remove('ATP6.hmm')
    command = copy_str+' ATP6_CAR.hmm ATP6.hmm'
    os.system(command)
    special_ATP6 = True

   if 'Mytilidae' in lineage:
    print 'will be using fast-evolving ATP8 profile for Mytilidae..'
    os.remove('ATP8.hmm')
    command = copy_str+' ATP8_MYT.hmm ATP8.hmm'
    os.system(command)
    special_ATP8 = True

   if 'Ostreoida' in lineage:
    print 'will be using fast-evolving ATP8, ND2 & ND6  profiles for Ostreoida..'
    os.remove('ND2.hmm')
    os.remove('ND6.hmm')
    os.remove('ATP8.hmm')
    command = copy_str+' ND6_OST.hmm ND6.hmm'
    os.system(command)    
    command = copy_str+' ND2_OST.hmm ND2.hmm'
    os.system(command)
    command = copy_str+' ATP8_OST.hmm ATP8.hmm'
    os.system(command)
    special_ATP8 = True
    special_ND2 = True
    special_ND6 = True

   if 'Pterioidea' in lineage:
    print 'will be using fast-evolving ATP6, ND2 & ND6  profiles for Pterioidea.. and oyster ATP8'
    os.remove('ND2.hmm')
    os.remove('ND6.hmm')
    os.remove('ATP6.hmm')
    command = copy_str+' ND6_PTE.hmm ND6.hmm'
    os.system(command)    
    command = copy_str+' ND2_PTE.hmm ND2.hmm'
    os.system(command)
    command = copy_str+' ATP6_PTE.hmm ATP6.hmm'
    os.system(command)
    command = copy_str+' ATP8_OST.hmm ATP8.hmm'
    os.system(command)
    special_ATP8 = True
    special_ATP6 = True
    special_ND2 = True
    special_ND6 = True

   if 'Arcoidea' in lineage:
    print 'will be using fast-evolving ND2 & ND6  profiles for Arcoidea..'
    os.remove('ND2.hmm')
    os.remove('ND6.hmm')
    command = copy_str+' ND6_ARC.hmm ND6.hmm'
    os.system(command)    
    command = copy_str+' ND2_ARC.hmm ND2.hmm'
    os.system(command)
    special_ND2 = True
    special_ND6 = True

   if 'Platyhelminthes' in lineage:
    print 'will be using specific ND2, ATP8 and ND4L profiles for Platyhelminthes..'
    os.remove('ND2.hmm')
##Need reliable ATP8 and ND4L hmms for this group...
#    os.remove('ATP8.hmm')
#    command = copy_str+' ATP8_PLAT.hmm ATP8.hmm'
#    os.system(command)    
    command = copy_str+' ND2_PLAT.hmm ND2.hmm'
    os.system(command)
    command = copy_str+' ND4L_PLAT.hmm ND4L.hmm'
    os.system(command)
    special_ND2 = True
#    special_ATP8 = True
    special_ND4L = True


   if 'Sternorrhyncha' in lineage:
    print 'will be using specific ND2, ND6, ND4L and rRNA profiles for Aphid-like insects..'
    os.remove('ND2.hmm')
    os.remove('ND6.hmm')
    command = copy_str+' ND2_STER.hmm ND2.hmm'
    os.system(command)
    command = copy_str+' ND6_STER.hmm ND6.hmm'
    os.system(command)
    command = copy_str+' ND4L_APH.hmm ND4L.hmm'
    os.system(command)
    special_ND2 = True
    special_ND6 = True
    special_ND4L= True
# for cm profiles the indexes are used, not cm files: i1m i f p
    os.remove('rrna.cm.i1m')
    command = copy_str+' rrna_APHID.cm.i1m rrna.cm.i1m'
    os.system(command)
    os.remove('rrna.cm.i1i')
    command = copy_str+' rrna_APHID.cm.i1i rrna.cm.i1i'
    os.system(command)
    os.remove('rrna.cm.i1f')
    command = copy_str+' rrna_APHID.cm.i1f rrna.cm.i1f'
    os.system(command)
    os.remove('rrna.cm.i1p')
    command = copy_str+' rrna_APHID.cm.i1p rrna.cm.i1p'
    os.system(command)
    special_rna = True


   if not special_COX2:
    try:
     os.remove('COX2.hmm')
    except:
     pass
    command = copy_str+' COX2_REF.hmm COX2.hmm'
    os.system(command)
   if not special_ND2:
    try:
     os.remove('ND2.hmm')
    except:
     pass
    command = copy_str+' ND2_REF.hmm ND2.hmm'
    os.system(command)
   if not special_ND4L:
    try:
     os.remove('ND4L.hmm')
    except:
     pass
    command = copy_str+' ND4L_REF.hmm ND4L.hmm'
    os.system(command)
   if not special_ND6:
    try:
     os.remove('ND6.hmm')
    except:
     pass
    command = copy_str+' ND6_REF.hmm ND6.hmm'
    os.system(command)
   if not special_ATP8:
    try:
     os.remove('ATP8.hmm')
    except:
     pass
    command = copy_str+' ATP8_3.hmm ATP8.hmm'
    os.system(command)
   if not special_ATP6:
    try:
     os.remove('ATP6.hmm')
    except:
     pass    
    command = copy_str+' ATP6_REF.hmm ATP6.hmm'
    os.system(command)

### implement alternative (not addidional) rrna cm here
   if not special_rna:
    try:
     os.remove('rrna.cm')
    except:
     pass
    command = copy_str+' rrna_REF.cm rrna.cm'
    os.system(command)
    os.remove('rrna.cm.i1m')
    command = copy_str+' rrna_REF.cm.i1m rrna.cm.i1m'
    os.system(command)
    os.remove('rrna.cm.i1i')
    command = copy_str+' rrna_REF.cm.i1i rrna.cm.i1i'
    os.system(command)
    os.remove('rrna.cm.i1f')
    command = copy_str+' rrna_REF.cm.i1f rrna.cm.i1f'
    os.system(command)
    os.remove('rrna.cm.i1p')
    command = copy_str+' rrna_REF.cm.i1p rrna.cm.i1p'
    os.system(command)

##  ALL analytical tasks invoked here ##
   
   print 'arwen for tRNA'
   command ='python a2gff.py '+ fasta_name +' '+ttable_str
   os.system(command)
   print 'critica for coding'
   command ='python c2gff.py '+ fasta_name + ' '+ttable_str
   os.system(command)
   print 'glimmer for cds'
   command ='python g2gff.py '+ fasta_name + ' '+ttable_str
   os.system(command)
   print 'cmscan for tRNA'
   command ='python t2gff.py '+ fasta_name + ' '+ttable_str
   os.system(command)
   print 'cmscan for rRNA'
   command ='python r2gff.py '+ fasta_name
   os.system(command)
   print 'getorf for orf'
   command ='python o2gff.py '+ fasta_name + ' '+ttable_str
   os.system(command)
   print 'wise for mitochondrial genes'
   command ='python w2gff.py '+ fasta_name + ' '+ttable_str
   os.system(command)
   print 'vmatch for repeats'
   command ='python v2gff.py '+ fasta_name
   os.system(command)

if verbose:print 'Parsing', gff_name

rec = parse_gff(gff_name,rec)

bed_name=name[0]+'.bed'
if (os.path.isfile(bed_name)):
  print 'parsing optional bed file'
  rec = parse_bed(bed_name,rec)
                     
# OK, gff contents in, features created, now can do the filtering & adjustments
# save the raw data - for debugging and reference
SeqIO.write(rec, raw_name, "genbank")
print 'raw data written to :\033[96m', raw_name,  '\033[0m'

# new record with te same sequence and selective copy of feature
new_rec = SeqRecord(rec.seq)
new_rec.annotations = rec.annotations


if verbose:print name[0],': clustering features'

cds = []
trn = []
rna = []
for f in rec.features:
    if f.type == 'CDS':
        cds.append(f)
    if f.type == 'tRNA':
        trn.append(f)
    if f.type == 'rRNA':
        rna.append(f)

list_CDS = []
while len(cds) > 0:
 f = cds.pop(0)
 CDS_candidates = [f]
 while True:
  match =[]
  for f in CDS_candidates:
   for g in cds:
       if overlap(f,g) and in_frame(f,g):
          if deepverbose: print 'in-frame overlap:', f.location, f.qualifiers['product'],g.location,g.qualifiers['product']
          match.append(g)
   for g in match:
    if g in cds:
      cds.remove(g)
    if g not in CDS_candidates:
      CDS_candidates.append(g)
    if deepverbose: print len(CDS_candidates)
  if match == []:
    break
 if verbose and len(CDS_candidates) > 1 :
     for f in CDS_candidates: print f.qualifiers['product'],
     print ' <'
 list_CDS.append(CDS_candidates)
try:
 del match
except:
 pass 

if verbose: print 'orfs in',len(list_CDS), 'groups'

list_tRNA = []
while len(trn) > 0:
 f = trn.pop(0)
 trn_candidates = [f]
 while True:
  match =[]
  for f in trn_candidates:
   for g in trn:
       if overlap(f,g) and big_score(f,g) > 40:
          if deepverbose: print 'trn overlap:', f.location, f.qualifiers['product'],g.location,g.qualifiers['product']
          match.append(g)
   for g in match:
    if g in trn:
      trn.remove(g)
    if g not in trn_candidates:
      trn_candidates.append(g)
    if deepverbose: print len(trn_candidates)
  if match == []:
    break

 trn_candidates=sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse=True)
#   trn_candidates=sorted(trn_candidates, key = lambda x: x.qualifiers['source'])
 list_tRNA.append(trn_candidates)
 if verbose:# and len(trn_candidates) > 1 :
     for f in trn_candidates: print f.qualifiers['label'],f.qualifiers['source'][0],
     print '<-',len(trn_candidates),trn_candidates[0].location

try:
 del match
except:
 pass 

if verbose: print 'tRNA candidates in',len(list_tRNA),'clusters'

list_16s = []
list_12s = []
for f in rna:
  if f.qualifiers['product'] == '12S':
    if deepverbose: print f.qualifiers['product'],len(f), len(f.qualifiers['structure'])
    list_12s.append(f)
  if f.qualifiers['product'] == '16S':
    if deepverbose: print f.qualifiers['product'],len(f), len(f.qualifiers['structure'])
    list_16s.append(f)

if verbose: print 'rRNA in',len(list_16s)+len(list_12s),'groups'

# orf sanitization DEBUGGING needed?
orf_list =[]
good_list =[]
while len(list_CDS) > 0:
  i = list_CDS.pop()
  support = {}
  if verbose and len(i)>1: print len(i),'long group'
  for f in i:
      support.update({f.qualifiers['source']:f.qualifiers['score']})
  if len(i) > len(support): # potential chimera, sanitize same source features
    if verbose: print len(support),'types represented'  
    j=[]
    if verbose:
      for f in i:
       print 'member', f.location, f.qualifiers['source'],f.qualifiers['product']
    while len(i) > 0:
     f=i.pop()
     f_is_good = True
     if verbose: print '->',f.qualifiers['product']
     for g in i+j:
      if f.qualifiers['source'] == g.qualifiers['source'] and (overlap(f,g)  and (len(f) < len(g))): 
         ## need testing/tuning but relaxing here should minimize false chimera calling and naming problems 
         f_is_good = False
         break
     if f_is_good:
       if verbose: print 'is good'
       j.append(f)
     else:
       if verbose: print 'is redundant'
    i = j
    if verbose :
     for f in i: print f.qualifiers['product'],
     print ' sanitized', len(i),len(support)
  if len(support) == 1 and i not in orf_list:
    i = sorted(i, key = lambda x: len(x), reverse = True)  
    orf_list.append(i)
  elif len(support) > 1 and i not in good_list:
    i = sorted(i, key = lambda x: x.location.start, reverse = False)
    good_list.append(i)
list_CDS = good_list

if verbose: print len(list_CDS),'good orf groups, non-clusterable orfs: ', len(orf_list)

# feature lists created. Now process them separately.

# start with CDS.

for CDS_candidates in list_CDS:
   #consider only multiple hits, the rest should be in orf_list anyway. Simplified.
   new_qual = {}
   for f in CDS_candidates:
      new_qual.update({f.qualifiers['source']:f.qualifiers['score']})
   if len(new_qual) == len(CDS_candidates): #definitely not chimera, proceed as usual
    new_feature = annotate_cds(CDS_candidates)
    if len(new_feature) < 85: #additional filtering of very short debries, they should still be traceable
      new_feature.type = 'orf'
    new_rec.features.append(new_feature)
   else: # chimera check here. Could be real chimera or two same kind hits... 
         # due to false positives (low score) or frameshift (same product)...
         # how to discriminate? Size of overlap and relative score?
         # if overlap larger than threhold or either scores much bigger than the other -> drop the false positive
         # if same product -> suggest frameshift and proceed as usual ##TO DO: frameshift corection???
         # otherwise try to divide the chimera.
         
       if verbose: print 'potential chimera detected in group containing',CDS_candidates[0].qualifiers['product']
       #first identify the two same-source features
       paired_list =[]
       unpaired_list =[]
       while len(CDS_candidates) > 0:
        f = CDS_candidates.pop()
        f_is_paired = False
#pairing algorithm is faulty. Should not rely on different orfs from the same group.
# if they end up here something must be gluing different orfs together
# consider removing the glue instead
#definitely exclude getorf results for chimera separation seeds. Use wise only?
        if f.qualifiers['source'] =='GeneWise':
         for g in CDS_candidates+paired_list+unpaired_list:
          if f.qualifiers['source'] == g.qualifiers['source']: 
           f_is_paired = True                                                                 
           break
        if f_is_paired and f not in paired_list:
          paired_list.append(f)
        if not f_is_paired and f not in unpaired_list:
          unpaired_list.append(f)
       
       if verbose: print 'number of paired and single features', len(paired_list), len(unpaired_list)
       paired_list = sorted(paired_list, key = lambda x: float(x.qualifiers['score']), reverse = True)
#       paired_list = sorted(paired_list, key = lambda x: len(x), reverse = True)
       if len(paired_list) == 0:
         CDS_candidates = unpaired_list
         new_feature = annotate_cds(CDS_candidates)
         if len(new_feature) < 85: #additional filtering of very short debries, they should still be traceable
          new_feature.type = 'orf'
         new_rec.features.append(new_feature)
         if verbose: print 'Not a a true chimera,', new_feature.qualifiers['product']
#       elif len(paired_list)> 2: 
#         print '\n Chimara separation failure possible, paired list:'
#         for xx in paired_list:
#          print xx.qualifiers['source'],xx.qualifiers['score']
#         for xx in unpaired_list:
#          print xx.qualifiers['source'],xx.qualifiers['score']
       else:         
        f = paired_list[0]
        g = paired_list[1]
       ### TO DO: what if there are more members than two? 
       #   The list should be sorted so that this does not matter. By score? By source? By length: longer first?
        f_score = float(f.qualifiers['score'])/len(f)
        g_score = float(g.qualifiers['score'])/len(g)
        if verbose: print f_score, g_score, big_score(f,g), '<-individual scores and overlap'
        if overlap(f,g) or f_score > 10* g_score or g_score > 10*f_score: #tuning needed
         if verbose:print f.qualifiers['product'], g.qualifiers['product'], 'are NOT part of chimeric orf, overlap or weak score'
         if f_score > g_score:
           if verbose:print 'will use', f_score, len(f), f.qualifiers['product'],f.location
           if verbose:print 'NOT     ', g_score, len(g), g.qualifiers['product'],g.location
           paired_list.remove(g)
         else:
           if verbose:print 'will use', g_score, len(g), g.qualifiers['product'],g.location
           if verbose:print 'NOT     ', f_score, len(f), f.qualifiers['product'],f.location
           paired_list.remove(f)
        else:
         if f.qualifiers['product']==g.qualifiers['product']:
           print 'multiple \033[91mFRAMESHIFT\033[0m errors suspected in',f.qualifiers['product']
         else:
           if verbose: print 'Chimeric ORF containing:', f.qualifiers['product'], g.qualifiers['product'], f.qualifiers['source']
           ## Try obvious thing first: cut the complete orf with inside trn.
           ## if there is one...
           for h in unpaired_list:
             if h.qualifiers['source'] == 'getorf': break
           ## longest is not always orf... so use getorf source... test it!
           if verbose: print 'chimeric orf is at:', h.location
           for t in list_tRNA:
               accepted = False
               if contain(t[0],h) and offset_contain(t[0],h,100) and offset_contain(t[0],h,-100) and len(t)>1:
               ## a bit more sophistication would be nice
                   for r in t:
                       if r.strand == h.strand:
                           accepted = True
                           break
                   if accepted: break
           if accepted:
            if verbose: print 'tRNA separating chimeric orf identified'
            fx=over3(r,h)
            gx=over5(r,h)
            if overlap(f,fx):
               CDS_candidates = [f,fx]
            if overlap(f,gx):
               CDS_candidates = [f,gx]
            new_feature = annotate_cds(CDS_candidates)
            new_feature.qualifiers['note']='resolved chimera'
            new_rec.features.append(new_feature)
            if overlap(g,fx):
               CDS_candidates = [g,fx]
            if overlap(g,gx):
               CDS_candidates = [g,gx]
            new_feature = annotate_cds(CDS_candidates)
            new_feature.qualifiers['note']='resolved chimera'
            new_rec.features.append(new_feature)
            paired_list =[]
            ##Control by list lengths is only temporary.
           else:
               if verbose: print 'no suitable tRNA found'
               ## The following does not have to work always, need more flexibility in finding the correct start codon...
               #  under3 coded... The result should be conservative, ea. there still could be overlap between orfs - may be resolved later
               paired_list = sorted([f,g], key = lambda x: x.location.start)
               strand = f.strand
               if strand == 1:
                   f = paired_list[0]
                   g = paired_list[1]
               else:
                   f = paired_list[1]
                   g = paired_list[0]
               if verbose: print 'attempting to separate', f.location, f.qualifiers['product'], g.location, g.qualifiers['product']
               gx = over5(f,h)
               fx,OK = under3(g,h)
               if OK:
                 for CDS_candidates in [[f,fx],[g,gx]]:
                   working_set = CDS_candidates
                   for x in unpaired_list:
                     if x.qualifiers['source']=='GeneWise' and overlap(x,working_set[0]):
                       working_set.append(x)
                   new_feature = annotate_cds(working_set)                    ##!!!blindly adding unpaired list twice is a risky workaround wise2 bug.
#                   new_feature = annotate_cds(CDS_candidates+unpaired_list) ## will adding only some members help?
                                                                            ##   may create problems in true vs. bogus chimeras
                   if verbose: print 'separated component:', new_feature.location, new_feature.qualifiers['product']
                   new_feature.qualifiers['note']='resolved chimera'
                   new_rec.features.append(new_feature)
                   paired_list = []
               else:
                 if verbose: print 'separation failed.'
                 #re-creating original list
                 paired_list.extend(unpaired_list)
                 unpaired_list = [] 
                 # test separation algoritm               
        if len(paired_list) > 0:
         CDS_candidates = paired_list + unpaired_list
         new_feature = annotate_cds(CDS_candidates)
         if len(unpaired_list) == 0:
           new_feature.qualifiers['product'] = f.qualifiers['product']+'+'+g.qualifiers['product']
         if verbose: print 'problematic feature annotated as', new_feature.location, new_feature.qualifiers['product']
         new_rec.features.append(new_feature)
                 
if verbose: print len(new_rec.features), ' features annotated (after CDS)'

#Now process tRNA, add congruent sets first:
# consider strict aa consensus here, use the best scoring infernal as prototype
# this is ensured by prior sorting by score and by source 
# so the first one ([0]) should be arwen, 
# the second ([1], but [0] after popping the first one) 
# should be best scoring infernal

#this is very arbitrary, may not give the best structure/location.
#consider sorting by score only and selecting the location with best score.

##it is still necessary to identify duplicates and delay their assignement
# therefore trn accounting should start even earlier.
# however perhaps separate list should be made for this purpouse by scanning the data twice...
del g
count_trn = {}
for trn_candidates in list_tRNA:
 if len(trn_candidates) > 0:
   f = trn_candidates[0]
   congruent = True
   for g in trn_candidates[1:]:
    if f.qualifiers['label'] <> g.qualifiers['label']:
      congruent = False
      break
   if congruent:
    try:
      count_trn[f.qualifiers['label']] += 1
    except KeyError:
      count_trn[f.qualifiers['label']] = 1

left_tRNA = []
while len(list_tRNA) > 0:
 trn_candidates = list_tRNA.pop(0)
 if len(trn_candidates) > 1:
    f = trn_candidates.pop(0)
    congruent = True
    for g in trn_candidates:
      if f.qualifiers['label'] <> g.qualifiers['label']:
        congruent = False
        break
    if congruent and count_trn[f.qualifiers['label']] == 1:
     g = trn_candidates[0]
     if verbose: print 'congruent annotation of', f.location, f.qualifiers['label']
     new_feature = annotate_trna(f)
     new_feature.qualifiers[g.qualifiers['source']] = g.qualifiers['score']
     new_rec.features.append(new_feature)
    else: #non-congruent or more than one candidate
     trn_candidates.append(f) #give it back
     left_tRNA.append(trn_candidates)
 else: # just one candidate/method
     left_tRNA.append(trn_candidates)
list_tRNA = left_tRNA

##Consider adding semi-congruence as a separate step:
# take into account clusters ==3, with 2 sources supporting the same specificity

# not really good... but does not seem to do much harm either... needs more testing
# at the very least should not create duplicates ...
# definitely with non-perfect-congruence duplicates should not be allowed!!!
# It does create problems: can discard the only existing possibility (fly and shit case)
# to fix it, all the trn accounting must be done earlier...

aa_in = []
for t in new_rec.features:
    if t.type == 'tRNA':
         aa_in.append(t.qualifiers['label'])
## to avoid aliasing copy element-by-element
aa_desired = []
for x in aa_list:
 aa_desired.append(x)

for x in aa_in:
    while x in aa_desired:
         aa_desired.remove(x)
# aa_desired contains a list of missing specificities. 
# Don't add anything not needed and update the list accordingly.

"""
left_tRNA = []
while len(list_tRNA) > 0:
 trn_candidates = list_tRNA.pop(0)
 if len(trn_candidates) == 3:
    trn_candidates = sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse = True)
    trn_candidates = sorted(trn_candidates, key = lambda x: x.qualifiers['source'])
    if verbose:
      for x in trn_candidates:
        print x.qualifiers['label'],x.qualifiers['source'],x.qualifiers['score'],
      print'.'
    f = trn_candidates.pop(0)
    congruent = False
    for g in trn_candidates:
      if f.qualifiers['label'] == g.qualifiers['label']:
        congruent = True
        if verbose: print 'match', f.qualifiers['label'], g.qualifiers['label']
        break
    if congruent and g.qualifiers['label'] in aa_desired:
#     g = trn_candidates[0]
     if verbose: print 'semi-congruent annotation of', f.location, f.qualifiers['label']
     new_feature = annotate_trna(f)
     new_feature.qualifiers[g.qualifiers['source']] = g.qualifiers['score']
     new_rec.features.append(new_feature)
     aa_desired.remove(g.qualifiers['label'])
    else: #non-congruent
     trn_candidates.append(f) #give it back
     left_tRNA.append(trn_candidates)
 else: # just one candidate ??? OR MORE THAN 3!!!
     left_tRNA.append(trn_candidates)
list_tRNA = left_tRNA
"""

if verbose: print len(list_tRNA), 'tRNA clusters left after consensus calling',len(aa_desired), 'still desired'

if verbose: print len(new_rec.features), 'features annotated after first rounds of tRNA'

## With single cm-scanning alternative annotation of rRNA genes is needed. 
## Without Amalgamation but with structure. Almost all work done in r2gff.py
"""
if len(list_12s) >0:
   list_12s = amalgamate(list_12s,tolerance=200)
   for loc in list_12s:
    new_feature = SeqFeature(loc, 'rRNA',{})
    new_feature.qualifiers['product'] = 'small rRNA subunit (12S)'
    new_feature.qualifiers['label'] = '12S'
    new_feature.qualifiers['source'] = 'rRNA'
    new_feature.qualifiers['score'] = 1
    new_rec.features.append(new_feature)
if len(list_16s) >0:
   list_16s = amalgamate(list_16s,tolerance=800)
   for loc in list_16s:      
    new_feature = SeqFeature(loc, 'rRNA',{})
    new_feature.qualifiers['product'] = 'large rRNA subunit (16S)'
    new_feature.qualifiers['label'] = '16S'
    new_feature.qualifiers['source'] = 'rRNA'
    new_feature.qualifiers['score'] = 1
    new_rec.features.append(new_feature)
"""
new_rec.features.extend(list_12s+list_16s) # all real work done in r2gff.py, trust it!
## but keep track of further changes, make sure modification of location is reflected in structure
## TO DO: HERE optional cutting with pA should be placed HERE. Currently pA within rRNA are ignored!
## It should be relatively simple - with shortening of the rRNA at 3' if pA is there... 
## Very similar to CDS cutting, but without frame checking
## Difference: structure must be adjusted. There is some potential for validation of the cut in structure.

if verbose: print len(new_rec.features), 'features annotated after rRNA calling'

print 'filtering...'

#make a sorted circular list, check for overlap
#if substantial overlap decide which one to keep.
# here the winner takes it all, no compromises
# implemented criteria:
#  -inclusion or nearly inclusion
#  -name and score
## repeat procedure until gene count remains the same - seems to work but check some more 
## TO DO: simplify this filter and revise factors - too many conflicting conditions to tune

while True:
 cdss =[]
 for f in new_rec.features:
  if f.type in ['CDS','rRNA']: #!! rRNA is in loop. Does it get modified? It should not!
   cdss.append(f)
 gene_count = len(cdss)

 cdsss = sorted(cdss, key= lambda x: x.location.start)
 circle_cds = cycle(cdsss)
 rec3 = SeqRecord(rec.seq)
 rec3.annotations = new_rec.annotations
 modifier = {}
 for gene in cdsss:
  gene.qualifiers['mark'] = False
  modifier.update({gene:len(gene)})
 gene1 = next(circle_cds)
 while True:
  gene1.qualifiers['mark'] = True
  gene2 = next(circle_cds)
  delta = min(len(gene1), len(gene2))
  delta = delta /2
  aln_score= big_score(gene1,gene2)
  lng = modifier[gene1]
  modifier.update({gene1:lng-aln_score})
  if verbose: print gene1.qualifiers['product'],modifier[gene1], gene2.qualifiers['product'],modifier[gene2]
  olap = (aln_score > delta)
  bigger2 = (len(gene2) > len(gene1))
  bigger1 = (len(gene1) > len(gene2))
  very_small1= len(gene1) < 85
  very_small2= len(gene2) < 85
  contained1 = contain(gene1,gene2)
  contained2 = contain(gene2,gene1)
  acceptable_score = 0.08 
##still not quite good... 
##may not be possible to accomodate all cases simultaneously.... 
##think of better diagnostics...
  bad_name1 = (not (gene1.qualifiers['product'] in mitochondria or gene1.type == 'rRNA')) or (gene1.qualifiers['score'] < acceptable_score)
  bad_name2 = (not (gene2.qualifiers['product'] in mitochondria or gene2.type == 'rRNA')) or (gene2.qualifiers['score'] < acceptable_score)
  muchsmaller1 = (len(gene2) > 3*len(gene1))
  muchsmaller2 = (len(gene1) > 3*len(gene2))
  bad_condition_for_first = (olap  and bad_name1) or contained1
  bad_condition_for_second = (olap  and bad_name2) or contained2
  if bad_condition_for_first and bad_condition_for_second:
    try:
     limit1 = length_limit[gene1.qualifiers['product']]
    except:
     limit1 = 85
    try:
     limit2 = length_limit[gene2.qualifiers['product']]
    except:
     limit2 = 85

    if muchsmaller1 or modifier[gene1]-big_score(gene1,gene2) < limit1: #drop 1
      if verbose: print 'g1 dropped, worse'
    elif muchsmaller2 or modifier[gene2]-big_score(gene1,gene2) < limit2: #drop 2 but retain 1
      rec3.features.append(gene1)
      gene2 = next(circle_cds)
      if verbose: print 'g2 dropped, worse'
    else: #retain both, can't decide
      rec3.features.append(gene1)
      if verbose: print 'OK, both problematic'
  elif bad_condition_for_first:
      if verbose: print 'g1 dropped, g2 OK'
  elif bad_condition_for_second:
      rec3.features.append(gene1)
      gene2=next(circle_cds)
      if verbose: print 'g2 dropped, g1 OK'
  else:
      rec3.features.append(gene1)
      if verbose: print 'OK, both good'
  gene1=gene2
  if gene2.qualifiers['mark']:
    break  

 for f in new_rec.features:
  if f.type not in ['CDS','rRNA']:
   rec3.features.append(f)
 new_rec = rec3  

 cdss =[]
 for f in new_rec.features:
  if f.type in ['CDS','rRNA']:
   cdss.append(f)
 new_gene_count = len(cdss)
 delta = gene_count - new_gene_count 
 if delta == 0: break #repeat until no change in gene count...

if verbose:print  new_gene_count,'non-tRNA genes retained'


## Dominance evaluation done after orf filtering and first round of trn annotations
## 
dominant = False
plus_l = []
minus_l = []
current_strand = None
for f in new_rec.features:
 if f.strand == 1: plus_l.append(f)
 if f.strand ==-1: minus_l.append(f)
pluss = len(plus_l)
minuss = len(minus_l) #number of genes counted
if pluss > minuss*5:
  dominant = 1
  current_strand = 1
if minuss > pluss*5:
  dominant = 1
  current_strand = -1 

global_strand = current_strand
global_dominant = dominant

if not global_dominant: #check the span of genes -> local strandness
 if verbose: print ' plus before:',len(plus_l)
 if len(plus_l) > 1: plus_l = amalgamate(plus_l,tolerance=650)
 if plus_l[-1].parts > 1:
  a = plus_l.pop()
  plus_l.extend(a.parts)
 if verbose: print ' plus after:',len(plus_l)
 if verbose: print ' minus before:',len(minus_l)
 if len(minus_l) > 1: minus_l = amalgamate(minus_l,tolerance=650)
 if minus_l[-1].parts >1:
  a = minus_l.pop()
  minus_l.extend(a.parts)
 if verbose: print ' minus after:',len(minus_l)
 if len(plus_l) > 1:
  plus_f = SeqFeature(CompoundLocation(plus_l),'misc',{})
 elif len(plus_l) > 0:
  plus_f = SeqFeature(plus_l[0],'misc',{})
 else: plus_f = None
 if len(minus_l) >1:
  minus_f = SeqFeature(CompoundLocation(minus_l),'misc',{})
 elif len(minus_l) > 0:
  minus_f = SeqFeature(minus_l[0],'misc',{})
 else: minus_f = None
#if verbose: print len(plus_f), len(minus_f), ' plus minus lengths'

# for the next round of tRNA, again consider the need and compatibility
aa_in = []
for t in new_rec.features:
    if t.type == 'tRNA':
         aminoacid = t.qualifiers['product'][-3:]
         A1 = t.qualifiers['label']
         aa = aminoacid + '-' + t.qualifiers['anticodon'][-4:-1]
         aa_in.append(A1)

## to avoid aliasing copy element-by-element
## WHY this is needed? Or is it?
aa_desired = []
for x in aa_list:
 aa_desired.append(x)

for x in aa_in:
    while x in aa_desired:
         aa_desired.remove(x)
#         print aa_desired, aa_in,aa_list

# First get rid of contained pseudo-tRNAs - by filtering the list.
# Need to check first representative only.

if verbose: print len(list_tRNA), ' tRNA clusters available before containment filtering'
tmp_list=[]
while len(list_tRNA)>0:
 t_group = list_tRNA.pop()
 f = t_group[0]
 good = True
 for h in new_rec.features:
    if h.type in ['CDS', 'rRNA'] and contain(f,h) and f.strand==h.strand:
        if h.type == 'CDS':
              hx = Oyala(f,h) # make sure strandness does not matter (use cds strand in over procedures)
# but strandness DOES matter! It is possible to have trn on complementary strand in compact genomes!
# Moreover, Oyala is not appropriate to check such cases as CDS will not be shortened/disrupted
# the only way to avoid is to include strandness in containment filtering and either
# - do no filtering for wrong strand
# - or make separate procedure for that, not using Oyala.
# However, this implies that all trn genes of a cluster have the same strandness, which may not be the case
# for this to work they should, so strandness would have to matter at clustering stage.
              certified_product = h.qualifiers['product']
              if len(certified_product) > 4:
                 actual_product = certified_product.split('+',1)
                 actual_limit = length_limit[actual_product[0]]+length_limit[actual_product[1]]
              else:
                 actual_limit = length_limit[certified_product]
              if (big_score(f,hx) > len(f)-4) or (len(hx) < actual_limit):
#                  delta = len(h) - len(hx)
                  good = False
        else: #rRNA case
              delta_step = min(len(f)/4,len(h)/20) ## still needs tuning
              for delta in range(0,len(h),delta_step):
               if not(offset_contain(f, h, delta) and offset_contain(f, h, -1*delta)):
                  break
              if len(h)-delta < length_limit[h.qualifiers['product']]:
                  good = False
                  if verbose: print 'rRNA would be shorten by', delta
        if not good:        
                  if verbose: print 'excessive overlap', f.qualifiers['label'],f.location,h.qualifiers['product'],h.location
                  break
 if good:
   tmp_list.append(t_group)
list_tRNA = tmp_list
if verbose: print len(list_tRNA), ' tRNA clusters available after containment filtering'
#End filtering tRNA for containment

if verbose: print 'Start singleton assignment' #label in one cluster only => certain assignment

while True:
 no_of_clusters = len(list_tRNA)
 #control the internal loop by checking out all aa_desired
 tmp_list = []
 while len(aa_desired) > 0:
  if verbose: print 'missing aa:', aa_desired
  accepted = False
  t_list = []
  for trn_c in list_tRNA:
   label_set=[]
   for t in trn_c:
    if t.qualifiers['label'] not in label_set:
     label_set.append(t.qualifiers['label']) 
   t_list.extend(label_set)
  if verbose: print 'current flatlist',t_list
  aa_d = aa_desired.pop(0)
  cases = t_list.count(aa_d)
  if aa_d == 'M2' and cases > 0:
    cases = 'delayed'
  if verbose: print cases, aa_d
  if cases == 0:
    if verbose: print 'no candidate tRNA for', aa_d
  elif cases == 1:

#this should be implemented differently - with dictionary access to the cluster, not with iteration.
#with such implementation considering all scenarios should be easier.

   for i in range(0,len(list_tRNA)):
     trn_candidates = list_tRNA[i]
     trn_candidates=sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse=True)
     for f in trn_candidates: 
#why is the list not sorted properly? Because of popping..Need to re-sort each time? Yes...
      aa = f.qualifiers['label']
      if verbose:print 'considering available', aa,'for desired:', aa_d
      if aa == aa_d:
       if verbose: print f.qualifiers['label'], f.qualifiers['source'], f.location, 'accepted'
       #the only option, no check needed
       new_trna = annotate_trna(f)
       new_rec.features.append(new_trna)
       accepted = True
       break
     if accepted:
       del list_tRNA[i]
       break
  else: #append the popped aa back to the list for cases >1
   tmp_list.append(aa_d)
   if verbose: print 'delayed assignment for',aa_d
 aa_desired = tmp_list
 if len(list_tRNA) == no_of_clusters:
   break

if verbose: print len(list_tRNA), ' tRNA clusters available after singleton assignment'

pAsites=[]
for f in rec.features:
 if f.qualifiers['source'] == 'bed':
   pAsites.append(f)


# singleton assignment should be repeated until no new genes are called...
# consider complicating this by permuting the aa_desired list and comparing the results...
# if the list is not too long this should be managable... 
# But would require significant complication of the code.

## At this point we should have only true multiple choices in aa_desired or M2 or empty aa_desired
## the assignment should be done with additional criteria. 
## Which are of course generally unknown/uncertain.
## Help to choose the real one with pA from bed? Could use that for forcing...
aa_left = []
while len(aa_desired) > 0:
 if verbose: print 'missing aa:', aa_desired
 accepted = False
 trn_candidates_flat = []
 for trn_c in list_tRNA:
   trn_candidates_flat.extend(trn_c)
 t_list = []
 for t in trn_candidates_flat:
      t_list.append(t.qualifiers['label'])
 if verbose: print 'current flatlist',t_list
 aa_d = aa_desired.pop(0)
 cases = t_list.count(aa_d)
 if verbose: print cases, aa_d
 if cases == 0:
   if verbose: print 'no candidate tRNA for', aa_d
 elif cases == 1:
   for i in range(0,len(list_tRNA)):
     trn_candidates = list_tRNA[i]
     trn_candidates=sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse=True)
     for f in trn_candidates:
      aa = f.qualifiers['label']
      if verbose:print 'considering available', aa,'for desired:', aa_d
      if aa == aa_d:
       #need to check condition, even for a singlet? Strandness at least?
        if global_dominant:
          dominant = global_dominant
          current_strand = global_strand
        else:
          dominant,current_strand = local_dominant(f,plus_f,minus_f)
        if ((dominant and (current_strand == f.strand)) or not dominant):# and float(f.qualifiers['score']) > 19:
          if verbose: print f.qualifiers['label'], f.qualifiers['source'], f.location, 'accepted'
          new_trna = annotate_trna(f)
          new_rec.features.append(new_trna)
          accepted = True
          break
     if accepted:
       del list_tRNA[i]
       break
 else: # multiple choice cases processed here
   # collect indexes to possibilities of aa_d
   candidate_indexes = []
   for i in range(0,len(list_tRNA)):
     trn_candidates = list_tRNA[i]
     trn_candidates=sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse=True)
     for f in trn_candidates:
      aa = f.qualifiers['label']
      if aa == aa_d and i not in candidate_indexes:
        candidate_indexes.append(i)
   # the list of indexes to lists with potential candidates.
   # store the index and some criteria, then iterate over the sorted flatlist
   flatlist_of_candidates = []
   for i in candidate_indexes:
    cluster = list_tRNA[i]
    for f in cluster:
     aa_OK = f.qualifiers['label'] == aa_d
     if global_dominant:
         dominant = global_dominant
         current_strand = global_strand
     else:
         dominant,current_strand = local_dominant(f,plus_f,minus_f)
     strandness_OK = ((dominant and (current_strand == f.strand)) or not dominant)
     if aa_OK:
      flatlist_of_candidates.append([f, i, strandness_OK])
                           #sorted(list_olf, key = lambda x: len(x.location), reverse = True)

   cluster_list = []
   for member in flatlist_of_candidates:
     cluster_list.append(member[1])
#   print cluster_list
   for member in flatlist_of_candidates:
     member.append(cluster_list.count(member[1]))
     pA_list = []
     for site in pAsites:
      if site.location.end in member[0] and site.strand == member[0].strand: #add strandness...
       pA_list.append(site)
     member.append(len(pA_list))
#count polyA overlaps and store them too. Then sort by the number of pA covered.


   flatlist_of_candidates = sorted(flatlist_of_candidates, key = lambda x: float(x[0].qualifiers['score']), reverse = True)
   flatlist_of_candidates = sorted(flatlist_of_candidates, key = lambda x: x[3], reverse = True) 
   flatlist_of_candidates = sorted(flatlist_of_candidates, key = lambda x: x[4], reverse = True) 
#Considering partial concordance... More hits from the same "i" are more important than score
#x[3] holds the number per cluster...
#x[4] holds the number of same strand polyA sites ending within the candidate 
#so the ones with pA are considered first
   if verbose:
     print 'flatlist:'
     for member in flatlist_of_candidates:
       print ' member:',member[0].qualifiers['label'], 'from cluster', member[1],', total# from the same cluster:', member[3],', pA sites:',member[4]
   for member in flatlist_of_candidates:
     f = member[0]
     i = member[1]
     if verbose: print 'considering', f.location, f.qualifiers['label'],'cluster', i
     strandness_OK = member[2]
     if strandness_OK:# and f.qualifiers['anticodon'][-4] in ['g','t']: #need tuning & debugging. 
       new_trn = annotate_trna(f)
       new_rec.features.append(new_trn)
       accepted = True
       if verbose: print f.qualifiers['label'], f.qualifiers['source'], f.location, 'accepted'
       del list_tRNA[i]
       accepted = True
       break
     else:
       if verbose: print 'not accepting', f.qualifiers['label'], f.location, strandness_OK,f.qualifiers['anticodon'][-4:-1]
 if not accepted:
     aa_left.append(aa_d)

# re-analyzing plus-minus here... May be important in the context of subsequent ncr filling

plus_l = []
minus_l = []
for f in new_rec.features:
 if f.strand == 1: plus_l.append(f)
 if f.strand ==-1: minus_l.append(f)
if not global_dominant: #check the span of genes -> local strandness
 if verbose: print ' plus before:',len(plus_l)
 plus_l = amalgamate(plus_l,tolerance=650)
 if plus_l[-1].parts > 1:
  a = plus_l.pop()
  plus_l.extend(a.parts)
 if verbose: print ' plus after:',len(plus_l)
 if verbose: print ' minus before:',len(minus_l)
 minus_l = amalgamate(minus_l,tolerance=650)
 if minus_l[-1].parts >1:
  a = minus_l.pop()
  minus_l.extend(a.parts)
 if verbose: print ' minus after:',len(minus_l)
 if len(plus_l) > 1:
  plus_f = SeqFeature(CompoundLocation(plus_l),'misc',{})
 elif len(plus_l) > 0:
  plus_f = SeqFeature(plus_l[0],'misc',{})
 else: plus_f = None
 if len(minus_l) >1:
  minus_f = SeqFeature(CompoundLocation(minus_l),'misc',{})
 elif len(minus_l) > 0:
  minus_f = SeqFeature(minus_l[0],'misc',{})
 else: minus_f = None


## The following is probably very rarely (if at all) used.
## NO, it is used and must be kept.

aa_desired = aa_left
while len(aa_desired) > 0:
 if verbose: print 'last attempt to assign missing aa:', aa_desired
 accepted = False
 trn_candidates_flat = []
 for trn_c in list_tRNA:
   trn_candidates_flat.extend(trn_c)
 t_list = []
 for t in trn_candidates_flat:
      t_list.append(t.qualifiers['label'])
 if verbose: print 'current flatlist',t_list
 aa_d = aa_desired.pop(0)
 cases = t_list.count(aa_d)
 if verbose: print cases, aa_d
 if cases == 0:
   if verbose: print 'no candidate tRNA for', aa_d
 elif cases == 1:
   for i in range(0,len(list_tRNA)):
     trn_candidates = list_tRNA[i]
     trn_candidates=sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse=True)
     for f in trn_candidates:
      aa = f.qualifiers['label']
      if verbose:print 'considering available', aa,'for desired:', aa_d
      if aa == aa_d:
       #need to check condition, even for a singlet?
        if global_dominant:
          dominant = global_dominant
          current_strand = global_strand
        else:
          dominant,current_strand = local_dominant(f,plus_f,minus_f)
        if ((dominant and (current_strand == f.strand)) or not dominant) and float(f.qualifiers['score']) > 0:
          if verbose: print f.qualifiers['label'], f.qualifiers['source'], f.location, 'accepted'
          new_trna = annotate_trna(f)
          new_rec.features.append(new_trna)
          accepted = True
          break
     if accepted:
       del list_tRNA[i]
       break
 else: # multiple choice cases processed here
   # collect indexes to possibilities of aa_d
   candidate_indexes = []
   for i in range(0,len(list_tRNA)):
     trn_candidates = list_tRNA[i]
     trn_candidates=sorted(trn_candidates, key = lambda x: float(x.qualifiers['score']), reverse=True)
     for f in trn_candidates:
      aa = f.qualifiers['label']
      if aa == aa_d and i not in candidate_indexes:
        candidate_indexes.append(i)
   # the list of indexes to lists with potential candidates.
   # store the index and some criteria, then iterate over the sorted flatlist
   flatlist_of_candidates = []
   for i in candidate_indexes:
    cluster = list_tRNA[i]
    for f in cluster:
     aa_OK = f.qualifiers['label'] == aa_d
     if global_dominant:
         dominant = global_dominant
         current_strand = global_strand
     else:
         dominant,current_strand = local_dominant(f,plus_f,minus_f)
     strandness_OK = ((dominant and (current_strand == f.strand)) or not dominant)
     if aa_OK:
      flatlist_of_candidates.append([f, i, strandness_OK])
   flatlist_of_candidates = sorted(flatlist_of_candidates, key = lambda x: float(x[0].qualifiers['score']), reverse = True)
                           #sorted(list_olf, key = lambda x: len(x.location), reverse = True)

   for member in flatlist_of_candidates:
     f = member[0]
     i = member[1]
     strandness_OK = member[2]
#     if (f.qualifiers['anticodon'][-4] in ['g','t']) or (f.qualifiers['anticodon'][-4:-1]=='cat'): 
#need tuning & debugging, consider removing this condition alltogether
     if True:
       new_trn = annotate_trna(f)
       new_rec.features.append(new_trn)
       accepted = True
       if verbose: print f.qualifiers['label'], f.qualifiers['source'], f.location, 'accepted'
       del list_tRNA[i]
       break
     else:
       if verbose: print 'not accepting', f.qualifiers['label'], f.location, strandness_OK,f.qualifiers['anticodon'][-4:-1]


trn_left_list =[]
for x in list_tRNA:
  for i in x: 
   if i not in trn_left_list:
    trn_left_list.append(i)

if verbose: print ' potential tRNA still left: ', len(trn_left_list)

#For Oyala procedure folow trna genes and check overlap with cds for each. 
#store each cds in rec3. Multiple modifications of the same cds by differet trns are not possible,
# the changes should be saved before going for next trn/rrn

print "Adjusting overlaps..."

OY=[] #Holds list of RNA genes overlapping with any CDS.
for f in new_rec.features:
 if f.type in ['tRNA','rRNA']:
   for g in new_rec.features:
    if g.type =='CDS' and overlap(f,g) and f.strand == g.strand:
     OY.append(f)
     break

for r in OY:
 rec3 = SeqRecord(rec.seq)
 rec3.annotations = new_rec.annotations
 for f in new_rec.features:
  fx = f
  if fx.type =='CDS' and overlap(r,fx) and r.strand == fx.strand:
   if verbose: print 'adjusting ', fx.qualifiers['product']
   if verbose: print fx.location,'before Oy'
   fx = Oyala(r,fx)
   if verbose: print fx.location,'after Oy'
  rec3.features.append(fx)
 new_rec = rec3
# something is not quite right here... Are r genes modified by Oyala? No... So print fx, not r location...


list_ncrs = collect_ncrs(new_rec)

master_copy_of_ncrs_sans_rep = list(list_ncrs)

# this will only be used for adding trn pseudogenes
# which could be located within repeats
 
## repeats - before NCRs!
## kmer-type dictionary approach retired
## use external vmatch (orders of magnitude faster)
## repeats are already annotated in the gff (v2gff.py)
word_size = 14 # fixed size k-mer tends to work better
#word_size = sequence_len / 1100

replist1 = []
for f in rec.features:
  if f.type =='repeat_region':
    replist1.append(f)

print len(replist1), 'repeats, filtering...'
replist2 = []
assoc ={'p':'inverted','supermax':'direct','tandem':'tandem'}
while len(replist1) >0:
 r1=replist1.pop(0)
 for r2 in replist1:
  if r1.qualifiers['ID'] == r2.qualifiers['ID']:
   if r1.location <> r2.location:
    de_virtualized_list = circularize([r1.location,r2.location])
#one small problem left: adjoinig parts. tbl2asn apparently insists on separating adjacent parts
#not sure how to deal with this. Seems more like a bug in tbl2asn than a solvable problem.
# added separating code to v2gff.py
    new_rep=SeqFeature(CompoundLocation(de_virtualized_list),'repeat_region',{})
    # sometimes repeats at the begining of the circle were getting zero-indexed location in gb output
    # Should be fixed now (compound location at the circle end-> circularize problem).
    new_rep.qualifiers['note'] = assoc[r1.qualifiers['repeat_type']]
    replist2.append(new_rep)

# filter duplicates and very short (arbitrary, smaller than two times word_size) repeats
replist1 = []
while len(replist2) >0:
 r1 = replist2.pop(0)
 mark = False
 for r2 in replist2:
  if r2.location == r1.location:
    r2.qualifiers['note'] = 'tandem' 
# is this always true that only tandem are duplicated? Or is this code really needed?
    mark = True
    break
 if not mark and len(r1) > word_size*2:
  replist1.append(r1)

print len(replist1),'repeats left'
new_rec.features.extend(replist1)
del replist1, replist2, assoc
#all repeats annotated

print "adding optional features in NCRs"
# Supplementary trn/CDS additions.
# possible implementation of some statistical measure of orf length expectancy?
# not bad, check and tune as needed. Extend the compatible scoring to glimmer(?)
  ### Rules implemented:
    # - tRNA pseudogenes significantly (?) overlapping ncr added first
    # - orfs prioritized based on length
    # - more overlap with ncr than with any other already present feature
    #   after Oyala with overlapping rna (?) features

# Need separate NCRs for ORFans and pseudo-trns?
# Trns can be added more aggressively, disregarding repeat regions.
# So for trn NCR list before repeat assignment should be stored (and not used for anything else).

#trns, 
# make sure the modification of list_ncrs inside the loop is safe.
#  There is a risk of skipping some elements and/or running out of index
#  conservatively safe way would be to re-start list trnsversal after each modification
# new implementation: (TO DO: check& debug)

list_ncrs = collect_ncrs(new_rec)
master_copy_of_ncrs = list(list_ncrs)

list_ncrs = list(master_copy_of_ncrs_sans_rep)

while True:
 # first filter out too short ncrs
 unprocessed = []
 while len(list_ncrs) > 0:
  ncr = list_ncrs.pop() #got the last element
  len_ncr = len(ncr)
  if len_ncr > 40:
   unprocessed.append(ncr)
   if deepverbose: print 'listing ncr', ncr.location
 list_ncrs = unprocessed

 #then compare coordinates of each with left trn candidates
 if verbose: print len(list_ncrs), '<-# long ncr locations'
 if len(list_ncrs) > 1:
  len_before = len(CompoundLocation(f2list(list_ncrs)))
 elif len(list_ncrs) > 0:
  len_before = len(list_ncrs[0])
 else: len_before = 0
 if verbose: print ' (start of trn loop) non-tivial NCRs length:', len_before
 unprocessed = []
 while len(list_ncrs) >0:
  ncr = list_ncrs.pop() #got the last element
  list_olf = []
  if verbose: print ncr.location
  for f in trn_left_list: #use flatlist? or keep using list_tRNA? Either way prepare a flatlist of overlaps
   if overlap(ncr,f):
    list_olf.append(f)
    if verbose: print '=>',f.location, len(list_olf)
  if len(list_olf) > 0: #OK, there are some candidates
   accepted = False
   if verbose: print len(list_olf), '<- #candidates'
   for f in list_olf:
     gain = big_score(ncr,f)
#     if gain > (len(f) - gain): #overlap bigger than the rest of it
     if gain > 0: #alternative, very relaxed condition
       if global_dominant:
           dominant = global_dominant
           current_strand = global_strand
       else: 
           dominant,current_strand = local_dominant(f,plus_f, minus_f)
       if not dominant or (dominant and (f.strand == current_strand)):
#         if verbose: print 'attempting addition of', f.qualifiers['label'],f.location
         new_pseudo = annotate_trna(f)
         new_pseudo.qualifiers['pseudo'] = None
         new_pseudo.qualifiers['product']=new_pseudo.qualifiers['product'][-3:]
         new_pseudo.type = 'misc'
         new_rec.features.append(new_pseudo) 
         accepted = True
         if verbose:print '->adding', new_pseudo.qualifiers['label'],new_pseudo.qualifiers['product'],new_pseudo.location
         new_ncr_loc = subtract(ncr,f)
         for loc in new_ncr_loc:
           new_ncr = SeqFeature(loc,'misc',{})
           unprocessed.append(new_ncr)
#         unprocessed.append(ncr) # to make orfans independent of pseudo - trn
         break # so that no overlapping features are added
   #Here either something was accepted
   # - if so leftovers were added to the unprocessed list for the second round
   #but if nothing was accepted we have to store the ncr back (for orfans)
   if not accepted:
     unprocessed.append(ncr)
  else: # no candidates, also store for orfans
   unprocessed.append(ncr)
 list_ncrs = unprocessed
 if len(list_ncrs) > 1:
  len_after = len(CompoundLocation(f2list(list_ncrs)))
 elif len(list_ncrs) >0:
  len_after = len(list_ncrs[0])
 else: len_after = 0
 if verbose: print ' (end of trn loop) left NCRs length:', len_after

 if len_after == len_before : break 
 # this will not work without elimination of ncrs... 
 # Otherwise at most one pseudo trn will be added per ncr
 # revert to the previous algo but store the copy of the original list and restore it here
 # watch out for aliases!!!

#ORFans 
#TO DO: needs major overhauling... Stats are not reliable enough.

list_ncrs = list(master_copy_of_ncrs)
while True:
 added = False
 for ncr in list_ncrs:
  if len(ncr) > 85: #consider only ncrs long enough to host true orf
   if verbose: print 'ncr at',ncr.location,len(ncr)
   list_olf = []
   for orf in orf_list:
    f = orf[0] #strange construct because list is not flat. Need to care about first (0) element only (because of prior sorting?)
    if overlap(ncr,f):
     list_olf.append(f)
   if len(list_olf) > 0:
    list_olf =  sorted(list_olf, key = lambda x: len(x.location), reverse = True)
    for f in list_olf:
        if verbose: print 'considering orf',f.location,len(f),f.qualifiers['score']
        gain = big_score(ncr,f)
        
        if (gain > (len(f) - gain)) or gain == len(ncr): #overlap bigger than the rest of orf
            if verbose: print 'significant gain ', gain, len(f)-gain,'extra'
            OK = 1
            fx = SeqFeature(f.location,f.type,{}) #make a copy to avoid aliasing
            fx.qualifiers=f.qualifiers
            for g in new_rec.features: ## iterating over ALL features seems excessive, 
                                        # consider limiting it to nearby features of certain kind
#              before_cut = len(f)
              if overlap(fx,g) and ('RNA' in g.type) and (f.strand == g.strand):
                fx = Oyala(g,fx)
                if big_score(fx,g) > 3: #failed Oyala
                  OK = 0
                  if verbose: print 'Oyala failed'
                  break
                else:
                  gain = big_score(fx,ncr) # Oyala OK, new feature, new length, new gain
              if big_score(fx,g) > gain/2: # this still needs tuning and testing...
                OK = 0
                if verbose: print 'too big overlap with neighbours'
                break
            #Add local strandness
            if not global_dominant:
             dominant,current_strand = local_dominant(f,plus_f,minus_f)
            else:
                dominant = global_dominant
                current_strand = global_strand
            if dominant and (f.strand <> current_strand):
               OK = 0
               if verbose: print 'wrong strand'
# adopt getorf scores here. Score > 20 indicates less than ~5% chance of obtaining such long orf by chance.
            if OK and (float(fx.qualifiers['score']) > 20  or ((gain > len(ncr)*0.7) and (len(fx)-gain < len(ncr)/2))): 
##tuning still needed, if these are just potential orfs, maybe relax?
             if verbose: print 'attempting addition of', f.qualifiers['product'],fx.location, 'gain:',gain,'ncr:',len(ncr),'left:',len(fx)-gain,'1%:', sequence_len*0.01
             new_location = fx.location
             new_feature = SeqFeature(new_location, 'CDS', {})
             new_feature.qualifiers=fx.qualifiers
             new_feature.qualifiers['product'] = 'orf'
             new_rec.features.append(new_feature)
             orf_list.remove(orf)
             added = True
             new_ncr_loc = subtract(ncr,fx)
             list_ncrs.remove(ncr)
             for loc in new_ncr_loc:
               new_ncr = SeqFeature(loc,'misc',{})
               list_ncrs.append(new_ncr)
             break
            else:
             if verbose: print 'too short',f.qualifiers['product'],fx.location, 'gain:',gain,'ncr:',len(ncr),'left:',len(fx)-gain,'1%:', sequence_len*0.01
  if added: break
 if not added: break   

cdss = []
for f in new_rec.features:
  if f.type in ['CDS']:
   f.qualifiers['mark']=False
   cdss.append(f)

# copy additional annotations of bed source (polyA sites). 
# Consider more fancy filtering.
# used in CDS refinement

if verbose and len(pAsites) >0: print "adjusting CDS to match pA"

max_score=0
for pA in pAsites:
 if pA.qualifiers['score'] > max_score:
   max_score=pA.qualifiers['score']
while len(pAsites)>0:
 i=pAsites.pop(0)
 if i.qualifiers['score']>=max_score/10:
  loc=FeatureLocation(i.location.start,i.location.start+1,i.strand)
  pA=SeqFeature(loc,'polyA site',{})
  pA.qualifiers['score']=str(i.qualifiers['score']/max(max_score,1)) #?consider something more fancy
  new_rec.features.append(pA)
  for f in cdss:
   if (pA.location.start in f) and (pA.strand == f.strand):
    if verbose: print 'polyA in CDS', f.qualifiers['product'],f.location,pA.location
    if len(f.location.parts) == 1:
      coordinates=f.location
      start=pA.location.start
    else:
      coordinates,start=pair_virt(f.location,pA.location)
      start=start.start
    distance_from_begin=start-coordinates.start
    distance_from_end=coordinates.end-start    
    if pA.strand == 1 and distance_from_begin < distance_from_end:
     cutter=SeqFeature(offset_loc(pA.location,(-1)),'misc',{})
     cutter.qualifiers['product']='pA'
     g = over5(cutter, f)
     f.location = g.location
     if verbose: print f.location
     
    elif pA.strand == -1 and distance_from_begin > distance_from_end: #untested! may not work as expected!
     cutter=SeqFeature(offset_loc(pA.location,1),'misc',{})
     cutter.qualifiers['product']='pA'
     g = over5(cutter, f)
     f.location = g.location
     if verbose: print f.location
     
    elif pA.strand == 1 and distance_from_begin > distance_from_end: 
     cutter = SeqFeature(pA.location,'misc',{})
     cutter.qualifiers['product']='pA'
     g = over3(cutter, f)
     f.location=g.location
     if verbose: print f.location
     
    elif pA.strand ==-1 and distance_from_begin < distance_from_end:
     cutter=SeqFeature(pA.location,'misc',{})
     cutter.qualifiers['product']='pA'
     g = over3(cutter, f)
     f.location=g.location
     if verbose: print f.location
   if deepverbose: print 'pA processed' 
try:
 del pAsites, max_score, i, g, start, distance_from_begin, distance_from_end, cutter, coordinates
except:
 pass
# using polyA information to adjust location.

# FINAL round of cds overlap detection, this time to reduce overlap length
# WITHOUT dropping anything, not touching non-CDS features
# should not be needed if pA info was reliable

if verbose: print 'cds overlap polishing procedure'

gene_count = len(cdss)
cdsss = sorted(cdss, key= lambda x: midpoint(x.location))
circle_cds = cycle(cdsss)
gene1 = next(circle_cds)
while True:
  if verbose: print 'next->'
  gene1.qualifiers['mark'] = True
  gene2 = next(circle_cds)
  bs=big_score(gene1,gene2)
  delta = 12
  if verbose:  print gene1.qualifiers['product'], gene2.qualifiers['product'],bs
  if bs > delta and gene1.strand == gene2.strand:
   if verbose: print ' attempting over5'
   if gene1.strand == 1:
    g2a = over5(gene1, gene2)
    gene2.location = g2a.location
   else:
    g1a = over5(gene2,gene1)
    gene1.location = g1a.location
   bs=big_score(gene1,gene2)
   if verbose: print ' past over5',gene1.qualifiers['product'], gene2.qualifiers['product'],bs
   if bs > delta and gene1.strand == gene2.strand:
    if verbose: print ' attempting over3'
    if gene1.strand == -1:
     g2a = over3(gene1, gene2)
     gene2.location = g2a.location
    else:
     g1a = over3(gene2,gene1)
     gene1.location = g1a.location
    bs=big_score(gene1,gene2)
    if verbose: print ' past over3', gene1.qualifiers['product'], gene2.qualifiers['product'],bs
    if bs > delta and gene1.strand == gene2.strand:
     #still unresolved... Try cutting in the middle of  the overlap and find valid combination of start codon (over5)
     #and matching virtual stop codon (over3). Use virtual cutting gene, with offset_loc
     if verbose: print ' attempting midpoint sweep'
     if gene1.strand == -1:
      cutter = SeqFeature(offset_loc(gene2.location,bs/2),'misc',{})
      cutter.qualifiers['product']='mid'
      g1a = over5(cutter,gene1)
      g2a = over3(g1a, gene2)
     else:
      cutter = SeqFeature(offset_loc(gene1.location,(-1)*bs/2),'misc',{})
      cutter.qualifiers['product']='mid'
      g2a = over5(cutter,gene2)
      g1a = over3(g2a,gene1)
     gene1.location = g1a.location
     gene2.location = g2a.location
     bs=big_score(gene1,gene2)
     if verbose: print ' past midpoint sweep', gene1.qualifiers['product'], gene2.qualifiers['product'],bs
  gene1=gene2 #move to the next gene
  if gene2.qualifiers['mark']: # exit the loop if the next one was processed
     break

mincdslen=1000
for f in new_rec.features:
   if f.type == 'CDS' and f.qualifiers['product'] in mitochondria:
     culen = len(f)
     if culen < mincdslen:
       mincdslen = culen

for f in new_rec.features:
   try:
    del f.qualifiers['mark']
    del f.qualifiers['evidence']
   except:
    pass
#cleanup gb file

#ORFs will be added later anyway.
#if short CDS are not desired, just dont create them.

   if f.type=='CDS' and len(f) < mincdslen: #move very short CDS features to optional type
                                    # shown along with prototype getorf features BUT not saved in tbl
     f.type='orf'
     if verbose: print 'diminishing CDS to orf',f.location

#Clear rRNA - tRNA overlaps. Structures need to be modified!!!
## By aposteriori comparison of f.location before and f.location after.
## should work fine, but the structure is not guaranteed to be correct 
##- may need sanitization. Which could be a problem becouse this works on RNA
## so would have to be done in correct strand (and for rRNA these are already flipped)

for f in new_rec.features:
  if f.type == 'rRNA':
    for g in new_rec.features:
      if f.strand == g.strand and overlap(f,g) and g.type == 'tRNA':
        before=copy.copy(f.location)
        if verbose: print f.qualifiers['product'],'\nloc:',before,'len:',len(f),'slen:', len(f.qualifiers['structure'])
        loclist = sorted(subtract(f,g), key = lambda x: len(x), reverse = True)
        if len(loclist) > 0:
         f.location = loclist[0]
        x1 = f.location.start-before.start
        x2 = f.location.end-before.start
        f.qualifiers['structure']=f.qualifiers['structure'][x1:x2]
        if verbose: print 'adjusted len:',len(f), 'slen:',len(f.qualifiers['structure'])
# syntax reminder; sorted(new_rec.features, key = lambda x: len(x), reverse = False):

onscreen(new_rec)

# orf identification and annotation here 
do_membrane_domains = True
phobius_rec =SeqRecord(new_rec.seq)
phobius_rec.annotations = new_rec.annotations
phobius_rec.description = new_rec.description
mark = False
if do_membrane_domains:
# mark is for tracking changes: turns true if any of orfs is named, triggering conditional onscreen
 for f in new_rec.features:
  if f.type == 'CDS':
    fx = SeqFeature(f.location,f.type,{})
    fx.qualifiers = f.qualifiers
    nt = f.extract(new_rec.seq)
    phas = len(nt) % 3
    if phas ==1:
      nt = nt+'AA'
    if phas ==2:
      nt = nt+'A'
    prot = nt.translate(table=ttable,to_stop=True,cds=True)
    prot_rec = SeqRecord(prot)
    prot_rec.annotations = new_rec.annotations
    prot_rec.id = str(f.qualifiers['product'])
    prot_rec.annotations.update({'topology':'linear'})
    prot_filename = name[0]+'_'+f.qualifiers['product']+'_at_'+str(f.location.start)
    prot_rec.desciption = 'translated'
    SeqIO.write(prot_rec,prot_filename+'.aa','fasta')

#    command = 'phobius.pl '+ prot_filename + ' 2> ' + redir + '| grep TRANSMEM > '+ prot_filename+'.pho'
#    os.system(command)

# calculate relative positions in the genome using feature location, virtualize if needed

    command = 'tmhmm '+ prot_filename + '.aa -noplot 2> ' + redir + '| grep TMhelix > '+ prot_filename+'.tmhmm'
    os.system(command)
# TO DO: figure out how to cleanup temporary directory left by tmhmm
#        what options control those? Apparently '-noplot' is enough...
    if len(f.location.parts) > 1:
      fx.location = virtualize(f.location)
    #simplified parser for filtered tmhmm output. Consider more sophistication (like moving all to gff3?)
    with open(prot_filename+'.tmhmm', 'r') as file_input:
     file_contents = file_input.read()
     lista_linii = file_contents.split("\n")
     for raw in lista_linii:
      col = raw.split()
      if len(col) > 3:
       start = int(col[-2])
       end = int(col[-1])
       f_type = 'domain'
       direct = fx.strand
       qual = {}
       ploc=FeatureLocation(start,end,strand=1)
### BUG ALERT PLEASE DEBUG THIS after switching to tmhmm from phobius
###    looks OK...domains do align at codon boundaries...
       if direct == 1:
        x_nt = start*3 + fx.location.start
        y_nt = end*3 + fx.location.start
       if direct == -1:
        x_nt = fx.location.end - end*3
        y_nt = fx.location.end - start*3
       loc = FeatureLocation(x_nt,y_nt,strand=direct)
       if loc.start > sequence_len:
         loc = loc + sequence_len*(-1)
       if loc.end > sequence_len: 
         loc = de_virtualize(loc)
       clean_feature = SeqFeature(loc, f_type, {})
       prot_feature=SeqFeature(ploc,f_type,{})
       clean_feature.qualifiers['source'] = 'tmhmm'
       clean_feature.qualifiers['note'] = 'transmembrane'
       prot_feature.qualifiers=clean_feature.qualifiers
       phobius_rec.features.append(clean_feature)
       prot_rec.features.append(prot_feature)
       try:
        os.remove(prot_filename + '.tmhmm')
       except:
        pass
    if str(fx.qualifiers['product']).upper() == 'ORF':
     if not (os.path.isfile(prot_filename+'.gff3')):
      command = 'python pfam2gff.py '+prot_filename+'.aa'
      os.system(command)
#      command = 'python super2gff.py '+work_aa_name+'.aa'
#      os.system(command)
#      command = 'python smart2gff.py '+work_aa_name+'.aa'
#      os.system(command)
#      command = 'python sf2gff.py '+work_aa_name+'.aa'
#      os.system(command)
#      command = 'python tigr2gff.py '+work_aa_name+'.aa'
#      os.system(command)
#      print ' aditional resources can be integrated through gff3 (run "db"2gff.py scripts with *.aa files)'
#    just store the protein records. Additional analyses can be done outside this script, if needed for selected files
#    most through hmmsearch with familiar domtblout format; this format can be parsed to gff with a separate script
#    in a very similar way as other annotation scripts. 
#    For performance reasons only pfam is left here.
     parse_gff3(prot_filename+'.gff3',prot_rec)
     for prot in prot_rec.features:
      try:
       if prot.qualifiers['name'] in mitochondria or float(prot.qualifiers['score']) > 20:
         print fx.qualifiers['product'],f.location,'identified as', prot.qualifiers['name']
         f.qualifiers['product'] = prot.qualifiers['name']
         mark = True
         break
      except:
       continue
#write all protein records to unique gb files. May consider adding these to submission for domain annotations.
    if len(prot_rec.features)>0:
     SeqIO.write(prot_rec,prot_filename+'.gb','genbank')
#     os.remove(prot_filename+'.aa')

# nice onscreen reporting, again but only if changed
print 'individual non-trivial protein records saved in gb files'
if mark:
 onscreen_cds(new_rec)


# add longest orfs...maybe better use their scores?
# TO DO: MAKE SURE orf locations were not modified (watch out for aliasing!)
# after that, turn the filtering back on.
longest_orfs =[]
for f in rec.features:
  if f.qualifiers['source'] == 'getorf':
   f.type='orf'
   longest_orfs.append(f)
#for f in new_rec.features:
# if f.type == 'CDS':
#  non_redundant_orfs = []
#  while len(longest_orfs)>0:
#   g=longest_orfs.pop()
#   if g.location <> f.location:
#    non_redundant_orfs.append(g)
#  longest_orfs = non_redundant_orfs

longest_orfs = sorted(longest_orfs, key = lambda x: int(x.qualifiers['score']), reverse = True)[:min(100,len(longest_orfs))]
new_rec.features.extend(longest_orfs)
if verbose: print '\n\n\norfs up to score', longest_orfs[-1].qualifiers['score']
#while len(longest_orfs) > 0:
#  i = longest_orfs.pop(0)
#  orf = SeqFeature(i.location, 'orf', {})
#  orf.qualifiers = i.qualifiers
#  if deepverbose: print i.location
#  new_rec.features.append(orf) 
#  if len(i) < sequence_len*0.013:break
#  if int(i.qualifiers['score']) < 10:break
del longest_orfs


#print 'checking for corrected', gbk_name
good_map = False

if (os.path.isfile(gbk_name)):
 print 'replacing data with edited',gbk_name
 new_rec_h = SeqIO.parse(gbk_name,"genbank")
 replacement = new_rec_h.next()
 try:
  taxid = replacement.annotations['organism']
 except:
  pass
 for f in replacement.features:
    for qual in f.qualifiers:
      f.qualifiers[qual] = ' '.join(f.qualifiers[qual])
      #gb parser is producing lists of strings whereas strings are expected here 
    if f.type in ['tRNA','rRNA']:
     try:
      f.qualifiers['structure'] = f.qualifiers['structure'].replace(' ','')
## more extensive structure sanitization needed (ea. in case of manual errors introduced)

     except:
      pass
 #check for shift and correct anticodon positions needed here
 if replacement.seq.upper() <> new_rec.seq.upper():
   if len(replacement.seq) <> sequence_len:
    print 'sequence length mismatch - this should not happen, correct initial fasta and start over'
    quit()
   middle = sequence_len/2
   half1 = new_rec.seq[:middle].upper()
   half2 = new_rec.seq[middle:].upper()
   pos1 = replacement.seq.upper().find(half1)
   pos2 = replacement.seq.upper().find(half2)
   if pos1 > 0:
     print pos1,'bp shift left detected, adding to embedded coordinates'
     ofset = pos1
   elif pos2 > 0:
     print middle-pos2, 'bp shift right detected, subtracting from embedded coordinates'
     ofset = pos2-middle
   else:
     print 'gross sequence mismatch - this should not happen, correct initial fasta and start over'
     quit()
   for f in replacement.features:
### Only "anticodon" qualifiers need this correction but the lists of types and qualifiers can be expanded.
### In particular "transl_except" in 'CDS' type could also be amended
### but currenlty it is created on-the-fly during tbl generation only and not stored in any gb

    if f.type in ['tRNA','misc']:
     try:
      for amendable in ['anticodon']:
       parts = f.qualifiers[amendable][1:-1].split(',',1)
       orig_loc = gb2loc(parts[0])
       if verbose: print '\n<-', f.qualifiers[amendable]
       new_loc = orig_loc + ofset
       if new_loc.start > sequence_len:
          new_loc = new_loc + sequence_len*(-1)
       if new_loc.start < 0:
          new_loc = new_loc + sequence_len
       if new_loc.end > sequence_len:
          new_loc = de_virtualize(new_loc)
       f.qualifiers[amendable] = '(pos:'+ loc2gb(new_loc)+','+parts[1]+')'
       if verbose: print '->', f.qualifiers[amendable]
     except:
       continue
   # regenerate vmatch indexes here - then dotplot can be moved to show replacement
   # overwrite original fasta? Seems to be the best choice but will cause problems...
   # orelse change the name... Will cause other problems... potentially...
   fasta_name = name[0]+'_mod.fa'
   SeqIO.write(replacement, fasta_name, "fasta")
   command ='python v2gff.py '+ fasta_name
   os.system(command)
 elif verbose:
  print 'same sequence, no shift in coordinates'
 
 new_rec = replacement
# SeqIO.write(new_rec,fasta_name+'.gb','genbank')
 good_map = True

 onscreen(new_rec)

# slice-in trna secondary structures.
zero_string = '.'*sequence_len
new_rec.letter_annotations['secondary_structure'] = zero_string

# since slicing can mess up annotations, prevent the list
# sorting provides nice opportunity to display them in order, as a side effect
# all annotations should be in new_rec at this point. Consider moving any downstream additions before this.
# ONE PROBLEM: if the feature was truncated the structure will not fit. 
## Make sure the truncation also truncates structure!! 
## Not so important for tRNA - they can be kept uncut (check!). But rRNA may be shortened!
## There is no way to tell which side was cut so keeping the track of these changes beforehead is the only way.

### BUG ALERT: multipart structural features!!! Slicing those in may change the length of the sequence!!!
### fixed but further debugging advised.

annotations = sorted(new_rec.features, key=lambda x: midpoint(x.location))
new_rec.features =[] 
print 'RNA structures'
for f in annotations:
    if 'RNA' in f.type:
     label = f.qualifiers['label']+'   '
     strand ='+'
     if f.location.strand ==-1:
      strand ='-'
     print strand,label[:3],
     if True: 
##??assume that structures stored in qualifiers are already correctly oriented. Make sure it is so.
      ss = f.qualifiers['structure']
      if verbose:
        print len(ss), len(f)
        
      dna_rec = slice_rec(new_rec,f)
      dna_rec.letter_annotations['secondary_structure'] = ss
      if len(f.location.parts) >1: 
         start = f.location.parts[0].start
         end = f.location.parts[-1].end
         if f.strand ==-1:
          start = f.location.parts[-1].start
          end = f.location.parts[0].end
         border = sequence_len-start
         new_rec = dna_rec[border:]+ new_rec[end:start]+dna_rec[:border]
      else:
         new_rec = new_rec[:f.location.start]+dna_rec+new_rec[f.location.end:]
#      if len(ss) < 100: may not work for fragments
      if f.type=='tRNA':
       ascore = '0'
       iscore = '0'
       try:
        ascore = f.qualifiers['arwen']
       except:
        pass
       try:
        iscore = f.qualifiers['infernal']
       except:
        pass
       score = max(float(ascore),float(iscore))*(-1)
       if f.strand == -1:
         ss=flip_struct(ss)
       ss=ss.replace('___',f.qualifiers['anticodon'][-4:-1].upper())
       print ss, "{:.2f}".format(score)#,shannon(f.extract(new_rec.seq).upper())
       if verbose:print strand,label[:3],dna_rec.seq
      else:
       score=f.qualifiers['score']
#       print 'rRNA'
       # simplified structure for printing. Reduce duplicates?
       i=0
       sscp=copy.copy(ss)
       while i < len(sscp)-1:
        if (not (sscp[i] in ['(',')','<','>','{','}','[',']'])) or sscp[i] == sscp[i+1]:
         sscp=sscp[:i]+sscp[i+1:]
        else:
         i+=1
       print sscp,score#, shannon(f.extract(new_rec.seq).upper())

new_rec.features = annotations

print 'ultrafast dotplot'

# prepare sorted list of features to be color coded on dotplot
# store tuples of: (start coordinate, length, color of shade, color on map, label for map)
# remember to calculate the first two for minus features correctly

bkf = sorted(new_rec.features, key=lambda x: x.location.start)
ckf =[]
cuco = cycle(['khaki','beige'])
for f in bkf:
 if f.type in ['rRNA','CDS','tRNA']:
    if len(f) > sequence_len /50: ## Tune dependence on length: similar to map width?
     if 'RNA' in f.type:
      lbl = f.qualifiers['label'][0:2]
     else:
      lbl = f.qualifiers['product'][0]+f.qualifiers['product'][-1]
    else:
      lbl = ''
    mapc = 'white'
    if 'RNA' in f.type:
       mapc = 'gray'
    x1 = int(f.location.parts[0].start)
#    delta = len(f) ## this will not work correctly on compound locations
    delta = len(f.location.parts[0]) # this will show only the first segment of each feature
                                     # consider iterating over all location parts
    if f.location.strand == -1:
       x1 = x1 + delta
       delta =-1*delta
    ckf.append((x1, delta, next(cuco),mapc,lbl))
del bkf, cuco, mapc, lbl, x1, delta   

   
kmer =14
results_f = fasta_name+'.'+ str(kmer)
if not (os.path.isfile(results_f)):
   command = 'vmatch -e 1 -p -d -allmax -l ' + str(kmer) + ' '+ fasta_name + ' > '+ results_f
   if verbose:print 'running', command
   os.system(command)
else:
   if verbose:print 're-using',results_f
xp =[]
yp = []
xr =[]
yr=[]
with open(results_f, 'r') as v_input:
   lista_linii = v_input.readlines()
   for linia in lista_linii:
       if linia[0] <>'#':    
	kol = linia.split()
        start1 = int(kol[2])
        size1 =  int(kol[0])
        end1 = start1+size1
        start2 = int(kol[6])
        size2 = int(kol[4])
        end2 = start2+size2
        x = []
        y = []
        step = max(size1/kmer,size2/kmer)
        for i in range(0,step):
         x.append(start1+i*kmer)
         y.append(start2+i*kmer)
        if kol[3] == 'D':
           yp.extend(y)
           xp.extend(x)
        if kol[3] == 'P':
           yr.extend(y)
           xr.extend(x)
           
#sf = sequence_len/2000
sf = 10

#mapw = 400
mapw = sequence_len / 50
origin = -0.75*mapw
pyplot.figure('dotplot',figsize=[sf,sf],dpi=100)

pyplot.gca().set_aspect('equal',adjustable='box')

pyplot.scatter(xp+yp,yp+xp, marker = '.', s=4, facecolor = 'black',alpha = 0.5, linewidth = 0)
pyplot.scatter(xr+yr,yr+xr, marker = '.', s=4, facecolor = 'red', alpha = 0.5, linewidth = 0)

pyplot.plot([0,sequence_len],[0,sequence_len],linewidth=0.5, zorder = 0.1)
for dat in ckf:
  x1 = dat[0]
  deltax = dat[1]
  cuco = dat[2]
  mapc = dat[3]
  lw = 0.01
  lbl = dat[4]
  pyplot.gca().add_patch(matplotlib.patches.Rectangle((x1,0),deltax, sequence_len, fc =cuco, alpha = 0.2, linewidth =0, zorder = 0))
  pyplot.gca().add_patch(matplotlib.patches.FancyArrow(x1,origin,deltax, 0, width = mapw, head_width = mapw, head_length = min(mapw/4,abs(deltax)-1), length_includes_head = True, ec = 'black',fc = mapc, alpha = 0.5, linewidth = lw, zorder = 10))
  
  pyplot.text(x1+deltax/2-mapw/20*copysign(1,deltax),origin-mapw/4,lbl, fontsize = 6,ha='center',zorder = 11)
  pyplot.text(origin-mapw/2.5,x1+deltax/2-mapw/8*copysign(1,deltax),lbl, fontsize = 6,va='center',zorder = 11)

  pyplot.gca().add_patch(matplotlib.patches.Rectangle((0,x1),sequence_len, deltax, fc =cuco, alpha = 0.2, linewidth =0, zorder = 0))
  pyplot.gca().add_patch(matplotlib.patches.FancyArrow(origin,x1,0,deltax, width = mapw, head_width = mapw, head_length = min(mapw/4,abs(deltax)-1), length_includes_head = True, ec = 'black',fc = mapc, alpha = 0.5, linewidth = lw, zorder = 10))
  
pyplot.xlabel(taxid)
pyplot.show()
pyplot.savefig(name[0]+'_dotplot.pdf')
pyplot.close('dotplot')

print 'dotplot written to:\033[96m',name[0]+'_dotplot.pdf\033[0m'

if good_map:
#  do_membrane_domains = True
  domain_features = []
  other_features = []
  while len(new_rec.features) > 0:
   f = new_rec.features.pop(0)
   if f.type == 'domain':
    domain_features.append(f)
   else:
    other_features.append(f)
  phobius_rec.features = domain_features
#  phobius_rec.features = []
  new_rec.features = other_features
# not sure if this separation of features is needed

#print new_rec.letter_annotations['secondary_structure']  
print '\npreparing GB submission\033[96m', tbl_name,'\033[0m'

allowed_tags=['product', 'anticodon']
allowed_types=['CDS','rRNA','tRNA','repeat_region']
with open(tbl_name, 'w') as tbl:
 tbl.write('>Feature ' + name[0]+'\n')
 for f in new_rec.features:
  if f.type in allowed_types:
    exception = None
    fx = f.location
    phase = len(fx) % 3
    if f.strand == -1:
     if f.type == 'CDS' and phase <> 0:
        if verbose: print f.qualifiers['product'], phase*f.strand
        if phase == 1:
          last_codon = FeatureLocation(fx.parts[-1].start,fx.parts[-1].start+1,fx.strand)  
        elif phase == 2:
          fx = diminute_start(f.location)                                                 ## Workaround tbl2asn bug: shorten the CDS by one
          last_codon = FeatureLocation(fx.parts[-1].start,fx.parts[-1].start+1,fx.strand) ## should be start+2, if not for the bug
#          last_codon = FeatureLocation(f.location.parts[-1].start,f.location.parts[-1].start+2,f.strand)
        if last_codon.end > sequence_len:
          last_codon = de_virtualize(last_codon)
        exception = '\t\t\ttransl_except\t(pos:' + loc2gb(last_codon) + ',aa:TERM)\n'
     tbl.write(str(fx.parts[0].end)+'\t'+str(fx.parts[0].start+1)+'\t'+str(f.type)+'\n')
     if len(f.location.parts) > 1:
      for i in range(1,len(fx.parts)):
            tbl.write(str(fx.parts[i].end)+'\t'+str(fx.parts[i].start+1)+'\n')
    else:
     if f.type == 'CDS' and phase <> 0:
        if verbose: print f.qualifiers['product'], phase*f.strand
        if phase == 1:
          last_codon = FeatureLocation(fx.parts[-1].end-1,fx.parts[-1].end, fx.strand)   
        elif phase == 2:
          fx = diminute_end(f.location)                                                        ## Workaround tbl2asn bug: shorten the CDS by one
          last_codon = FeatureLocation(fx.parts[-1].end-1,fx.parts[-1].end, fx.strand)         ## should be end-1, if not for the bug
#          last_codon = FeatureLocation(f.location.parts[-1].end-2,f.location.parts[-1].end, f.strand)  
        if last_codon.start < 0:
          last_codon = de_virtualize(last_codon + sequence_len)
        exception = '\t\t\ttransl_except\t(pos:' + loc2gb(last_codon) + ',aa:TERM)\n'
     tbl.write(str(fx.parts[0].start+1)+'\t'+str(fx.parts[0].end)+'\t'+str(f.type)+'\n')
     if len(fx.parts) > 1:
      for i in range(1,len(fx.parts)):
           tbl.write(str(fx.parts[i].start+1)+'\t'+str(fx.parts[i].end)+'\n')
    for key,val in f.qualifiers.iteritems():
      if key in allowed_tags:
       tbl.write('\t\t\t'+str(key)+'\t'+str(val)+'\n')
    if f.type == 'repeat_region':
      exception = '\t\t\tnote\t'+f.qualifiers['note']+'\n' 
    if exception:
      tbl.write(exception)

# running tbl2asn

tmp_rec = SeqRecord(new_rec.seq)
tmp_rec.name = name[0]
tmp_rec.id = name[0]
tmp_rec.description = '[organism='+taxid+'] [location=mitochondrion] [mgcode='+ttable_str+'] [topology=circular]'

SeqIO.write(tmp_rec, fsa_name, 'fasta')
try:
 os.remove(val_name,sqn_name,gbf_name)
except:
 pass 
command ='tbl2asn -as -i '+ fsa_name + ' -Vvb 2> '+redir
os.system(command)
if os.path.isfile(val_name):
 with open(val_name, 'r') as fh:
  errors = fh.read().split('\n')
  for L in errors:
   if verbose or (len(L)> 0 and 'Circular' not in L):
    print L
  print 'tbl2asn output written to files:\033[96m', sqn_name, gbf_name, val_name,'\033[0m'
  if verbose:  print 'for real submission use the tbl with appropriate fasta, sbt and -T option'
else:
    print 'tbl2asn failed to generate', val_name

if verbose: print 'preparing letter annotations'

have_coverage = False
cov_file=name[0]+'.coverage'
if os.path.isfile(cov_file):
 value_list=parse_csv(cov_file)
 if len(value_list)==sequence_len:
  new_rec.letter_annotations['coverage']=parse_csv(cov_file)
  have_coverage = True
  print '\nwill use coverage data in graphics\n'
 else:
  print 'length mismatch',len(value_list),sequence_len
# needed for stats filtered by codon position (1,2,3; 0 for non-CDS)
prepare_codon_stats = True
if prepare_codon_stats:
  # add the mask of codon positions
 in_codon_positions = [0]*sequence_len 
 for f in sorted(new_rec.features, key = lambda x: len(x), reverse = False):
## Sorting is a quick hack to "solve" overlaps.
    if f.type =='CDS':
            codon_pos = 0
            for pos in f:
              codon_pos += 1
              in_codon_positions[pos] = codon_pos
              if codon_pos == 3: codon_pos = 0
 new_rec.letter_annotations['codon'] = in_codon_positions

# degenerate:  4, 3, 2 or 1 aa possible, 0 for non-CDS 
# meaning: 0-fold (all 4 different), 2-fold(? 2 or 3 possibilities), 4-fold (just one)
# neutral are: definitely all 1, probably many 0 and maybe some 2 and 3 but not 4
# non-neutral are: definitely all 4, probably some 0 as well as some 3 and 2, but not 1
# CAUTION: too narrow filter will cause math problems

#no_graphics = False
#if no_graphics: quit()

prepare_degenerate_stats = True
if prepare_degenerate_stats:
  zero1 = [0]*sequence_len
#  for pos in range(0, sequence_len):
#    zero1.append(0)
  new_rec.letter_annotations['degenerate'] = list(zero1)
  new_rec.letter_annotations['prop']= list(zero1)
  if True:
   for i in range(0,sequence_len):
    if new_rec.letter_annotations['secondary_structure'][i] in ['<','>','[',']','(',')','{','}','_']:
#     print 'hit'
     new_rec.letter_annotations['prop'][i] = 4
#  print new_rec.letter_annotations['prop']
  new_rec.letter_annotations['gc_mild'] = list(zero1)
  new_rec.letter_annotations['at_mild'] = list(zero1)
  dublet = new_rec + new_rec
  for f in new_rec.features:
    sequ = dublet.seq                 
    # virtualization + rec duplication (for 'codon') + recalculation (=in-line de-virtualization)
    if f.type =='CDS':
       fx = SeqFeature(f.location,f.type,{})
       fx.qualifiers = f.qualifiers
#       print fx.location, f.location
       if len(fx.location.parts) > 1:
         fx.location = virtualize(fx.location)
#       print 'in feature:', fx.location, fx.qualifiers['product']
#       for pos in f:
#         print pos, dublet.letter_annotations['codon'][pos], ' ',
       for pos in fx:
         a_dict = {}
         g_dict = {}
         in_codon_number = dublet.letter_annotations['codon'][pos]
#         print pos,in_codon_number, '#in codon'
         if f.strand == 1:
           if in_codon_number == 1:
             triplet = sequ[pos:pos+3]
           if in_codon_number == 2:
             triplet = sequ[pos-1:pos+2]
           if in_codon_number == 3:
             triplet = sequ[pos-2:pos+1]
         if f.strand == -1:
           if in_codon_number == 1:
             triplet = sequ[pos-2:pos+1].reverse_complement()
           if in_codon_number == 2:
             triplet = sequ[pos-1:pos+2].reverse_complement()
           if in_codon_number == 3:
             triplet = sequ[pos:pos+3].reverse_complement()
#         print triplet, 'triplet'
         current_codon = Seq(str(triplet.upper()),IUPAC.ambiguous_dna) 
#         print current_codon
         for nt in ['A','C','T','G']:
           if in_codon_number == 1:
             mutated_codon = nt + current_codon[1:3]
           if in_codon_number == 2:
#             print 'cc', current_codon
             mutated_codon = current_codon[0]+ nt + current_codon[2]
           if in_codon_number == 3:
             mutated_codon = current_codon[0:2] + nt
           mutated_aa = Seq(str(mutated_codon)).translate(table=ttable)
           a_dict.update({str(mutated_aa):str(mutated_codon)})
           mutated_group = aa_prop[str(mutated_aa)]
           g_dict.update({mutated_group:mutated_aa})
         g_mld={}
         for nt in ['C','G']:
           if in_codon_number == 1:
             mutated_codon = nt + current_codon[1:3]
           if in_codon_number == 2:
#             print 'cc', current_codon
             mutated_codon = current_codon[0]+ nt + current_codon[2]
           if in_codon_number == 3:
             mutated_codon = current_codon[0:2] + nt
           mutated_aa = Seq(str(mutated_codon)).translate(table=ttable)
           mutated_group = aa_prop[str(mutated_aa)]
           g_mld.update({mutated_group:mutated_aa})
         a_mld={}
         for nt in ['A','T']:
           if in_codon_number == 1:
             mutated_codon = nt + current_codon[1:3]
           if in_codon_number == 2:
#             print 'cc', current_codon
             mutated_codon = current_codon[0]+ nt + current_codon[2]
           if in_codon_number == 3:
             mutated_codon = current_codon[0:2] + nt
           mutated_aa = Seq(str(mutated_codon)).translate(table=ttable)
           mutated_group = aa_prop[str(mutated_aa)]
           a_mld.update({mutated_group:mutated_aa})
         real_pos = pos
         if real_pos >= sequence_len: #why ==? Are the letter annotations really synchronized?
           real_pos = real_pos - sequence_len
#         print pos, real_pos, str(current_codon), in_codon_number
         new_rec.letter_annotations['degenerate'][real_pos] = len(a_dict)
         ## group assessment here; counting how many groups the aa's belong to (g_dict, prop).
         new_rec.letter_annotations['prop'][real_pos] = len(g_dict)
         ## additional letter-specific dicts, considering only specific changes. Not much use, probably not needed.
         ## could make another dict based on this second-level translation (codon->aa->group)
         new_rec.letter_annotations['gc_mild'][real_pos] = len(g_mld)
         new_rec.letter_annotations['at_mild'][real_pos] = len(a_mld)
#  print '\n',new_rec.letter_annotations['prop']
#for global phase tracking, currently not used
prepare_phase_stats = False
if prepare_phase_stats:
  zero2 = []
  for pos in range(0, sequence_len):
    zero2.append(0)
  new_rec.letter_annotations['phase'] = zero2
  for pos in range(0, sequence_len):
    new_rec.letter_annotations['phase'][pos] = pos % 3

if verbose: print 'letter annotations created'

##  test correctness of letter_annotations here:
#example_feature = new_rec.features[4]
#for example_feature in new_rec.features:
# if example_feature.type == 'CDS':
#  example_start = example_feature.location.parts[0].start
#  example_end = example_feature.location.parts[-1].end
#  if example_start > example_end:
#   sub_rec = new_rec[example_start:] + new_rec[:example_end]
#  else:
#   sub_rec = new_rec[example_start:example_end]
#  print example_feature.qualifiers['product'],sub_rec.letter_annotations['degenerate'],sub_rec.letter_annotations['codon']


#e_Nc & gene composition stuff here
# only GCn is calculated, but this could relatively easy be extended to any other stat
# all stats are stored per gene in csv.
# TO DO: implement per-strand stats but this should not be high priority.

if verbose: print "codon usage:"
ctable_name = name[0]+'_'+ttable_str+'_ECN.csv'
fh = open(ctable_name,'w')
fh.write(name[0]+',Nc,Nc_TM,NC_nonTM,GCn,GCn_TM,GCn_nonTM,len,len_TM,len_GCn,len_GCn_TM,pos\n')
print 'effective number of codons in CDS, TM and non-TM plus %GC at neutral sites:'
rox={}
genome=[]
ox=[]
for f in new_rec.features:
 if f.type =='CDS':
   frecord = slice_rec(new_rec,f)
   full_rec=FeatureLocation(0,len(frecord))
   bogusCDS=copy.copy(f)
   bogusCDS.location=full_rec
   bogusCDS.strand=f.strand
   frecord.features =[bogusCDS]
   cux = codon_usage(frecord)
   vecn = ECNa(cux)
   vgcn, filtered_len = r_con(frecord, [1], ['G','C'], 'prop')
   mrecord = slice_rec(phobius_rec,f)
   try:
    del tm_record
   except:
    pass
   for g in mrecord.features:
    g_record = slice_rec(frecord,g)
    try:
     tm_record += g_record
    except:
     tm_record = g_record
   try: #need to try because not all CDS have transmembrane domains
    full_rec=FeatureLocation(0,len(tm_record))
    bogusTM=copy.deepcopy(bogusCDS)
    bogusTM.location=full_rec
    bogusTM.strand=f.strand
    tm_record.features=[bogusTM]
    tm_cux=codon_usage(tm_record)
    mgcn, mfiltered_len = r_con(tm_record, [1], ['G','C'], 'prop')
    tm_len=len(tm_record)
   except:
    tm_cux={}
    mgcn=0
    mfiltered_len = 1 #need to do better? Store without division and set to zero?
    tm_len=0
   mecn = ECNa(tm_cux)
#   nontm_cux=dsub(copy.deepcopy(cux),tm_cux) #dsub modifies argument (unlike dsum)
   nontm_cux=dsub(cux,tm_cux)
   necn = ECNa(nontm_cux)
   wrl=[f.qualifiers['product'], str(vecn), str(mecn), str(necn), str(vgcn/float(filtered_len)), str(mgcn/float(mfiltered_len)), str((vgcn-mgcn)/float(filtered_len-mfiltered_len)), str(len(f)), str(tm_len) ,str(filtered_len), str(mfiltered_len), str(midpoint(f.location))]
   fh.write(",".join(wrl)+'\n')
   frecord.annotations['utable']=copy.deepcopy(cux)
   frecord.annotations['GCn']=copy.deepcopy((vgcn,filtered_len))
   frecord.annotations['ttable']=copy.deepcopy(tm_cux)
   frecord.annotations['GCnTM']=copy.deepcopy((mgcn,mfiltered_len))
   genome.append(frecord)
   for C in OXPHOS:
    if f.qualifiers['product'] in OXPHOS[C]:
     ox.append(frecord)
     try:
      rox[C].append(frecord)
     except KeyError:
      rox[C] = [frecord]
description={}
for C in rox:
 desc = 'Complex '+C
 description[C] = desc
cases =sorted(rox.keys())
rox['OXPHOS']=ox
description['OXPHOS']= 'OXPHOS'
rox['genome']=genome
description['genome']='genome'
cases.append('OXPHOS')
cases.append('genome')
for C in cases:
 total=get_stats(rox[C])
#could add some faisafe condition (?on total len) here
 total_cux=dsum(total['utable'])
 total_tm=dsum(total['ttable'])
 total_n=dsub(total_cux,total_tm)
 vecn=ECNa(total_cux)
 mecn=ECNa(total_tm)
 necn=ECNa(total_n)
 vgcn=sum(i for i, _ in total['GCn'])
 filtered_len=sum(i for _,i in total['GCn'])
 mgcn = sum(i for i, _ in total['GCnTM'])
 mfiltered_len=sum(i for _,i in total['GCnTM'])
# vgcn = total_vgcn / float(total_len)
 if len(rox[C])>1:
  wrl= [description[C],str(vecn),str(mecn),str(necn),str(vgcn/float(filtered_len)), str(mgcn/float(mfiltered_len)), str((vgcn-mgcn)/float(filtered_len-mfiltered_len))]
  fh.write(",".join(wrl)+'\n')
  print str(sum(total_cux.values()))+'\t'+description[C], " \t{0:.1f}\t{1:.1f}\t{2:.1f}\t{3:.1f}".format(vecn,mecn,necn,vgcn*100.0/filtered_len)###,"{:.1f}".format(EAA(cux))
fh.close()
print ' ECN written to:\033[96m',ctable_name,'\033[0m'

# may easily do codon usage for any gene or combination thereof - consider saving more than just genome-wide
ctable = total_cux #this depends on the order of cases, make sure genome is last for this to work
utable_name = name[0]+'_'+ttable_str+'codon_usage.csv'
fh = open(utable_name,'w')
utable = CodonTable.unambiguous_dna_by_id[ttable_id]
number_of_good_codons = 0
#number_of_good_codons = sum(ctable.values())
for xc in utable.forward_table.iteritems():
  xcs = str(xc[0])
  xcn = ctable.get(xcs,0)
  number_of_good_codons += xcn
  for_saving = xcs+','+str(xcn)+','+str(xc[1])+'\n'
  fh.write(for_saving)
#  if vecn <40 and xcn <5:
#    print xcs,xcn,xc[1]
fh.close()

del rox, ctable

print number_of_good_codons,'sense codons'
print ' codon usage table written to:\033[96m',utable_name,'\033[0m'

CDSr(new_rec)
trnas(new_rec)
print 'Individual genes written to:\033[96m',name[0]+'.CDS.fasta',name[0]+'.trn.fasta\033[0m'


print 'preparing graphics'

if good_map:
 presentation_rec1 = new_rec
 presentation_rec2 = phobius_rec
else:

 if verbose:print 'checking start'
 #check valid starting points for drawing..
 # make sure repeats don't cross the end
 # make a copy of the record and add bogus features covering the span of each repeat.
 # span with bogus feature is created at v2gff.py stage

 bogus_repeats=[]
 for f in rec.features:
   if f.type=='bogus':
     bogus_repeats.append(f)
 test_rec = copy.copy(new_rec)
# This is problematic for repeat-rich genomes. Need serious troubleshooting before production.
# turned off for now.
# test_rec.features.extend(bogus_repeats) 

 valid_loc = []
 test_f_list = sorted(new_rec.features, key = lambda x: ending(x), reverse = False)
 e_f_dict = {}
 e_list = []
 for f in test_f_list:
  e = ending(f)
  e_list.append(e)
  try:
   e_f_dict[e].append(f)
  except KeyError:
   e_f_dict[e] = [f]
 e_list = sorted(list(set(e_list)))
 for x in e_list:
  good = True
  for f in test_f_list:
   if f not in e_f_dict[x] and x in f:
    good = False
    break
  if good:
    valid_loc.append(FeatureLocation(x,x+1,None))  

 valid_loc = valid_loc + raw_ncr_locations(test_rec)
 valid_loc = amalgamate_variable_list(valid_loc)
 valid_loc = sorted(valid_loc, key = lambda x: x.start, reverse = False)
 if len(valid_loc) >1:
  valid_start = SeqFeature(CompoundLocation(valid_loc),'misc',{})
 elif len(valid_loc) == 1:
  valid_start = SeqFeature(valid_loc[0],'misc',{})
 else: # no valid starts, the one we have will have to do
  valid_start = SeqFeature(FeatureLocation(0,1,1),'misc',{})
  print 'failed to locate truly good drawing start'
# valid_loc = amalgamate([valid_start],0)
 del test_rec, bogus_repeats
 if verbose:
   print 'checked'
   print '# of valid start regions:',len(valid_loc)
 # rotate presentation if desired
 if 0 in valid_start:
  presentation_rec1 = new_rec
  presentation_rec2 = phobius_rec
 else:
  print 'start inside a feature, drawing will be rotated',
#  point_up = midpoint(valid_loc[0])
#  point_down = midpoint(valid_loc[-1])
  point_up = valid_loc[0].start
  point_down = valid_loc[-1].end - 1
# minimal rotation to allow seamless drawing
  if sequence_len - point_down < point_up:
   offset_point = point_down
   print sequence_len - point_down,'bp left, new start at', offset_point
   delta = sequence_len - point_down
  else:
   offset_point = point_up
   print offset_point, 'bp right'
   delta = -1*point_up
  presentation_rec1 = SeqRecord(new_rec.seq)
  presentation_rec2 = SeqRecord(new_rec.seq)
  presentation_rec1.annotations = new_rec.annotations
  presentation_rec1.letter_annotations = new_rec.letter_annotations
  presentation_rec1 = presentation_rec1[offset_point:] + presentation_rec1[0:offset_point]
  presentation_rec2 = presentation_rec2[offset_point:] + presentation_rec2[0:offset_point]

  for f in new_rec.features:
   if deepverbose: print f.type, f.location
   fx = SeqFeature(offset_universal(f.location,delta),f.type,{})
   fx.qualifiers.update(f.qualifiers)
   if deepverbose: print fx.location, 'copied from', f.location, fx.location.end - f.location.parts[-1].end, fx.location.start - f.location.parts[0].start
   presentation_rec1.features.append(fx)

  for f in phobius_rec.features:
   fx = SeqFeature(offset_loc(f.location,delta),f.type,{})
   fx.qualifiers.update(f.qualifiers)
   presentation_rec2.features.append(fx)
  if verbose:print 'rotated'

if verbose: print 'stats'

# per-domain stats in pseudo-window mode.

#tm_ECN =[]
#for f in presentation_rec2.features:
# coordinate=midpoint(f.location)
# rx=slice_rec(presentation_rec2,f)
# rx.features[0].type='CDS'
# rxtable=codon_usage(rx)
# print len(rxtable)
# ecn=ECNa(rxtable)
# print ecn
# tm_ECN.append((coordinate,ecn))
#print tm_ECN

### windowed ECN and other stats
fix_wdw = 160
if verbose: print 'preparing set',fix_wdw
rec_set = window_set(presentation_rec1, fix_wdw)
#re-create CDS features, including partials, from "codon" letter annotations, for each window
#this is needed for ECN
# since windowed ECN stats are rather unreliable consider removing this.
# Perhaps save stats earlier with coordinates 
# and make the window-independent graphs out of them before windowed stats.

for ea in rec_set: 
 ea.features=[]
 beg = 0
 while beg < len(ea)-3:
  frame = 0
  for i in range(beg, len(ea)-3):
   if ea.letter_annotations['codon'][i] == 1 and ea.letter_annotations['codon'][i+1] == 2:## and ea.letter_annotations['codon'][i+3] == 3:
    frame = +1
    break
   if ea.letter_annotations['codon'][i] == 3 and ea.letter_annotations['codon'][i+1] == 2:## and ea.letter_annotations['codon'][i+3] == 1:
    frame = -1
    break
  if frame:
   start = i
   for j in range(i, len(ea),3):
    if ea.letter_annotations['codon'][j] <> ea.letter_annotations['codon'][start]:
     break
   stop = j
   f = SeqFeature(FeatureLocation(start,stop,frame),'CDS',{})
   ea.features.append(f)
   beg = j+1
  else:
   beg = len(ea)+1
if verbose: print 'windowed stats'
## calculate some windowed stats with the same set of records (fix_wdw!)
## TODO: check if all are actually used
wdw_ECN = []
#wdw_05 = []
at_skew_bars = []
at_neutral = []
at_non_neutral = []
gc_skew_bars = []
gc_con = []
gc_non_neutral = []
gc_n = []
sch =[]


#if have_coverage:
# Lcov = []
for ea in rec_set:
 coordinate = int(ea.id)
 cux = codon_usage(ea)
 valueE = ECNa(cux)
 wdw_ECN.append((coordinate,valueE))
# wdw_ECN.append((coordinate,0))
#Coverage here: use rec_set (really, any will do).
# and either just store the letter annotation from the middle or do some averaging,
#based on the number of records in the set and the length of sequence.
# lets start with just the middle value.
# if have_coverage:
#   valueC = float(ea.letter_annotations['coverage'][len(ea)/2]) # this needs to be optimized
#   Lcov.append((coordinate, math.log(valueC)))
# scale5 = 0.5
# wdw_05.append((coordinate,scale5)) #needed for debugging only?

 A = ea.seq.upper().count('A')
 T = ea.seq.upper().count('T')
 G = ea.seq.upper().count('G')
 C = ea.seq.upper().count('C')
# valueAT_S = float(A - T)/(A + T)
 valueAT_S = safe_divide((A - T),(A + T))

 valueGC = float(G + C)/fix_wdw 
#This will not be accurate if many Ns...Consider implementation of a checkpoint.

 valueGC_S = safe_divide((G - C),(G + C))

 gc_con.append((coordinate, valueGC))
 at_skew_bars.append((coordinate, valueAT_S))
 gc_skew_bars.append((coordinate, valueGC_S))

 An,d = r_con(ea,[0,1],['A'],'prop')
 Tn,d = r_con(ea,[0,1],['T'],'prop')
 valueATN = safe_divide((An-Tn),(An + Tn))
 at_neutral.append((coordinate, valueATN))

 As,d = r_con(ea,[0,2],['A'],'codon')
 Ts,d = r_con(ea,[0,2],['T'],'codon')
 valueATS = safe_divide((As-Ts),(As + Ts))
 at_non_neutral.append((coordinate, valueATS))

 Gs,d = r_con(ea,[0,2],['G'],'codon')
 Cs,d = r_con(ea,[0,2],['C'],'codon')
 valueGCS = safe_divide((Gs-Cs),(Gs + Cs))
 gc_non_neutral.append((coordinate, valueGCS))

 Gn,d = r_con(ea,[0,1],['G'],'prop')
 Cn,d = r_con(ea,[0,1],['C'],'prop')
 valueGCN = safe_divide((Gn+Cn),d)
 gc_n.append((coordinate, valueGCN))

 sch.append((coordinate,shannon(ea.seq)))

## Some of these could be saved here as raw data for further analysis/alternate plotting.
##  augmented with matching GCn, easy to add GC3 if needed/desired!

## some clever merging, not sure if needed...First convert both to dict, 
 # then iterate over the first one and retrieve both values with the same key.
 #finally sort and save
if verbose: print 'saving windowed ECN'
d1=dict(gc_n)
d2=dict(wdw_ECN)
dataf=[]
for it in d1:
 dataf.append([it,d1[it],d2[it]])
dataf=sorted(dataf)
fh=open(name[0]+'.wdw_ECN.csv','w')
for it in dataf:
 line = ','.join(str(x) for x in it)
 fh.write(line+'\n')
fh.close()
if verbose: print 'saved'

if verbose: print 'modifiyng wdw ECN'
## NOW, back to wdw_ECN: replace zeroes with the median of the rest!
# first calculate the statistics. Bulid a list of non-zero values:

x_zero = min([x[1] for x in wdw_ECN])
non_zero_list =[]
for x in wdw_ECN:
 if x[1] > x_zero:
  non_zero_list.append(x[1])
# calculate median (and other stats):
x_med = median(non_zero_list)
x_min = min(non_zero_list)
x_max = max(non_zero_list)
#replace minimals with median:
new_ECN = []
while len(wdw_ECN) >0:
 x = wdw_ECN.pop()
 coordinate = x[0]
 valueE = x[1]
 if x[1] == x_zero:
   valueE = x_med
 new_ECN.append((coordinate,valueE))
wdw_ECN = sorted(new_ECN)


## two other windowed sets, no dependence on the features implemented
mid_wdw = 500
if verbose: print 'preparing set', mid_wdw
at_neu2 = []
gc_neu2 = []
rec_set = window_set(presentation_rec1, mid_wdw)

if verbose: print 'windowed stats'
for ea in rec_set:
 coordinate = int(ea.id)
 A,d = r_con(ea,[0,3],['A'],'codon')
 T,d = r_con(ea,[0,3],['T'],'codon')
 valueAT = float(A-T)/(A + T)
 at_neu2.append((coordinate, valueAT)) #AT_skew at third codon position
 G,d = r_con(ea,[0,3],['G'],'codon')
 C,d = r_con(ea,[0,3],['C'],'codon')
 valueAT = float(G-C)/(G + C)
 gc_neu2.append((coordinate, valueAT)) #GC_skew at third codon position

lrg_wdw = 1000
if verbose: print 'preparing set', lrg_wdw
rec_set = window_set(presentation_rec1, lrg_wdw)
at_neu3 = []
gc_neutral = []
gc_neu3 = []
at_neutral_long = []

if verbose: print 'windowed stats'
for ea in rec_set:
 coordinate = int(ea.id)
 A,d = r_con(ea,[0,1],['A'],'prop')
 T,d = r_con(ea,[0,1],['T'],'prop')
 G,d = r_con(ea,[0,1],['G'],'prop')
 C,d = r_con(ea,[0,1],['C'],'prop')
 at_neu3.append((coordinate, float(C-T)/(C+T))) #CT_skew at neutral sites, consider renaming
 gc_neutral.append((coordinate, float(G-C)/(G+C))) #GC_skew at neutral sites
 gc_neu3.append((coordinate, float(G-A)/(G+A))) #GA_skew (?) at neutral sites, consider renaming
 at_neutral_long.append((coordinate, float(A-T)/(A+T))) #AT_skew at neutral sites

if verbose: print 'stats prepared'

if verbose: print 'normalizing'
scaleAT = scaling_factor([at_skew_bars,at_neutral,at_non_neutral,at_neu2,at_neu3])
at_skew_bars = normalize_values(at_skew_bars,scaleAT)
at_neutral = normalize_values(at_neutral,scaleAT)
at_non_neutral = normalize_values(at_non_neutral, scaleAT)
at_neu2 = normalize_values(at_neu2, scaleAT)
at_neu3 = normalize_values(at_neu3, scaleAT)
at_neutral_long = normalize_values(at_neutral_long, scaleAT)
#at_05 = normalize_values(wdw_05,scaleAT)

scaleGC = scaling_factor([gc_skew_bars,gc_neutral,gc_non_neutral,gc_neu2,gc_neu3])
gc_skew_bars = normalize_values(gc_skew_bars,scaleGC)
gc_neutral = normalize_values(gc_neutral,scaleGC)
gc_non_neutral = normalize_values(gc_non_neutral, scaleGC)
gc_neu2 = normalize_values(gc_neu2,scaleGC)
gc_neu3 = normalize_values(gc_neu3,scaleGC)

max_gc = scaling_factor([gc_con,gc_n])/3 ###!!
gc_con = normalize_values(gc_con,max_gc,gc_c/sequence_len)
gc_n = normalize_values(gc_n,max_gc,gc_c/sequence_len)

max_ECN = (x_max-x_min)/2
wdw_ECN = normalize_values(wdw_ECN,max_ECN,x_med)
#tm_ECN = normalize_values(tm_ECN,50,30)

if verbose: print 'drawing'
#two diagrams: circular full map and "pub figure (clean)", with only selected objects

diagram = GenomeDiagram.Diagram('map')
pub_diagram = GenomeDiagram.Diagram('figure')

s_track = diagram.new_track(4, name='Scale-only', scale_smalltick_interval = 200, scale_smallticks = -1.0, height = 8, scale_largetick_interval = 1000, scale_largetick_labels = True, scale_format = 'SInt',scale_fontsize=9)
f_track = diagram.new_track(6, name='Features', scale = 0, height = 8)
a_track = diagram.new_track(3, name='AT-skew',scale = 0, height = 40)
g_track = diagram.new_track(2, name='GC-skew',scale = 0, height = 40)
c_track = diagram.new_track(1, name='GC-contents',scale = 0, height = 40)
r1_track = diagram.new_track(5, name='repeats',scale = 0, height =2)
r2_track = diagram.new_track(7, name='repeats',scale = 0, height =0)
#r_empty = diagram.new_track(9, name='repeats',scale = 0, height =20)
ss_track = pub_diagram.new_track(4, name='Scale-only', scale_smalltick_interval = 200, scale_smallticks = -1.0, height = 4, scale_largetick_interval = 1000, scale_largetick_labels = True, scale_format = 'SInt',scale_fontsize=10)
ff_track = pub_diagram.new_track(6, name='Features', scale = 0, height = 6)
aa_track = pub_diagram.new_track(3, name='AT-skew',scale = 0, height = 20)
gg_track = pub_diagram.new_track(2, name='GC-skew',scale = 0, height = 20)
cc_track = pub_diagram.new_track(1, name='GC-contents',scale = 0, height = 15)
oo_track = pub_diagram.new_track(7, name='empty',scale = 0, height = 3)
rr_track = pub_diagram.new_track(8, name='bogus',scale=0, height=0)

blue_rainbow = cycle(['aquamarine','turquoise','lightseagreen','teal','darkcyan','darkslategray'])
red_rainbow = cycle(['rosybrown','lightcoral','indianred','brown','firebrick','sienna', 'saddlebrown','maroon','darkred'])
dark1 = colors.Color(0.5,0.5,0.5,alpha = 0.25)
good_color=colors.Color(0.8,0.8,0.8,alpha=1)
bad_color=colors.Color(1,1,1,alpha=0) #just dont draw it, there must be a better way to do it... 
dark = colors.CMYKColor(0,0,0,0.3)
darker = colors.CMYKColor(0,0,0,0.9)
darkest = colors.Color(0,0,0)
reddish = colors.CMYKColor(0,0.2,0.1,0.1)
blueish = colors.CMYKColor(0.9,0.5,0,0.2)
ff_set = ff_track.new_set() 
# this is a set primarily for holding annotations but pA will be under/over drawn on it as a graph
# and so will be repeat boxes. Thats why it is defined here already (repdraw follows)...

##REPDRAW
bogus_set= rr_track.new_set()
c_backset = r1_track.new_set()
for f in presentation_rec1.features:
    if f.type == 'repeat_region':
        nice = 0
        if f.qualifiers['note'] == 'inverted':
           nice = next(red_rainbow)
        elif f.qualifiers['note'] == 'direct':
           nice = next(blue_rainbow)
        for loc in f.location.parts:
         loc.strand=None
         g= SeqFeature(loc)
         c_backset.add_feature(g, label=False, color=nice, sigil='BOX')
#it is better to add each part of the repeat as a separate feature, 
#apparently the ones crossing start are not drawn correctly.
#
# For clean draw two crosslinks for each repeat.
# first create bogus location between them.
# create crosslinks from this location on empty outer track to both location parts of the repeat
# nothing crossing start is drawn anyway, so this is not going to be acceptable for presentation
# unless, of course the condition is recognized and taken care of by rotation
# use color to indicate reverse repeats

        if len (f.location.parts[0]) > 2*word_size:
         rep_color = dark
         if f.qualifiers['note']== 'inverted':
           rep_color = reddish
         ff_set.add_feature(SeqFeature(f.location.parts[0],'misc',{}), label=False, color = bad_color, border = rep_color, sigil='BOX')
         ff_set.add_feature(SeqFeature(f.location.parts[-1],'misc',{}), label=False, color = bad_color, border = rep_color, sigil='BOX')
         rep_len=max(len(f.location.parts[0]),len(f.location.parts[-1]))
         rep_separation = f.location.parts[-1].start - f.location.parts[0].end
# draw crosslink only if parts separated by a distance
         if rep_separation > 0 and rep_separation > rep_len/2:
          x1 = midpoint(f.location.parts[0])
          x2 = midpoint(f.location.parts[-1])
          bogus = (x2+x1)/2
          pub_diagram.cross_track_links.append(CrossLink((rr_track,bogus-3,bogus),(ff_track, f.location.parts[0].end-3, f.location.parts[0].end),bad_color,rep_color))
          pub_diagram.cross_track_links.append(CrossLink((rr_track,bogus,bogus+3),(ff_track, f.location.parts[-1].start, f.location.parts[-1].start+3),bad_color,rep_color))

# consider making several layers of crosslinks if the list is long.


#r2_track.add_set(c_backset)

##polyA sites presentation
pA_sites=[]
for f in presentation_rec1.features:
 if 'polyA' in f.type:
  pA_sites.append(f)
if len(pA_sites)>0:
 pA_set=f_track.new_set('graph')
 pA_tbl_plus=[]
 pA_tbl_minus=[]

 m=1.95
 
# how much sticks out of the track (1=perfect fit, may not be even seen in concordant config). 
## Primary reason of doing this as a this convoluted graph. ##

 for f in pA_sites:
  X=int(f.location.start)
  if X > 0:
   X1 = X - 1
  else:
   X1 = X - 1 + sequence_len
  if X + 1 < sequence_len:
   X2 = X + 1
  else:
   X2 = X + 1 - seqence_len
  Y=-1*f.strand # change if want it on the other side.
#  Y=f.strand
  if Y==1:
   pA_tbl_plus.append((X1,-1))
   pA_tbl_plus.append((X,Y*m))
   pA_tbl_plus.append((X2,-1))
  else:
   pA_tbl_minus.append((X1,1))
   pA_tbl_minus.append((X,Y*m))
   pA_tbl_minus.append((X2,1))
 if len(pA_tbl_plus)>0:
  pA_set.new_graph(pA_tbl_plus,'',style='bar',colour=blueish,altcolour=bad_color,centre=0)
 if len(pA_tbl_minus)>0:
  pA_set.new_graph(pA_tbl_minus,'',style='bar',colour=bad_color,altcolour=blueish,centre=0)
 ff_track.add_set(pA_set) #copy pA to clean display.
del pA_sites

f_set = f_track.new_set()

nice = 0 # not sure why does it work. Causes black outline and white fill. 
R,G,B,A=1,0.65,0,0.5
orange=colors.Color(R,G,B,A)
K=1-max(R,G,B)
C=(1-R-K)/(1-K)
M=(1-G-K)/(1-K)
Y=(1-B-K)/(1-K)
orange1=colors.CMYKColor(C,M,Y,K)
del A,R,G,B,C,M,Y,K

## Domains (very simple but efficient)
for f in presentation_rec2.features:
  f_set.add_feature(f,label=False,colour=orange,border = False,sigil='OCTO')
  ff_set.add_feature(f,label=False,colour=orange1,border=False, sigil='OCTO')

## Other features
#consider changing all colors to non-transparent CMYK for easier eps conversion.
# obligatory for clean graph
# TO DO: full contol of both fill and stroke color and on the order of drawing.
for f in presentation_rec1.features:
 if f.type in ['CDS', 'tRNA', 'rRNA', 'misc']:
  shape = 'ARROW'
  bcolor=colors.Color(0,0,0)
  nice = bad_color
  try:
   lbl = f.qualifiers['product']
  except:
    try:
     lbl = f.qualifiers['label']
    except:
#     print f
     pass
  lbl_pos = 'middle'
  if 'RNA' in f.type:
   nice =  1
   bcolor=colors.CMYKColor(0,0,0,0.8)
   lbl = f.qualifiers['label']
   try:
    lbl=lbl+f.qualifiers['note']
   except:
    pass
  xxx = lbl
  if f.type == 'misc':
    xxx= f.qualifiers['label']
  if len(f) < 130:
   lbl_pos = 'bottom'
  small_label=9
  big_label=11
  f_set.add_feature(f, label=True, label_size = small_label, name = lbl, label_position = lbl_pos, label_strand = 1, color=nice, border=bcolor, sigil=shape, arrowshaft_height=1.0)
#  if f.type <> 'misc':
  ff_set.add_feature(f, label=True, label_size = big_label, name = xxx, label_position = lbl_pos, label_strand = 1, color=nice, border=bcolor, sigil=shape, arrowshaft_height=1.0)

## AT-skew stats
a_set = a_track.new_set('graph')
aa_set = aa_track.new_set('graph')

a_set.new_graph(at_skew_bars,'AT skew', style = 'bar',centre=0)
a_set.new_graph(at_non_neutral,'AT non-neutral', style = 'line',colour = colors.red,centre=0)
a_set.new_graph(at_neutral,'AT neutral', style = 'line',colour = colors.black,centre=0)
a_set.new_graph(at_neu2,'AT neutral2', style = 'line',colour = colors.green,centre=0)
a_set.new_graph(at_neu3,'AT neutral3', style = 'line',colour = colors.blue,centre=0)

aa_set.new_graph(at_skew_bars,'AT skew', style = 'bar',centre=0)
aa_set.new_graph(at_non_neutral,'AT non-neutral', style = 'line',colour = colors.red,centre=0)
aa_set.new_graph(at_neutral_long,'AT neutral', style = 'line',colour = colors.black,centre=0)

## GC-skew stats
g_set = g_track.new_set('graph')
gg_set = gg_track.new_set('graph')

first_color = colors.lightblue
second_color = colors.blue
if skewc < 0:
 first_color = colors.blue
 second_color = colors.lightblue
g_set.new_graph(gc_skew_bars, 'GC skew', style = 'bar', colour=first_color, altcolour=second_color,centre=0)

g_set.new_graph(gc_non_neutral,'GC non-neutral', style = 'line', colour=colors.red,centre=0)
g_set.new_graph(gc_neutral,'GC neutral', style = 'line', colour=colors.black,centre=0)
g_set.new_graph(gc_neu2,'GC neutral2', style = 'line',colour = colors.green,centre=0)
g_set.new_graph(gc_neu3,'GC neutral3', style = 'line',colour = colors.blue,centre=0)

gg_set.new_graph(gc_skew_bars, 'GC skew', style = 'bar', colour=first_color, altcolour=second_color,centre=0)
gg_set.new_graph(gc_neutral,'GC neutral2', style = 'line',colour = colors.green,centre=0)
#gg_set.new_graph(gc_5,'gray 0.5',style='line',colour = colors.gray,centre =0)

## supplementary orf for map added to scale track in map only!
s_set = s_track.new_set()
for f in presentation_rec1.features:
 if f.type == 'orf':
   s_set.add_feature(f, label=False, colour=dark1,border=False, sigil='ARROW')

## shannon entropy and ECN added to composition track in map only (dirty hack).
#  otherwise - composition (GC content) and GCn added to both map and figure, same track
# TO DO: consider simlifying by adding the same set to both tracks.

s_gr = c_track.new_set('graph')
sch_n = normalize_values(sch,0.1,1.2)
shannon_graph= c_track.new_set('graph')
shannon_graph.new_graph(sch_n,'entropy',style='line',colour=dark1,centre=0)

c_set = c_track.new_set('graph')
cc_set = cc_track.new_set('graph')

c_set.new_graph(gc_con,'GC', style = 'bar', colour=colors.lightcoral, altcolour=colors.red,centre=0)
cc_set.new_graph(gc_con,'GC', style = 'bar', colour=colors.lightcoral, altcolour=colors.red,centre=0)

s_gr.new_graph(wdw_ECN,'ECN', style = 'line',colour = colors.darkblue,centre = 0)
c_set.new_graph(gc_n,'GCn', style = 'line', colour = colors.green,centre = 0)
cc_set.new_graph(gc_n,'GCn', style = 'line', colour=colors.green, centre=0)

# adding graphics and feature sets to the same track is possible: f_track.add_set(c_set)
# sets can be created first, added to track later. No sure how to define their type then.
# anyway, separate set for polyA sites can be created as graphics and overlaid with main features.
# for that need simple list of tuples consisting of coordinates and arbitrary value (+/-1)
# may do the same trick for coverage.

#Coverage here: use rec_set (really, any will do).
# and either just store the letter annotation from the middle or do some averaging,
# based on the number of records in the set and the length of sequence (?) 
# hint: look in wdw procedure for step
# lets start with just the middle value.
# consider scale in the background

if have_coverage:
 Lcov = []
 mid_idx=len(rec_set[0])/2
 step=(int(rec_set[2].id)-int(rec_set[1].id))/2 # should be = sequence_len/2200, so in the range of 10
 x1=mid_idx-step
 x2=mid_idx+step+1
 step=2*step
 for ea in rec_set:
  coordinate = int(ea.id)
  valueC = math.fsum(ea.letter_annotations['coverage'][x1:x2])/step
#  valueC = ea.letter_annotations['coverage'][mid_idx]/step
  Lcov.append((coordinate, math.log10(valueC)))

# syntax reminder:
# def normalize_values(value_list,max_value,ofset=0):(y-ofset)/max
 max_cov = scaling_factor([Lcov])
# now prep the scale. This is log10 so one order of magnitude should be adequate tick 
 scale1=int(max_cov)
 s1=[]
 s2=[]
 s3=[]
 for coordinate,value in Lcov:
#  s1.append((coordinate,scale1))
  s2.append((coordinate,scale1-1))
  s3.append((coordinate,scale1-2))
 Lcov= normalize_values(Lcov,max_cov/3.0,max_cov/1.7)
#use the same transformation for scale
# s1 = normalize_values(s1,max_cov/3.0,max_cov/1.7)
 s2 = normalize_values(s2,max_cov/3.0,max_cov/1.7)
 s3 = normalize_values(s3,max_cov/3.0,max_cov/1.7)
 cov_set = f_track.new_set('graph')
 cov_graph = cov_set.new_graph(Lcov,'coverage',style='line',colour=blueish,centre=0)
 zz_track = pub_diagram.new_track(9, name='outer',scale=0, height=25)

 scale_set = zz_track.new_set('graph')
# back_scale = scale_set.new_graph(s1,'',style='line',colour=dark,centre=0)
 back_scale = scale_set.new_graph(s2,'',style='line',colour=dark,centre=0)
 back_scale = scale_set.new_graph(s3,'',style='line',colour=dark,centre=0)
 zz_track.add_set(cov_set) #order of adding sets to tracks determines order of drawing (last added, last drawn)

aa,bb= diagram.range() ## NOT sure why this no longer works... 
#if bb == sequence_len: #Anyway - this was a debugging failsafe, may switch it off
if True:
 diagram.draw(format='circular', circular = True, circle_core = 0.3)
 pub_diagram.draw(format='circular', circular = True, circle_core = 0.35)
 diagram.write(pdf_name, 'PDF')
 pub_diagram.write(clean_pdf_name, 'PDF')
 print 'graphics written to:\033[96m', pdf_name, clean_pdf_name, '\033[0m'
else:
 print 'format problem, invalid graph not saved'

if not good_map: # this is to prevent recurring shifts, no new gb will be created if edited gbk was used
 new_rec.id = name[0]
 new_rec.name = name[0]
 new_rec.annotations['organism']='similar to '+taxid
 new_rec.annotations['taxonomy']=lineage
 new_rec.annotations['topology']='circular'
 new_rec.annotations['data_file_division']=rec.annotations['data_file_division']
 new_rec.annotations['date']=date.today().strftime("%d-%b-%Y").upper()
 new_rec.features.extend(phobius_rec.features)
 SeqIO.write(new_rec, output_name, "genbank")
 print 'annotated genome written to:\033[96m', output_name,'\033[0m'

coding = ['CDS','rRNA','tRNA']
with open(name[0]+'.crex','w') as fh:
 fh.write(">"+new_rec.id+'\n')
 for f in new_rec.features:
           if f.type not in coding:
             continue
           else:
            if f.type == 'CDS':
              allowed_tags = ['product']
            else:
              allowed_tags = ['label']
            for key, vals in f.qualifiers.iteritems():
              if key not in allowed_tags:
                continue
              else:
                for v in vals:
                   if f.strand == 1:
                     fh.write(v)
                   else:
                     fh.write("-"+v)
                fh.write(' ')
 fh.write('\n')             
print 'gene order for CREx written to:\033[96m', name[0]+'.crex\033[0m'
if verbose:print 'the end'

