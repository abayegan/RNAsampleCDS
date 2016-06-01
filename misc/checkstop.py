#! /usr/bin/env python
import os,sys
from aminoAcidAndGeneticCodes import *

def parsefasta(fastafile):
	d={}
	fin = open(fastafile,'r')
	line = fin.readline().strip()
	if not line:
		print "Error in fasta file!"
		sys.exit(1)
	seqid = line[1:]
	line = fin.readline().strip().replace('-','').replace('T','U')
	seq=''
	while line:
		if line[0]!='>':
			seq += line
		elif line[0]=='>':
			if seqid in d:
				print 'sequences in the fasta file must have unique ids!'
				sys.exit(1)
			d[seqid]=seq.upper().replace('-','').replace('T','U')
			seq=''
			seqid = line[1:]
		line = fin.readline().strip()
	d[seqid] = seq.upper().replace('-','').replace('T','U')
	return d
	
def translateMRNA(mrna,k):
  #WARNING: mrna is 0-indexed
  aaSeq = ""
  i     = 0
  while i<len(mrna)-1:
    codon  = mrna[i:i+3].replace('T','U')
    if len(codon)<3:
      #print aaSeq
      #print len(aaSeq)
	  #print "Untranslated nucleotides: %s" % codon
      break
    gcode = genCode[codon].upper()
    #check for premature abortion of translation, which occurs in PR55
    #gag polyprotein
    if gcode == 'STOP':
      #print gcode,i,mrna[i:i+3]
      #print aaSeq
      print 'stop codon in %s!' %k
      break
    aa     = aaCode[gcode]
    aaSeq += aa
    i     += 3
  return aaSeq

d = parsefasta(sys.argv[1])
for k,v in d.iteritems():
	translateMRNA(v,k)
	translateMRNA(v[1:],k)
	
