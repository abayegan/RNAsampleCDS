#! /usr/bin/python2.6
#A. Bayegan
import os,sys,random
from aminoAcidAndGeneticCodes import *

def printusage(exe):
	txt ="""
	Usage:
	
	./RNAsampleCDS.py -p peptideFile -n numSamples -c constraintFile -m BlosumMatrix -t threshold -gc GC_content -w -pssm pssmOutFile -v

	Options:
	
		-p 
		Define the file name containing the peptides to be coded by the output RNA sequence. 	
		
		-n
		Define an integer as the number of RNAs printed in the output.  
	
		-c 
		Define nucleotide constraints on the coding RNA sequences. It is possible to define IUPAC code nucleotide constraint of length <= 3n+2 for the given peptides of length n. If the constraint length is m < 3n+2 the program assumes all the nucleotides from m to 3n+2 can code for any nucleotide(IUPAC code 'N').   
	
		-m 
		If the RNA sequences are not required to code exactly the given peptides, it is possible to define a similarity matrix between the amino acids and use it as in input to the program. Blosum62 similarity matrix is included in the program in a file titled "blosum62.txt". If it is used the output RNAs code peptides with blosum similarity greater or equal than 'threshold'(defined in -t flag) to the given peptides.
		
		-t
		Define an integer as the minimum threshold for the input similarity matrix. It must be used together with -m flag and guarantees the output RNAs code peptides with similarity greater or equal than 'threshold' to the given peptides.
	
		-gc
		Define a single integer value or an interval for GC-content. The input can be a single value or a range of format gcl:gcu where gcl and gcu are inclusive lower bounds.
	
		-w
		use weighted sampling 
	
		-pssm
		Define a file name to save the output logo plot. The logo plot is with computed for RNA sequences satisfying all other input options.
	
		-v
		For verbose output"""
	print (txt)
	sys.exit(1)			

def iscompatibleAAWithIupac(aa,iupac):
	#iupac codon has the regular 24 AA iupac codes one extra letters 'O'.
	# O represents any aa or 'stop'
	iscompat = False
	if aa == iupac:
		iscompat = True
	elif iupac=='B':
		if aa in ['D','N']:
			iscompat=True
	elif iupac=='X':
		if aa!='*':
			iscompat = True
	elif iupac=='Z':
		if aa in ['E','Q']:
			iscompat= True
	elif iupac=='O':
		iscompat=True
	return iscompat

def aa2codonList(filename=None,threshold=1):
	if filename!=None:
		Blosum62,AAs = readBlosumMatrix(filename)
	#initialize dictionary
	D = {}
	for aa in singleLetterAAcodes:
		D[aa] = []
	#create lists of codons for each aa
	for x in NUCL:
		for y in NUCL:
			for z in NUCL:
				codon = x+y+z
				if genCode[codon].upper() != 'STOP':
					aa    = aaCode[genCode[codon].upper()]
					if filename == None:
						D[aa].append(codon)
					else: #use Blosum similarity
						for aaBis in singleLetterAAcodes:
							if Blosum62[aa][aaBis]>=threshold:
								D[aaBis].append(codon)
	return D

def createListOfCodons():
	#create list L of all codons that are not STOP codons
	L = []
	for a in NUCL:
		for b in NUCL:
			for c in NUCL:
				codon = a+b+c
				if codon not in stopCodon:
					L.append(codon)
	return L

def createAllCodons():
	#create list L of all codons
	L = []
	for a in NUCL:
		for b in NUCL:
			for c in NUCL:
				codon = a+b+c
				L.append(codon)
	return L

def nuc2int(nuc):
	if nuc == 'A': return 0
	elif nuc == 'C': return 1
	elif nuc == 'G': return 2
	elif nuc == 'U' or nuc == 'T' : return 3
	else:
		print "error in nuc2int! Sequences must only contain A,C,G,U,T!"
		sys.exit(1)

def readBlosumMatrix(filename):
	file = open(filename)
	line = file.readline()
	AAs  = ['$']+line.split()  #Add leading bogus symbol to ensure 1-indexed
	#initialize dictionary
	D    = {}
	for aa in AAs: D[aa]={}
	line = file.readline()
	n    = 0
	while line:
		n    += 1
		words = line.split()
		aa    = words[0]
		assert (AAs[n]==aa)
		for k in range(1,len(words)):
			aa = AAs[k]
			D[words[0]][aa] = int(words[k])
		line = file.readline()
	file.close()
	return D,AAs

#given a dictionary as input sample from keys with probabilities in values
def rouletteWheel(d):
	if(sum(d.values())-1 > 0.001):
		print "error in roulette wheel! probability list must add up to 1."
		sys.exit(1)
	randnum = random.random()
	i=0;cumsum=0;n=len(d)
	for k,val in d.iteritems():
		cumsum += val
		if(randnum<cumsum): 
			x=k;
			break
	return x

def translateMRNA(mrna):
  #WARNING: mrna is 0-indexed
  aaSeq = ""
  i     = 0
  while i<len(mrna)-1:
    codon  = mrna[i:i+3].replace('T','U')
    if len(codon)<3:
      print aaSeq
      print len(aaSeq)
      print "Untranslated nucleotides: %s" % codon
      sys.exit(1)
    gcode = genCode[codon].upper()
    #check for premature abortion of translation, which occurs in PR55
    #gag polyprotein
    if gcode == 'STOP':
      #print gcode,i,mrna[i:i+3]
      #print aaSeq
      return aaSeq
    aa     = aaCode[gcode]
    aaSeq += aa
    i     += 3
  return aaSeq

def GCcont(seq):
	sum=0
	for ch in seq:
		if ch=='G' or ch=='C':
			sum += 1
	return sum

def consistent(fourMer1,fourMer2):
	if fourMer1[-1]==fourMer2[0]:
		return True 
	else:
		return False

def printRNA(RNA):
	for rna in RNA:
		print rna

def printFASTA(RNA):
	for k in range(len(RNA)):
		print "> %s" % (k+1)
		print RNA[k]
