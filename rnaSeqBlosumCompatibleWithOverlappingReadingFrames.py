#! /usr/bin/env python

# rnaSeqBlosumCompatibleWithOverlappingReadingFrames.py
# P.Clote


# This program implements breadth first search BFS to compute the following.
#
# Given peptide1,peptide2 this program returns all mRNA that APPROXIMATELY translate peptide1 (reading frame 0)
# and APPROXIMATELY translate peptide2 (reading frame +1). In contrast to the related program
#        rnaSeqCompatibleWithOverlappingReadingFrames.py
# this program returns mRNAs that translate into peptides, such that each amino acid of the peptide has Blosum62
# similarity at least the user-specified threshold with the input peptide1 [resp. peptide2]



import sys,os,copy
from aminoAcidAndGeneticCodes import genCode,aaCode, singleLetterAAcodes
from extractHIVgagPolJunction import translateMRNA

NUCL     = ['A','C','G','U']


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

def createListOf4mersForAllPairsOfResidues(filename=None,threshold=1):
  C = aa2codonList(filename,threshold)
  #D[dimer] is list of codon pairs that agree on overlap
  D = {}
  for aa1 in singleLetterAAcodes:
    for aa2 in singleLetterAAcodes:
      L  = []
      for x in C[aa1]:
        for y in C[aa2]:
          if x[1:]==y[:2]:
            L.append( x+y[-1] )
      D[aa1+aa2] = L
  return D

def consistent(fourMer1,fourMer2):
  if fourMer1[-1]==fourMer2[0]:
    return True 
  else:
    return False

def bfs(peptide1,peptide2,filename=None,threshold=1):
  #WARNING: This assumes that peptide1 in reading frame 0 and 
  #peptide2 in reading frame 1.
  #WARNING: Assume len(peptide1) = len(peptide2) 
#  assert (len(peptide1) == len(peptide2))
  lenPeptide = min(len(peptide1),len(peptide2)) 
  D     = createListOf4mersForAllPairsOfResidues(filename,threshold)
  dimer = peptide1[0]+peptide2[0]
  #terminate with list of RNAsequences of len 3*lenPeptide+1
  RNA   = D[dimer]
  for i in range(1,lenPeptide):
    dimer = peptide1[i]+peptide2[i]
    L     = D[dimer] 
    RNAbis= []
    for rna in RNA:
      for fourMer in L:
        if rna[-1]==fourMer[0]:
          RNAbis.append(rna+fourMer[1:])
    RNA = copy.deepcopy(RNAbis)
  return RNA       

def check(peptide1,peptide2,RNA):
  for rna in RNA:
    if translateMRNA(rna[:-1]) != peptide1:
      return False
    if translateMRNA(rna[1:]) != peptide2:
      return False
  return True

def printRNA(RNA):
  for rna in RNA:
    print rna

def main(peptide1,peptide2,filename=None,threshold=1):
  minLen = min(len(peptide1),len(peptide2))
  RNA    = bfs(peptide1,peptide2,filename,threshold)
  print "Num seq: %s\tLen seq: %s" % (len(RNA),len(RNA[0]))
  printRNA(RNA)
  print "Check translation of gag and pol"
  for i in range(1,len(RNA)+1):
    rna = RNA[i-1]
    print "%s\t%s" % (translateMRNA(rna),translateMRNA(rna[1:]))
#  print check(peptide1,peptide2,RNA)


if __name__ == '__main__':
  if len(sys.argv) < 3:
    print "Usage: %s peptide1 peptide2 [Blosum62matrix [threshold]]" % sys.argv[0]
    print "  1) Peptides given in IUPAC single letter codes"
    print "  2) Filename of amino acid (Blosum62) similarity matrix (default None)"
    print "  3) Threshold in Blosum62 for codons that translate to residues with >= threshold similarity (default +1)"
    sys.exit(1)
  peptide1 = sys.argv[1]
  peptide2 = sys.argv[2]
  if len(sys.argv)>3:
    filename = sys.argv[3]
  else:
    filename = None
  if len(sys.argv)>4:
    threshold = int(sys.argv[4])
  else:
    threshold = 1
  main(peptide1,peptide2,filename,threshold)

