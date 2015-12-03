#! /usr/bin/env python

# sampleRNAsWithSequenceConstraintsWithOverlappingReadingFrames.py
# P.Clote

# This program differs from the program
#    sampleRNAsWithOverlappingReadingFramesBlosumCompatibleWithGivenPeptides.py
# in several ways: (1) sequence constraints are supported, (2) only forward
# partition function ZF and statistically rigorous sampling are supported.
# The program
#    sampleRNAsWithOverlappingReadingFramesBlosumCompatibleWithGivenPeptides.py
# does NOT support sequence constraints, but it does compute the PSSM, supports
# sampling from PSSM as well as statistically rigorous sampling, computes the
# number of sequence that code in overlapping reading frames 0 and 1 (useful
# to compute the probability that a sequence code anything in 2 overlapping 
# reading frames).

import sys,os,copy,random
from math import log,exp
from aminoAcidAndGeneticCodes import genCode,aaCode, singleLetterAAcodes
from extractHIVgagPolJunction import translateMRNA
from rnaSeqBlosumCompatibleWithOverlappingReadingFrames import readBlosumMatrix
from rnaSeqBlosumCompatibleWithOverlappingReadingFrames import readBlosumMatrix

NUCL   = ['A','C','G','U']
STOP   = ['UAA', 'UGA', 'UAG']
DEBUG  = False
PRINT  = False #set False if use main as function to return list


def createListOfCodons():
  #create list L of all codons that are not STOP codons
  L = []
  for a in NUCL:
    for b in NUCL:
      for c in NUCL:
        codon = a+b+c
        if codon not in STOP:
          L.append(codon)
  return L


def createListOfOverlappingCodons():
  #create list L of ALL 4-tuples corresponding to overlapping codons
  CODONS = createListOfCodons
  L      = []
  for a in NUCL:
    for b in NUCL:
      for c in NUCL:
        codon1 = a+b+c
        if codon1 not in STOP:
          for d in NUCL:
            codon2 = codon1[1:]+d 
            if codon2 not in STOP:
              overlappingCodons = a+b+c+d
              L.append(overlappingCodons)
  return L


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

def createListOf4mersForAllPairsOfResidues0(filename=None,threshold=1):
  #C[aa] is list of all residues with Blosom similarity >= threshold to aa
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


def createListOf4mersForAllPairsOfResidues(filename=None,threshold=1):
  if filename!=None:
    Blosum,AAs = readBlosumMatrix(filename)
  #create list L of ALL 4-tuples corresponding to overlapping codons
  L = createListOfOverlappingCodons()
  #D[dimer] is list of codon pairs that agree on overlap
  D = {}
  for aa1 in singleLetterAAcodes:
    for aa2 in singleLetterAAcodes:
      Tuples = [] 
      for tuple in L:
        aa1bis = aaCode[genCode[tuple[:3]].upper()]
        aa2bis = aaCode[genCode[tuple[1:]].upper()]
        if filename==None:
          if aa1==aa1bis and aa2==aa2bis:
            Tuples.append(tuple)
        else: #use Blosum similarity; recall if aa==aabis then Blosum similar
          if Blosum[aa1][aa1bis]>=threshold and  Blosum[aa2][aa2bis]>=threshold:            Tuples.append(tuple)
      D[aa1+aa2] = Tuples
  return D

def consistent(fourMer1,fourMer2):
  if fourMer1[-1]==fourMer2[0]:
    return True 
  else:
    return False

def computeZF(peptide1,peptide2,D,constraint=None): 
  #forwards "partition function"
  # Here, peptide1 = p_{0}...p_{n-1} and peptide2 = q_{0}...q_{n-1}
  # Function computes the ZF[(k,ch)] = number of RNA sequences
  #    s = a_{0}...a_{3k} where len(s)=3*k+1  and such that
  # (1) translateMRNA(s) = peptide1[:k] 
  # (2) translateMRNA(s[1:]) = peptide2[:k]
  # (3) s[-1] = ch 
  # The final coding sequences are of the form a_{0} ... a_{3n}.
  #
  # Parameter explanation:
  # D is dictionary of 4-mer RNAs corresponding to peptide1 and peptide2
  # D is returned by function
  #    createListOf4mersForAllPairsOfResidues(filename,threshold)
  #
  #WARNING: This assumes that peptide1 is in reading frame 0 and 
  #         peptide2 is in reading frame 1.
  #WARNING: Current code assumes len(peptide1) = len(peptide2) 
  #         Future extension can easily remove this assumption
  #
  p = peptide1; q = peptide2
  assert (len(p) <= len(q)) 
  n  = min(len(p),len(q)) 
  #------------Initialize ZF---------------------
  ZF = {}
  for k in range(0,n+1):
    for ch in NUCL:
      ZF[(k,ch)] = 0.0
  if constraint==None or constraint[0]=='N':
    for ch in NUCL: ZF[(0,ch)] = 1.0
  else: #first position is constrained     
    ch = constraint[0]
    ZF[(0,ch)] = 1.0
  #------------Fill matrix ZF-------------------
  if constraint==None:
    for k in range(1,n+1): #k in [1,...,n] 
      dimer = p[k-1]+q[k-1]
      L     = D[dimer]
      for s in L:
        ZF[(k,s[-1])] += ZF[(k-1,s[0])]
  else: #must satisfy sequence constraint
    for k in range(1,n+1): #k in [1,...,n] 
      dimer            = p[k-1]+q[k-1]
      L0               = D[dimer]
      if len(constraint)<=3*(k-1):
        constraint4tuple = ''
      else:
        constraint4tuple = constraint[3*(k-1):3*k+1]
      L                = []
      for tuple in L0:
        ok = 1
        for i in range(len(constraint4tuple)):
          ch = constraint4tuple[i]
          if ch!='N' and tuple[i]!=ch:
            ok = 0
            #break
        if ok: 
          L.append(tuple)
      if DEBUG:
        print "computeZF ",k,L
      for s in L:
        ZF[(k,s[-1])] += ZF[(k-1,s[0])]
  return ZF


def printRNA(RNA):
  for rna in RNA:
    print rna

def sample(peptide1,peptide2,seqLen,numSamples,D,constraint):
  #D[4-tuple] is list of "similar" 4-tuplesaa1+aa2]
  #Here, the key 4-tuple codes aa1 in frame 0 and aa2 in frame 1
  #4-tuples in list L=D[4-tuple] code aa1' in frame 0 and aa2' in frame 1
  #where Blosum(aa1,aa1')>= threshold value, and Blosum(aa2,aa2')>= value.
  #WARNING: Do NOT need to worry about constraints here, since ZF handles this
  assert( len(peptide1)==len(peptide2) )
  SAMPLES = []
  ZF      = computeZF(peptide1,peptide2,D,constraint)
  for i in range(numSamples):
    #initialize empty sequence and set k=n
    k  = seqLen; rna = ""
    #--------choose random nucleotide for which ZF[(k,ch)]>0 
    ch = random.choice(NUCL)
    while ZF[(k,ch)]==0: ch = random.choice(NUCL)
    #--------construct sample
    while k>0:
      #create list of 4-tuples having similarity to p[k],q[k] in frame 0,1
      #and whose last character is ch
      Tuples0 = []
      dimer   = peptide1[k-1]+peptide2[k-1] #peptides are 0-indexed
      for tuple in D[dimer]:
        if tuple[-1]==ch: Tuples0.append(tuple)
      #-----------------------------------------------
      #remove tuples that don't satisfy constraints
      #-----------------------------------------------
      Tuples = []
      if constraint==None: 
        Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
      else:
        if len(constraint)<=3*(k-1):
          constraint4tuple = ''
        else:
          constraint4tuple = constraint[3*(k-1):3*k+1]
        for tuple in Tuples0:
          ok = 1
          for i in range(len(constraint4tuple)):
            ch = constraint4tuple[i]
            if ch!='N' and tuple[i]!=ch:
              ok = 0
              break
          if ok: 
            Tuples.append(tuple)
      #-----------------------------------------------
      #choose random tuple with first char ch1 for which ZF[(k-1,ch1)]>0
      #-----------------------------------------------
      tuple = random.choice(Tuples)
      while ZF[(k-1,tuple[0])]==0: tuple = random.choice(Tuples)
#      print (k-1),tuple,tuple[0],ZF[(k-1,tuple[0])],"+++"
      rna = tuple[1]+tuple[2]+tuple[3]+rna
      ch  = tuple[0] #set ch for next round
      k   = k-1
    rna = ch+rna #append the first character, which remains to be appended
    SAMPLES.append(rna)
  return SAMPLES
      


def printFASTA(RNA):
  for k in range(len(RNA)):
    print "> %s" % (k+1)
    print RNA[k]
 
def main(peptide1,peptide2,numSamples,filename=None,threshold=1,constraint=None):
  D      = createListOf4mersForAllPairsOfResidues(filename,threshold)
  seqLen = min(len(peptide1),len(peptide2))
  ZF     = computeZF(peptide1,peptide2,D,constraint)
  if DEBUG:
    for ch in NUCL:
      sys.stdout.write( ch+'\n' )
      for k in range(3):
        sys.stdout.write( "\tZF[(%s,%s)] = %s\n  " % (k,ch,ZF[(k,ch)]) )
  ZFval  = 0.0
  for ch in NUCL: ZFval += ZF[(seqLen,ch)]
  RNA    = sample(peptide1,peptide2,seqLen,numSamples,D,constraint)
  if PRINT:
    print "%f" % ZFval
    printRNA(RNA)
  return ZFval,RNA


if __name__ == '__main__':
  if len(sys.argv) < 3:
    text = """
Usage: %s peptide1 peptide2 [numSamples [constraint[[Blosum62matrix [threshold]]]]
 1) peptides given in IUPAC single letter codes and must have same length
 2) numSamples is desired number of samples to generate from PSSM defined
    from RNAs a_{0}...a_{3n} which translate peptide1 in reading frame 0
    and peptide2 in reading frame 1 (default 10)
 3) filename of amino acid (Blosum62) similarity matrix (default None)
 4) threshold in Blosum62 for codons that translate to 
    residues with >= threshold similarity (default +1)""" % sys.argv[0]
    print text
    text = """
This code can be used to do following:
a) Sample RNA sequences that translate into peptides that are either
   identical to, or Blosum similar to, peptide1 and peptide2. Such output
   is used with http://weblogo.berkeley.edu/ to generate a web logo.
b) In addition or instead of sampling, code can compute the PSSM, either
   by DP as in (a), or by generating a list of all RNA sequences that satisfy
   the given constraints. The latter is done by breadth first search, and
   may take exponential time, since the number of sequences grows exponentially.
   This portion of the code was used as well to test (a). Note that the
   DP algorithm takes time/space that is linear in the peptide length.
   Program can be extended (future work) to allow user to stipulate for EACH
   residue position whether an exact tranlated residue is necessary or a
   position specific Blosum-similar one.
Goal of program is to determine, in conjunction with RNAiFold, the extent
to which evolutionary pressure constrains the amino acids translated versus
the RNA secondary structure -- in HIV-1 there is a stem-loop programmed
-1 ribosomal frameshift signal' at the gag/pol framshift site, where pol
starts at position 1631 in GenBank AF033819.3. See also Figure 1 of
    Ofori, L.O., Hilimire, T.A., Bennett, R.P. et al.
    J. Med. Chem. 57(3), 723--732, February 2014"""
    print text
    sys.exit(1)
  peptide1 = sys.argv[1]
  peptide2 = sys.argv[2]
  #default settings that maybe changed by flags below
  numSamples = 10
  filename   = None
  threshold  = 1
  constraint = None
  args     = sys.argv[3:]
  #sequence constraint indicated by -c flag
  for i in range(len(args)):
    if args[i]=='-c':
      constraint = args[i+1]
      del[args[i+1]]; del[args[i]]
      break
    else:
      constraint = None
  #number of samples, with -n flag
  for i in range(len(args)):
    if args[i]=='-n':
      numSamples = int(args[i+1])
      del[args[i+1]]; del[args[i]]
      break
    else:
      numSamples = 10
  #similarity matrix filename, with -m flag
  for i in range(len(args)):
    if args[i]=='-m':
      filename = args[i+1]
      del[args[i+1]]; del[args[i]]
      break
    else:
      filename = None
  #blosum threshold, with -t flag
  for i in range(len(args)):
    if args[i]=='-t':
      threshold = int(args[i+1])
      del[args[i+1]]; del[args[i]]
      break
    else:
      threshold = 1
  main(peptide1,peptide2,numSamples,filename,threshold,constraint)

