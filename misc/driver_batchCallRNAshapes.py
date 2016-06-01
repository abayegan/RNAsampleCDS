#! /usr/bin/env python 

#driver_batchCallRNAshapes.py
#P.Clote

#Program computes a number of measures including mean optimization,
#Z-score of probability to fold into stem-loop structure []
#for rna coding BLOSUM62 +1 similar peptides p,q to those coded by Rfam RNA 
#
#WARNING: It took me a long time to realize that the SlipperySeq should be
#defined to be the 7-mer prefix of an input RNA from RF00480, rather than
#the sequence from Ofori


import sys,os,random
from stats import getSampleStats
import sampleRNAsWithSequenceConstraintsWithOverlappingReadingFrames
from altschulEriksonDinuclShuffle import dinuclShuffle

RANDOM_RNA         = False
#NUMSAMPLES         = 100000
NUMSAMPLES         = 50
BLOSUM_THRESHOLD   = 1
SHAPES             = ['[]']
PRINT              = False

#SlipperySeq        = 'UUUUUUA'
#the slippery sequence is now a variable with value the 7-nt
#prefix of an input HIV FSS

#During debugging phase, can fix random seed
#random.seed(321)


def generateRandomRNAs(n,length,SlipperySeq):
  #generates n many RNAs each of size 'length' 
  RNAs = []
  for i in range(n):
    #test if initial sequence constraint present
    if SlipperySequence == '':
      L = []
    else:
      L  = [SlipperySeq] #slippery sequence
    #process
    for k in range(length-len(SlipperySeq)):
      x = random.random()
      if x<0.25:
        L.append('A')
      elif x<0.5:
        L.append('C')
      elif x<0.75:
        L.append('G')
      else:
        L.append('U')
    s = "".join(L)
    RNAs.append(s)
  return RNAs



def aux0(filename):
  #Parse a file containing two 17-mers peptides in each line for Pol resp Gag 
  #peptides in HIV-1 Gag-Pol overlap region and return list of peptides p,q
  L    = []
  file = open(filename)
  line = file.readline()
  while line:
    words = line.split()
    FASTA = words[0]
    p     = words[1]
    q     = words[2]
    rna   = words[3]
    assert( len(p)==len(q) )
    L.append( (FASTA,p,q,rna) )
    line = file.readline()
  file.close()
  return L

def runAndParseRNAshapesForListOfRNAs(RNAs):
  n   = 1  #used to create FASTA comment in temporary file
  tmp = open('./tmp','w')
  for rna in RNAs:
    tmp.write('> %s\n%s\n' % (n,rna))
    n += 1
  tmp.close()
  cmd = "RNAshapes -q -m '%s' < ./tmp" % '[]'
  cmd = 'RNAshapes -q < ./tmp'
  file = os.popen(cmd)
#  cmd = "RNAshapes -q -m '%s' < ./tmp > ./tmpOut" % '[]'
#  os.system(cmd)
#  file = open('./tmpOut')
  line = file.readline()
  L    = [] #L is list of (prob,energy,rna,shrep) 
  while line:
    if line[0]=='>':
      FASTA = line[1:].strip() #FASTA comment
      rna   = file.readline().strip() #RNA sequence
      words = file.readline().split() #data from RNAshapes -q
      if words[-1] in SHAPES:
        prob   = float(words[-2])
        energy = float(words[0])
        shrep  = words[1]
        L.append( (prob,energy,rna,shrep) )
        line = file.readline()
      else:
        #Line is "shape [] not found."
        line = file.readline()
        line = file.readline()
      #WARNING: this part of code needs reworking if SHAPES contains more
      #than one shape. Code works correctly if SHAPES contains only one shape
      #This is the case for application with shape []
    else:
      line = file.readline()
  file.close()
  Probs = []; ShrepEnergies = []
  for word in L:
    prob,energy,rna,shrep = word
    Probs.append(prob); ShrepEnergies.append(energy)
  probMean,probStdev,max,min = getSampleStats(Probs)
  shrepEnergyMean,shrepEnergyStdev,max,min = getSampleStats(ShrepEnergies)
  if PRINT:
    for word in L:
      prob,energy,rna,shrep = word
      print "%s\t%s\t%s\t%s" % (prob,energy,rna,shrep)
  return probMean,probStdev, shrepEnergyMean,shrepEnergyStdev


def runAndParseRNAshapes(rna):
  #This function makes individual call on RNAshapes rather than a batch call
  #Using this function is SLOW! and better to use
  #     function runAndParseRNAshapesForListOfRNAs(RNAs)
  #This outputs probSum and avgEnergy for structures of rna having
  #a shape in SHAPES
  #
  #However, unlike the previous function, this function works correctly
  #if SHAPES contains more than one shape
  probSum = 0; avgEnergy = 0
  cmd  = "echo %s | RNAshapes -q" % rna
  L    = [] #L is list of (prob,energy,rna,shrep) for the target shape(s)
  file = os.popen(cmd)
  line = file.readline()                          #skip first ine
  assert( line.strip()==rna ),"%s,%s"%(rna,line) #first line is rna
  ok   = False                                   #found target SHAPE ?
  line = file.readline()                         #data from RNAshapes -q
  while line:
    words = line.split()
    if words[-1] in SHAPES:
      #found target shape
      ok     = True
      prob   = float(words[-2])
      energy = float(words[0])
      shrep  = words[1]
      L.append( (prob,energy,rna,shrep) )
    line  = file.readline()
  file.close()
  if ok: #Note: size of L bounded by num shapes in SHAPES
    energySum = 0.0
    for i in range(len(L)):
      prob,energy,rna,shrep = L[i]
      probSum   += prob
      energySum += prob*energy  #weighted energy
    avgEnergy = energySum/probSum     #normalize by sum of probabilities
  return probSum,avgEnergy

def aux1(rna,p,q):
  #INPUT : rna that codes p,q in frames 0,1
  #OUTPUT: 
  #Compute Z-score of prob that rna folds into target shape (eg. [])
  #and Z-score of shrep free energies. Note that several shapes can
  #be considered, such as _[[]_] and _[_[]_] with RNAshapes -t 1
  #----------------------------------------------------------------
  #Compute probability of stem-loop using RNAshapes for input rna
  rna0          = rna
  prob0,energy0 = runAndParseRNAshapes(rna)
  #----------------------------------------------------------------
  #Compute probability of stem-loop using RNAshapes for each rna in L
  Probs = []; ShrepEnergies = []; Optimization = []; L = []
  if RANDOM_RNA:
    ZFval = (4.0)**len(rna)
    RNAs  = generateRandomRNAs(NUMSAMPLES,len(rna),rna[:7])
            #initial 7-nt prefix of rna is the slippery sequence
  else:
    #if RANDOM_CODING_RNA then set BLOSUM = -100 
    #otherwise want BLOSUM=1 or some reasonable value
    SlipperySeq= rna0[:7] #initial 7-nt prefix is slippery sequence
    ZFval,RNAs = sampleRNAsWithSequenceConstraintsWithOverlappingReadingFrames.main(p,q,NUMSAMPLES,'blosum62.txt',BLOSUM_THRESHOLD,SlipperySeq) 
  a,b,c,d = runAndParseRNAshapesForListOfRNAs(RNAs)
  probMean         = a
  probStdev        = b
  shrepEnergyMean  = c
  shrepEnergyStdev = d
  return probMean,probStdev,shrepEnergyMean,shrepEnergyStdev,ZFval,prob0,energy0

def mfe(rna):
  cmd    = "echo %s | RNAfold " % rna
  file   = os.popen(cmd)
  line   = file.readline() #rna
  energy = file.readline().split()[-1].replace('(','').replace(')','')
  return float(energy)

def mfeAndProb(rna): #function not used in sequel
  cmd    = "echo %s | RNAshapes -q -m '[]' " % rna
  file   = os.popen(cmd)
  line   = file.readline() #rna, i.e.  GGGGGGCCCCCCCCUUUUUUUAAAAAAA
  line   = file.readline() #-11.70  (((((....)))))..............  0.8528473  []
  words  = line.split()
  energy = float(words[0])
  prob   = float(words[2])
  return energy,prob

def main0(filename):
  #used to determine prob of [] and shrep energy for all RNAs that 
  #code in frames 0,1 regardless of peptide coded, or just of random RNA
  #L is a list of FASTA, p, q, and rna seq for the Gag-Pol RNA from HIV-1
  print "FASTA \t optMean \t optStdev \t probMean \t probStdev \
\t shrepEnergyMean \t shrepEnergyStdev \t ZFval \t prob0 \
\t energy0 \t zScoreEnergy \t zScoreProb"
  L    = aux0(filename)
  FASTA,p,q,rna = L[0]
  a1,a2,a3,a4,a5,a6,a7 = aux1(rna,p,q)
  #rna0 is the current rna from list L
  probMean   = a1 #avg value of P(rna has shape [])
  probStdev  = a2 #stdev of P(rna has shape [])
  shrepEnergyMean   = a3 #avg value of energy of shrep for []
  shrepEnergyStdev  = a4 #stdev of energy of shrep for []
  ZFval             = a5 #number of rnas that code p,q in frames 0,1
  prob0             = a6 #prob that rna0 has shape []
  energy0           = a7 #energy of shrep for [] for rna0
  #energy0,prob0    = mfeAndProb(rna) #already computed
  zScoreEnergy      = (energy0 - shrepEnergyMean)/shrepEnergyStdev
  zScoreProb        = (prob0 - probMean)/probStdev
  #print FASTA,a1,a2,a3,a4,a5,a6,a7,zScoreEnergy,zScoreProb 
  sys.stdout.write("%s %s %s %s %s %s %s %s %s %s\n" % \
      (FASTA,a1,a2,a3,a4,a5,a6,a7,zScoreEnergy,zScoreProb))
  sys.stdout.flush()

def main(filename):
  print "FASTA \t optMean \t optStdev \t probMean \t probStdev \
\t shrepEnergyMean \t shrepEnergyStdev \t ZFval \t prob0 \
\t energy0 \t zScoreEnergy \t zScoreProb"
  L    = aux0(filename)
  for FASTA,p,q,rna in L:
    a1,a2,a3,a4,a5,a6,a7 = aux1(rna,p,q)
    #rna0 is the current rna from list L
    probMean   = a1 #avg value of P(rna has shape [])
    probStdev  = a2 #stdev of P(rna has shape [])
    shrepEnergyMean   = a3 #avg value of energy of shrep for []
    shrepEnergyStdev  = a4 #stdev of energy of shrep for []
    ZFval             = a5 #number of rnas that code p,q in frames 0,1
    prob0             = a6 #prob that rna0 has shape []
    energy0           = a7 #energy of shrep for [] for rna0
    #energy0,prob0    = mfeAndProb(rna) #already computed
    zScoreEnergy      = (energy0 - shrepEnergyMean)/shrepEnergyStdev
    zScoreProb        = (prob0 - probMean)/probStdev
    #print FASTA,a1,a2,a3,a4,a5,a6,a7,zScoreEnergy,zScoreProb 
    sys.stdout.write("%s %s %s %s %s %s %s %s %s %s\n" % \
        (FASTA,a1,a2,a3,a4,a5,a6,a7,zScoreEnergy,zScoreProb))
    sys.stdout.flush()


if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s filename" % sys.argv[0]
    print "  Each line of file contains 4 items separated by spaces or a tab:"
    print "  1) FASTA comment (one word, no spaces, with/without '>' symbol)"
    print "  2) 17-mer Pol peptide from reading frame 0"
    print "  3) 17-mer Gag peptide from reading frame 1"
    print "  4) 52-nt RNA which translates into Pol [resp Gag] 17-mer"
    print "     in reading frame 0 [resp 1]."
    print "WARNING: Code makes popen call of RNAshapes -q -m '[]' "
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)


