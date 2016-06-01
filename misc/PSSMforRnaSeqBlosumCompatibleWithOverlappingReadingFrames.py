#! /usr/bin/env python

# PSSMforRnaSeqBlosumCompatibleWithOverlappingReadingFrames.py
# P.Clote

# Given peptide1,peptide2 this program returns all mRNA that APPROXIMATELY translate peptide1 (reading frame 0)
# and APPROXIMATELY translate peptide2 (reading frame +1). In contrast to the related program
#        rnaSeqCompatibleWithOverlappingReadingFrames.py
# this program returns mRNAs that translate into peptides, such that each amino acid of the peptide has Blosum62
# similarity at least the user-specified threshold with the input peptide1 [resp. peptide2]



import sys,os,copy,random
from aminoAcidAndGeneticCodes import genCode,aaCode, singleLetterAAcodes
from extractHIVgagPolJunction import translateMRNA

PRINT = False
NUCL  = ['A','C','G','U']


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

def computeZF(peptide1,peptide2,D): #forwards "partition function"
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
	for k in range(1,n+1):
		for ch in NUCL:
			ZF[(k,ch)] = 0.0
	for ch in NUCL: ZF[(0,ch)] = 1.0
	#------------Fill matrix ZF-------------------
	for k in range(1,n+1): #k in [1,...,n] 
		dimer = p[k-1]+q[k-1]
		L     = D[dimer]
		for s in L:
			ZF[(k,s[-1])] += ZF[(k-1,s[0])]
	ZFval=0.0
	for ch in NUCL: ZFval += ZF[(n,ch)]			
	print '%f'%ZFval
	return ZF

def computeZB(peptide1,peptide2,D): #backwards "partition function"
	# Here, peptide1 = p_{0}...p_{n-1} and peptide2 = q_{0}...q_{n-1}
	# Function computes the ZB[(k,ch)] = number of RNA sequences
	#    s = a_{3k}...a_{3n} where len(s)=3*(n-k)+1  and such that
	# (1) translateMRNA(s) = peptide1[k:] 
	# (2) translateMRNA(s[1:]) = peptide2[k:]
	# (3) s[0] = ch
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
	#------------Initialize ZB---------------------
	ZB = {}
	for k in range(n):
		for ch in NUCL:
			ZB[(k,ch)] = 0.0
	for ch in NUCL: ZB[(n,ch)] = 1.0
	#------------Fill matrix ZB-------------------
	for k in range(n-1,-1,-1): #k in [n-1,n-2,...,0] 
		dimer = p[k]+q[k]
		L     = D[dimer]
		for s in L:
			ZB[(k,s[0])] += ZB[(k+1,s[-1])]
	ZFval=0.0
	for ch in NUCL: ZFval += ZB[(0,ch)]			
	print '%f'%ZFval
	return ZB

def createPSSM(peptide1,peptide2,D,ZF,ZB):
	# Here, peptide1 = p_{0}...p_{n-1} and peptide2 = q_{0}...q_{n-1}
	# Function computes PSSM[(k,ch) = frequency of ch in position k of the
	# multiple alignment of ALL sequences
	#    s = a_{0}...a_{3n}  such that
	# translateMRNA(s) = peptide1, and translateMRNA(s[1:]) = peptide2
	#
	# ZF[(k,ch)] is the number of RNA sequences
	#    s = a_{0}...a_{3k} such that
	# (1) translateMRNA(s) = peptide1[:k] 
	# (2) translateMRNA(s[1:]) = peptide2[:k]
	# (3) s[-1] = ch 
	#
	# ZB[(k,ch)] is the number of RNA sequences
	#    s = a_{3k}...a_{3n}such that
	# (1) translateMRNA(s) = peptide1[k:] 
	# (2) translateMRNA(s[1:]) = peptide2[k:]
	# (3) s[0] = ch
	#
	# PSSM[(3*k,ch)] = ZF[(k,ch)]*ZB[(k,ch)]
	# by using ZF and ZB 
	#
	p = peptide1; q = peptide2; n = min(len(p),len(q))
	assert (len(p) <= len(q)) 
	denom1 = 0.0; denom2 = 0.0
	for ch in NUCL:
		denom1 += ZF[(n,ch)]
		denom2 += ZB[(0,ch)]
	assert (denom1 == denom2)
	#----------Initialize PSSM ---------------
	PSSM = {}
	for i in range(3*n+1):
		for ch in NUCL:
			PSSM[(i,ch)] = 0.0
	#----------Compute PSSM(3*k+i,ch) for i=0,1,2 ---------------
	for k in range(n+1):
		for ch in NUCL:
			PSSM[(3*k,ch)] =  ZF[(k,ch)]*ZB[(k,ch)]
	for k in range(n):
		dimer = p[k]+q[k]
		L     = D[dimer]
		for s in L:
			PSSM[(3*k+1,s[1])] += ZF[(k,s[0])]*ZB[(k+1,s[-1])]
			PSSM[(3*k+2,s[2])] += ZF[(k,s[0])]*ZB[(k+1,s[-1])]
	#----------Normalize PSSM ---------------
	for i in range(3*n+1): #i in [0,...,3n]
		for ch in NUCL:
			PSSM[(i,ch)] = PSSM[(i,ch)]/denom1
	return PSSM

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

def checkPSSM(PSSM):
	keys = PSSM.keys()
	keys.sort() # keys is list [(0,'A'),...,(3n,'U')] 
							# where n = min(len(peptide1),len(peptide2))
	N = len(keys)/4
	for i in range(N):
		As = PSSM[(i,'A')]
		Cs = PSSM[(i,'C')]
		Gs = PSSM[(i,'G')]
		Us = PSSM[(i,'U')]
		print "%d:\t%.4f\t%.4f\t%.4f\t%.4f" % (i,As,Cs,Gs,Us)
		assert (As+Cs+Gs+Us == 1)


def computePSSMfromList(RNA):
	if RNA==[]:
		return None
	#------------Initialize PSSM ------------------------------------
	n    = len(RNA[0])
	PSSM = {}
	for i in range(n):
		for ch in NUCL:
			PSSM[(i,ch)] = 0.0
	#------------Compute frequencies in RNA for PSSM  --------------
	for rna in RNA:
		for i in range(n):
			PSSM[(i,rna[i])] += 1
	#------------Normalize PSSM  ----------------------------------
	N = len(RNA)
	for i in range(n):
		for ch in NUCL:
			PSSM[(i,ch)] = PSSM[(i,ch)]/N
	return PSSM


def crossCheckPSSMvalues(PSSM,PSSMbis):
	assert (len(PSSM)==len(PSSMbis) )
	keys    = PSSM.keys();
	keysBis = PSSMbis.keys();
	keys.sort(); keysBis.sort()
	assert ( keys==keysBis )
	for k in keys:
		assert PSSM[k]==PSSMbis[k] 


def sample(seqLen,PSSM,numSamples):
	for i in range(numSamples):
		L = []
		for j in range(seqLen):
			z = random.random()
			if z<PSSM[(j,'A')]:
				L.append('A')
			elif z<PSSM[(j,'A')]+PSSM[(j,'C')]:
				L.append('C')
			elif z<PSSM[(j,'A')]+PSSM[(j,'C')]+PSSM[(j,'G')]:
				L.append('G')
			else:
				L.append('U')
		sampledString = "".join(L)
		return sampledString

def printFASTA(RNA):
	for k in range(len(RNA)):
		print "> %s" % (k+1)
		print RNA[k]
	

def main(peptide1,peptide2,numSamples,filename=None,threshold=1):
	D      = createListOf4mersForAllPairsOfResidues(filename,threshold)
	ZF     = computeZF(peptide1,peptide2,D)
	ZB     = computeZB(peptide1,peptide2,D)
	PSSM   = createPSSM(peptide1,peptide2,D,ZF,ZB)
	n      = min(len(peptide1),len(peptide2))
	RNA    = []
	for i in range(numSamples):
		RNA.append(sample(3*n+1,PSSM,numSamples))
	printRNA(RNA)
	if PRINT:
		RNA    = bfs(peptide1,peptide2,filename,threshold)
		PSSMbis= computePSSMfromList(RNA)
		crossCheckPSSMvalues(PSSM,PSSMbis)
		checkPSSM(PSSM)
		print "Num seq: %s\tLen seq: %s" % (len(RNA),len(RNA[0]))
		printRNA(RNA)
		print "Check translation of gag and pol"
		for i in range(1,len(RNA)+1):
			rna = RNA[i-1]
			print "%s\t%s" % (translateMRNA(rna),translateMRNA(rna[1:]))
		print check(peptide1,peptide2,RNA)


if __name__ == '__main__':
	if len(sys.argv) < 3:
		text = """
Usage: %s peptide1 peptide2 [numSamples [Blosum62matrix [threshold]]]
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
	if len(sys.argv)>3:
		numSamples = int(sys.argv[3])
	else:
		numSamples = 10
	if len(sys.argv)>4:
		filename = sys.argv[4]
	else:
		filename = None
	if len(sys.argv)>5:
		threshold = int(sys.argv[5])
	else:
		threshold = 1
	main(peptide1,peptide2,numSamples,filename,threshold)


