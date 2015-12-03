#! /usr/bin/python2.6

# overcode.py
# P.Clote - A. Bayegan

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
PRINT  = True #set False if use main as function to return list


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

def nuc2int(nuc):
	if nuc == 'A': return 0
	elif nuc == 'C': return 1
	elif nuc == 'G': return 2
	elif nuc == 'U' or nuc == 'T' : return 3
	else:
		print "error in nuc2int! Sequences must only contain A,C,G,U,T!"
		sys.exit(1)
	
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

def GCcont(seq):
	sum=0
	for ch in seq:
		if ch=='G' or ch=='C':
			sum += 1
	return sum

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


def computeZF0(n,constraint):
	#Compute number of RNAs of len 3*n+1 that code in frames 0,1 and
	#respect constraint
	#
	L = createListOfOverlappingCodons()
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
			for s in L:
				ZF[(k,s[-1])] += ZF[(k-1,s[0])]
	else: #must satisfy sequence constraint
		for k in range(1,n+1): #k in [1,...,n] 
			if len(constraint)<=3*(k-1):
				constraint4tuple = ''
			else:
				constraint4tuple = constraint[3*(k-1):3*k+1]
			L0                = []
			for tuple in L:
				ok = 1
				for i in range(len(constraint4tuple)):
					ch = constraint4tuple[i]
					if ch!='N' and tuple[i]!=ch:
						ok = 0
						#break
				if ok: 
					L0.append(tuple)
			if DEBUG:
				print "computeZF ",k,L0
			for s in L0:
				ZF[(k,s[-1])] += ZF[(k-1,s[0])]
	return ZF

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

def computeZF_GC(peptide1,peptide2,D,constraint=None): 
	#forwards "partition function"
	# Here, peptide1 = p_{0}...p_{n-1} and peptide2 = q_{0}...q_{n-1}
	# Function computes the ZF[(k,x,ch)] = number of RNA sequences
	#    s = a_{0}...a_{3k} where len(s)=3*k+1  and such that
	# (1) translateMRNA(s) = peptide1[:k] 
	# (2) translateMRNA(s[1:]) = peptide2[:k]
	# (3) s[-1] = ch
	# (4) x number of G's and C's in s
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
	#ZF1={}
	for k in range(0,n+1):
		for ch in NUCL:
			#ZF1[(k,ch)]=0.0
			for x in range(0,3*k+6): #Warning: 3k+2 doesn't work.
				ZF[(k,x,ch)] = 0.0
	if constraint==None or constraint[0]=='N':
			dimer = p[0]+q[0]
			L     = D[dimer]
			#~ print "tuples at level 1:", L
			for s in L:
				ZF[(1,GCcont(s),s[-1])] += 1.0
				#ZF1[(1,s[-1])] +=1.0
			if(DEBUG):
				for a in NUCL:
					for x in range(0,5):
						print 'ZF[1,%d,%s]:%d' % (x,a,ZF[(1,x,a)])
	else: #first position is constrained     
		ch = constraint[0]
		ZF[(0,ch)] = 1.0
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(2,n+1):
			dimer = p[k-1]+q[k-1]
			L     = D[dimer]
			if(DEBUG):
				print "tuples at level %d:" %k , L
			for s in L:
				#ZF1[(k,s[-1])]+=ZF1[(k-1,s[0])]
				for x in range(0,3*k+2):
					#print s,x,GCcont(s[1:4])
					if(x-GCcont(s[1:])>=0):
						ZF[(k,x,s[-1])] += ZF[(k-1,x-GCcont(s[1:]),s[0])]
			if(DEBUG):
				for a in NUCL:
					for x in range(0,3*k+2):
						print 'ZF[%d,%d,%s]:%d' % (k,x,a,ZF[(k,x,a)])
		#ZFval  = 0.0
		#for ch in NUCL: ZFval += ZF1[(n,ch)]
		#print "%f" %ZFval
		
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
	#~ print "CG-ZF"
	#~ for k in sorted(ZF):
		#~ print k,ZF[k]
	return ZF

def printRNA(RNA):
	for rna in RNA:
		print rna

def getNucSampleProb(ZF,k,nuclist):
	probdict={}
	sum=0.0
	for nuc in nuclist:
		sum += ZF[k,nuc]
	for nuc in nuclist:
		probdict[nuc]=(ZF[k,nuc]/sum)
	return probdict

def sample(peptide1,peptide2,seqLen,numSamples,D,constraint,weighted):
	#D[4-tuple] is list of "similar" 4-tuplesaa1+aa2]
	#Here, the key 4-tuple codes aa1 in frame 0 and aa2 in frame 1
	#4-tuples in list L=D[4-tuple] code aa1' in frame 0 and aa2' in frame 1
	#where Blosum(aa1,aa1')>= threshold value, and Blosum(aa2,aa2')>= value.
	#WARNING: Do NOT need to worry about constraints here, since ZF handles this
	assert( len(peptide1)==len(peptide2) )
	SAMPLES = []
	ZF      = computeZF(peptide1,peptide2,D,constraint)
	#GCZF    =  computeZF_GC(peptide1,peptide2,D,constraint)
	
	for i in range(numSamples):
		#initialize empty sequence and set k=n
		k  = seqLen; rna = ""
		#--------choose random nucleotide for which ZF[(k,ch)]>0 
		if weighted==1:
			sampleDict = getNucSampleProb(ZF,k,NUCL)
			#print sampleDict
			ch4 = rouletteWheel(sampleDict)
			while ZF[(k,ch4)]==0: ch4 = rouletteWheel(sampleDict)
			#print 'selected ch4',ch4	
		else:
			ch4 = random.choice(NUCL)
			while ZF[(k,ch4)]==0: ch4 = random.choice(NUCL)
		#--------construct sample
		while k>0:
			#create list of 4-tuples having similarity to p[k],q[k] in frame 0,1
			#and whose last character is ch
			Tuples0 = []
			dimer   = peptide1[k-1]+peptide2[k-1] #peptides are 0-indexed
			for tuple in D[dimer]:
				if tuple[-1]==ch4: Tuples0.append(tuple)
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
			print 'tuples before filteration:'
			print Tuples
			if weighted==1:
				selTuples=[];nuclist=[]
				for tup in Tuples: #get the set of nucleotides in t[0]
						if nuclist==NUCL: break
						if tup[0] not in nuclist: nuclist.append(tup[0])
				print 'nuclist', nuclist
				sampleDict = getNucSampleProb(ZF,k-1,nuclist)
				print sampleDict
				ch1 = rouletteWheel(sampleDict)
				while ZF[(k-1,ch1)]==0: ch1 = rouletteWheel(sampleDict)
				print 'sampled ch1',ch1
				for tup in Tuples:
					#print tup
					if(tup[0]==ch1 and k!=1): 
						selTuples.append(tup)
					if(k==1): #Warning: since ZF[0,nuc]=1 (nuc in {A,C,G,U}) the first nucleotide does not need to be ch1 and can be anything 
						selTuples = Tuples
				print 'tuples after filteration:'
				print selTuples
				tuple = random.choice(selTuples)
				print "selected tuple:",tuple
			else:
				#print Tuples
				tuple = random.choice(Tuples)
				while ZF[(k-1,tuple[0])]==0: tuple = random.choice(Tuples)
				
				#print (k-1),tuple,tuple[0],ZF[(k-1,tuple[0])],"+++"
			#ZF[(k,tuple[3])] -= 1.0
			print "-1:",k,tuple[3]
			rna = tuple[1]+tuple[2]+tuple[3]+rna
			print rna
			ch4  = tuple[0] #set ch for next round
			k   = k-1
		rna = ch4+rna #append the first character, which remains to be appended
		for k in sorted(ZF):
			print k,ZF[k]
		SAMPLES.append(rna)
	return SAMPLES
			
def sample_GC(peptide1,peptide2,seqLen,numSamples,GC,D,constraint,weighted):
	#D[4-tuple] is list of "similar" 4-tuplesaa1+aa2]
	#Here, the key 4-tuple codes aa1 in frame 0 and aa2 in frame 1
	#4-tuples in list L=D[4-tuple] code aa1' in frame 0 and aa2' in frame 1
	#where Blosum(aa1,aa1')>= threshold value, and Blosum(aa2,aa2')>= value.
	#WARNING: Do NOT need to worry about constraints here, since ZF handles this
	assert( len(peptide1)==len(peptide2) )
	SAMPLES = []
	ZF      = computeZF_GC(peptide1,peptide2,D,constraint)
	for i in range(numSamples):
		#initialize empty sequence and set k=n
		k  = seqLen; rna = ""; possible=-1; gc=GC
		#--------choose random nucleotide for which ZF[(k,ch)]>0 
		random.shuffle(NUCL)
		print NUCL
		for ch4 in NUCL:
			if(ZF[k,gc,ch4]!=0):
				possible =1
				break
		if(possible==-1):
			print "no more sequences with GC-content %d found!" %gc
			return SAMPLES
		#gc -= GCcont(ch4)
		#--------construct sample
		while k>0:
			#create list of 4-tuples having similarity to p[k],q[k] in frame 0,1
			#and whose last character is ch
			Tuples0 = []
			dimer   = peptide1[k-1]+peptide2[k-1] #peptides are 0-indexed
			for tuple in D[dimer]:
				if tuple[-1]==ch4: Tuples0.append(tuple)
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
			print 'tuples before filteration:'
			print Tuples
			random.shuffle(Tuples)
			for tuple in Tuples:
				tupleGC = GCcont(tuple[1:])
				if((gc-tupleGC)>=0):
					if(ZF[k-1,gc-tupleGC,tuple[0]]!=0):
						break
			gc -= GCcont(tuple[1:])
			print (k-1),tuple,tuple[0],"+++"
			#ZF[(k,tuple[3])] -= 1.0
			print "-1:",k,tuple[3]
			rna = tuple[1]+tuple[2]+tuple[3]+rna
			print rna
			ch4  = tuple[0] #set ch for next round
			k   = k-1
		rna = ch4+rna #append the first character, which remains to be appended
		#for k in sorted(ZF):
			#print k,ZF[k]
		SAMPLES.append(rna)
	return SAMPLES

def printFASTA(RNA):
	for k in range(len(RNA)):
		print "> %s" % (k+1)
		print RNA[k]

def main0(n,constraint=None):
	ZF     = computeZF0(n,constraint)
	if DEBUG:
		for ch in NUCL:
			sys.stdout.write( ch+'\n' )
			for k in range(3):
				sys.stdout.write( "\tZF[(%s,%s)] = %s\n  " % (k,ch,ZF[(k,ch)]) )
	ZFval  = 0.0; seqLen = n
	for ch in NUCL: ZFval += ZF[(seqLen,ch)]
	print "%f" % ZFval
	return ZFval

def main(peptide1,peptide2,numSamples,weighted,gc,filename=None,threshold=1,constraint=None):
	D      = createListOf4mersForAllPairsOfResidues(filename,threshold)
	seqLen = min(len(peptide1),len(peptide2))
	if (gc==-1):
		ZF     = computeZF(peptide1,peptide2,D,constraint)
		if DEBUG:
			for ch in NUCL:
				sys.stdout.write( ch+'\n' )
				for k in range(3):
					sys.stdout.write( "\tZF[(%s,%s)] = %s\n  " % (k,ch,ZF[(k,ch)]) )
		ZFval  = 0.0
		for ch in NUCL: ZFval += ZF[(seqLen,ch)]
		RNA    = sample(peptide1,peptide2,seqLen,numSamples,D,constraint,weighted)
	else:
		ZF     = computeZF_GC(peptide1,peptide2,D,constraint)			
		ZFval  = 0.0
		for ch in NUCL:
			for x in range(3*seqLen+2):ZFval += ZF[(seqLen,x,ch)]	
		RNA    = sample_GC(peptide1,peptide2,seqLen,numSamples,gc,D,constraint,weighted)
	if PRINT:
		print "%f" % ZFval
		printRNA(RNA)
		for rna in (RNA):
			print GCcont(rna)
	return ZFval,RNA

def printusage(exe):
	text = """
Usage: %s peptide1 peptide2 -n numSamples -c constraint -m BlosumMatrix -t threshold -gc GC_content -w

1) peptides given in IUPAC single letter codes and must have same length
2) numSamples is desired number of samples to generate from PSSM defined
		from RNAs a_{0}...a_{3n} which translate peptide1 in reading frame 0
		and peptide2 in reading frame 1 (default 10)
3) filename of amino acid (Blosum62) similarity matrix (default None)
4) threshold in Blosum62 for codons that translate to 
		residues with >= threshold similarity (default +1)
5) GC-content of the sampled sequences
6) Representative sampling """ % exe
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

if __name__ == '__main__':
	if len(sys.argv) < 3:
		printusage(sys.argv[0])
	peptide1 = sys.argv[1]
	peptide2 = sys.argv[2]
	#default settings that maybe changed by flags below
	numSamples = 10
	filename   = None
	threshold  = 1
	weighted = -1
	gc=-1
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
	#GC content with -gc flag
	for i in range(len(args)):
		if args[i]=='-gc':
			gc = int(args[i+1])
			if(gc<0):
				print "\nGC-content should not be negative!"
				printusage(sys.argv[0])
			del[args[i+1]]; del[args[i]]
			break
	#weighted sampling, with -w flag
	for i in range(len(args)):
		if args[i]=='-w':
			weighted = 1
			break
	n   = len(peptide1)
#  main0(n,constraint)
	main(peptide1,peptide2,numSamples,weighted,gc,filename,threshold,constraint)

