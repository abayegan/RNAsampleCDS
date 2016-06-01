#! /usr/bin/python2.6

# P.Clote - A. Bayegan
"""This program samples RNA sequences that translate into peptide1 and peptide2 in reading frames 0 and 1 respectively.
Translation could be based on either exact or Blosum similar amino acids. In addition, RNA sequence constraints and weighted
sampling are supported. In weighted sampling nucleotide frequencies are representative of frequencies in the total partition function."""

import sys,os,copy,random
from math import log,exp
from utility import *

DEBUG  = False
VERBOSE= False
PRINT  = True #set False if use main as function to return list
	
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
						codon2 = codon1[1:][::-1]+d
						if codon2[::-1] not in STOP:
							overlappingCodons = a+b+c+d
							L.append(overlappingCodons)
	#print '# of overlapping codons',len(L)
	return L

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
				aa2bis = aaCode[genCode[tuple[1:][::-1]].upper()]
				if filename==None:
					if aa1==aa1bis and aa2==aa2bis:
						Tuples.append(tuple)
				else: #use Blosum similarity; recall if aa==aabis then Blosum similar
					if Blosum[aa1][aa1bis]>=threshold and  Blosum[aa2][aa2bis]>=threshold: Tuples.append(tuple)
			D[aa1+aa2] = Tuples
	return D

#~ def computeZB(peptide1,peptide2,D,constraint=None):
	#~ p = peptide1; q = peptide2
	#~ assert (len(p) <= len(q)) 
	#~ n  = min(len(p),len(q)) 
	#~ #------------Initialize ZF---------------------
	#~ ZB = {}
	#~ for k in range(0,n+1):
		#~ for ch in NUCL:
			#~ ZB[(k,ch)] = 0.0
	#~ for ch in NUCL: ZB[(n,ch)] = 1.0
	#~ #------------Fill matrix ZF-------------------
	#~ for k in range(n-1,-1,-1): #k in [1,...,n] 
		#~ dimer = p[k]+q[k]
		#~ L     = D[dimer]
		#~ print "tuples at level %s:"%k ,L
		#~ for s in L:
			#~ ZB[(k,s[0])] += ZB[(k+1,s[-1])]
	#~ ZFval  = 0.0
	#~ for ch in NUCL: ZFval += ZB[(0,ch)]			
	#~ #print '%f'%ZFval
	#~ for key in sorted(ZB):
		#~ print key,ZB[key]
	#~ return ZB
	
def computeZF(peptide1,peptide2,D,constraint=None): 
	#forwards "partition function"
	# Here, peptide1 = p_{0}...p_{n-1} and peptide2 = q_{0}...q_{n-1}
	# Function computes the ZF[(k,ch)] = number of RNA sequences
	#    s = a_{0}...a_{3k} where len(s)=3*k+1  and such that
	# (1) translateMRNA(s) = peptide1[:k] 
	# (2) translateMRNA(s[1:][::-1]) = peptide2[n-k:]
	# (3) s[-1] = ch 
	# The final coding sequences are of the form a_{0} ... a_{3n}.
	#
	# Parameter explanation:
	# D is dictionary of 4-mer RNAs corresponding to peptide1 and peptide2
	# D is returned by function
	#    createListOf4mersForAllPairsOfResidues(filename,threshold)
	#
	#WARNING: This assumes that peptide1 is in reading frame 0 and 
	#         peptide2 is in reading frame 1 of the negative strand.
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
		for ch in constraint[0]:
			ZF[(0,ch)]=1.0
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(1,n+1): #k in [1,...,n] 
			dimer = p[k-1]+q[n-k]
			#print dimer
			L     = D[dimer]
			#print L
			for s in L:
				ZF[(k,s[-1])] += ZF[(k-1,s[0])]
	else: #must satisfy sequence constraint
		for k in range(1,n+1): #k in [1,...,n] 
			dimer            = p[k-1]+q[n-k]
			L0               = D[dimer]
			constraint4tuple = constraint[3*(k-1):3*k+1]
			L                = []
			for tuple in L0:
				ok = 1
				for i in range(len(constraint4tuple)):
					ch = constraint4tuple[i]
					if tuple[i] not in iupacDict[ch]:
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
	DEBUG=0
	p = peptide1; q = peptide2
	assert (len(p) <= len(q)) 
	n  = min(len(p),len(q)) 
	#------------Initialize ZF---------------------
	ZF = {}
	for k in range(0,n+1):
		for ch in NUCL:
			for x in range(0,3*k+6): #Warning: 3k+2 doesn't work.
				ZF[(k,x,ch)] = 0.0
	if constraint==None or constraint[0]=='N':
		dimer = p[0]+q[n-1]
		L     = D[dimer]
		for s in L:
			ZF[(1,GCcont(s),s[-1])] += 1.0
		if(DEBUG==1):
			print "tuples at level 1:", L
			for a in NUCL:
				for x in range(0,5):
					print 'ZF[1,%d,%s]:%d' % (x,a,ZF[(1,x,a)])
	else: #first tuple is constrained
		dimer = p[0]+q[n-1]
		L0     = D[dimer]
		constraint4tuple = constraint[:4]
		L=[]
		for tuple in L0:
			ok = 1
			for i in range(len(constraint4tuple)):
				ch = constraint4tuple[i]
				if tuple[i] not in iupacDict[ch]:
					ok = 0
					#break
			if ok: 
				L.append(tuple)
		for s in L:
				ZF[(1,GCcont(s),s[-1])] += 1.0    
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(2,n+1):
			dimer = p[k-1]+q[n-k]
			L     = D[dimer]
			if(DEBUG==1):
				print "tuples at level %d:" %k , L
			for s in L:
				for x in range(0,3*k+2):
					if(x-GCcont(s[1:])>=0):
						ZF[(k,x,s[-1])] += ZF[(k-1,x-GCcont(s[1:]),s[0])]
			if(DEBUG==1):
				for a in NUCL:
					for x in range(0,3*k+2):
						print 'ZF[%d,%d,%s]:%d' % (k,x,a,ZF[(k,x,a)])
		#ZFval  = 0.0
		#for ch in NUCL: ZFval += ZF1[(n,ch)]
		#print "%f" %ZFval
		
	else: #must satisfy sequence constraint
		for k in range(1,n+1): #k in [1,...,n] 
			dimer            = p[k-1]+q[n-k]
			L0               = D[dimer]
			constraint4tuple = constraint[3*(k-1):3*k+1]
			L                = []
			for tuple in L0:
				ok = 1
				for i in range(len(constraint4tuple)):
					ch = constraint4tuple[i]
					if tuple[i] not in iupacDict[ch]:
						ok = 0
						#break
				if ok: 
					L.append(tuple)
			for s in L:
				for x in range(0,3*k+2):
					if(x-GCcont(s[1:])>=0):
						ZF[(k,x,s[-1])] += ZF[(k-1,x-GCcont(s[1:]),s[0])]
			if(DEBUG==1):
				for a in NUCL:
					for x in range(0,3*k+2):
						print 'ZF[%d,%d,%s]:%d' % (k,x,a,ZF[(k,x,a)])
	#print "CG-ZF"
	#for k in sorted(ZF):
		#print k,ZF[k]
	return ZF

def getNucSampleProb(ZF,k,nuclist):
	probdict={}
	sum=0.0
	for nuc in nuclist:
		sum += ZF[k,nuc]
	for nuc in nuclist:
		probdict[nuc]=(ZF[k,nuc]/sum)
	return probdict

def sample(peptide1,peptide2,seqLen,numSamples,replacement,D,constraint,weighted):
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
			possible=0
			shNUCL= NUCL
			random.shuffle(shNUCL)
			for ch4 in shNUCL:
				if ZF[(k,ch4)]>0:
					possible = 1
					break
			if not possible:
				print "no sequence coding the given peptides in 0 and 1 reading frames of the opposite strands found!"
				sys.exit(0)
		#--------construct sample
		while k>0:
			#create list of 4-tuples having similarity to p[k],q[k] in frame 0,1
			#and whose last character is ch
			Tuples0 = []
			dimer   = peptide1[k-1]+peptide2[n-k] #peptides are 0-indexed
			for tuple in D[dimer]:
				if tuple[-1]==ch4: Tuples0.append(tuple)
			#-----------------------------------------------
			#remove tuples that don't satisfy constraints
			#-----------------------------------------------
			Tuples = []
			if constraint==None: 
				Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
			else:
				constraint4tuple = constraint[3*(k-1):3*k+1]
				for tuple in Tuples0:
					ok = 1
					for i in range(len(constraint4tuple)):
						ch = constraint4tuple[i]
						if tuple[i] not in iupacDict[ch]:
							ok = 0
							break
					if ok:
						Tuples.append(tuple)
			#-----------------------------------------------
			#choose random tuple with first char ch1 for which ZF[(k-1,ch1)]>0
			#-----------------------------------------------
			if weighted==1:
				#sum = 0; #this toulette wheel does not work
				##print "k",k
				#r = (random.random())*ZF[k-1,ch4]
				##print ch4 , ZF[k-1,ch4] , r , Tuples
				#shTuples = Tuples
				#random.shuffle(shTuples)
				#for tuple in shTuples:
					##print tuple
					#sum += ZF[k-1,tuple[0]]
					##print sum	
					#if r < sum:
						##print 'selected for %d' %k
						#rna = tuple[1]+tuple[2]+tuple[3]+rna
						#ch4  = tuple[0] #set ch for next round
						#k   = k-1
						#break
						##print "not broke!!!!!"
				
				selTuples=[];nuclist=[]
				for tup in Tuples: #get the set of nucleotides in t[0]
						if nuclist==NUCL: break
						if tup[0] not in nuclist: nuclist.append(tup[0])
				#print 'nuclist', nuclist
				sampleDict = getNucSampleProb(ZF,k-1,nuclist)
				#print sampleDict
				ch1 = rouletteWheel(sampleDict)
				#print 'sampled ch1',ch1
				for tup in Tuples:
					#print tup
					if(tup[0]==ch1 and k!=1): 
						selTuples.append(tup)
					if(k==1): #Warning: since ZF[0,nuc]=1 (nuc in {A,C,G,U}) the first nucleotide does not need to be ch1 and can be anything 
						selTuples = Tuples
				#print 'tuples after filteration:'
				#print selTuples
				tuple = random.choice(selTuples)
				#print "selected tuple:",tuple
				rna = tuple[1]+tuple[2]+tuple[3]+rna
				ch4  = tuple[0] #set ch for next round
				k   = k-1
			else: #uniform sampling
				#print Tuples
				tuple = random.choice(Tuples)
				while ZF[(k-1,tuple[0])]==0: tuple = random.choice(Tuples)
				rna = tuple[1]+tuple[2]+tuple[3]+rna
				ch4  = tuple[0] #set ch for next round
				k   = k-1
		rna = ch4+rna #append the first character, which remains to be appended
		#~ for k in sorted(ZF):
			#~ print k,ZF[k]
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
		shNUCL = NUCL
		random.shuffle(shNUCL)
		for ch4 in shNUCL:
			if(ZF[k,gc,ch4]!=0):
				possible =1
				break
		if(possible==-1):
			print "no sequence with GC-content %d found!" %gc
			return SAMPLES
		#print gc,ZF[k,gc,ch4]
		#gc -= GCcont(ch4)
		#--------construct sample
		while k>0:
			#create list of 4-tuples having similarity to p[k],q[k] in frame 0,1
			#and whose last character is ch
			Tuples0 = []
			dimer   = peptide1[k-1]+peptide2[n-k] #peptides are 0-indexed
			for tuple in D[dimer]:
				if tuple[-1]==ch4: Tuples0.append(tuple)
			#-----------------------------------------------
			#remove tuples that don't satisfy constraints
			#-----------------------------------------------
			Tuples = []
			if constraint==None: 
				Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
			else:
				constraint4tuple = constraint[3*(k-1):3*k+1]
				for tuple in Tuples0:
					ok = 1
					for i in range(len(constraint4tuple)):
						ch = constraint4tuple[i]
						if tuple[i] not in iupacDict[ch]:
							ok = 0
							break
					if ok: 
						Tuples.append(tuple)
			#-----------------------------------------------
			#choose random tuple with first char ch1 for which ZF[(k-1,ch1)]>0
			#-----------------------------------------------
			#print 'tuples before filteration:'
			#print Tuples
			random.shuffle(Tuples)
			for tuple in Tuples:
				tupleGC = GCcont(tuple[1:])
				if((gc-tupleGC)>=0):
					#print "**",tuple,tupleGC,ZF[k-1,gc-tupleGC,tuple[0]]
					if k>1:
						if(ZF[k-1,gc-tupleGC,tuple[0]]!=0):
							break
					else:
						if(GCcont(tuple)==gc):
							break
			gc -= tupleGC
			#print (k-1),tuple,tuple[0],"+++"
			#ZF[(k,tuple[3])] -= 1.0
			rna = tuple[1]+tuple[2]+tuple[3]+rna
			ch4  = tuple[0] #set ch for next round
			k   = k-1
		rna = ch4+rna #append the first character, which remains to be appended
		#for k in sorted(ZF):
			#print k,ZF[k]
		SAMPLES.append(rna)
	return SAMPLES


def bfs_1(peptide1,peptide2,filename=None,threshold=1):
  #WARNING: This assumes that peptide1 in reading frame 0 and 
  #peptide2 in reading frame 1 of the opposite strands.
  #WARNING: Assume len(peptide1) = len(peptide2) 
#  assert (len(peptide1) == len(peptide2))
  lenPeptide = min(len(peptide1),len(peptide2)) 
  D     = createListOf4mersForAllPairsOfResidues(filename,threshold)
  dimer = peptide1[0]+peptide2[lenPeptide-1]
  #terminate with list of RNAsequences of len 3*lenPeptide+1
  RNA   = D[dimer]
  for i in range(1,lenPeptide):
    dimer = peptide1[i]+peptide2[lenPeptide-i-1]
    L     = D[dimer] 
    RNAbis= []
    for rna in RNA:
      for fourMer in L:
        if rna[-1]==fourMer[0]:
          RNAbis.append(rna+fourMer[1:])
    RNA = copy.deepcopy(RNAbis)
  return RNA

def main(peptide1,peptide2,numSamples,weighted,gc,replacement,filename=None,threshold=1,constraint=None):
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
		RNA    = sample(peptide1,peptide2,seqLen,numSamples,replacement,D,constraint,weighted)
	else:
		ZF     = computeZF_GC(peptide1,peptide2,D,constraint)			
		ZFval  = 0.0
		for ch in NUCL:
			for x in range(3*seqLen+2):ZFval += ZF[(seqLen,x,ch)]	
		RNA    = sample_GC(peptide1,peptide2,seqLen,numSamples,gc,D,constraint,weighted)
	if PRINT:
		if(VERBOSE): print "#total number of RNAs:%f" % ZFval
		printRNA(RNA)
		#~ for rna in (RNA):
			#~ print GCcont(rna)
	return ZFval,RNA

def parseConstraintFile(fname,n):
	if not (os.path.isfile(fname)):
		print "\nconstraint file does not exist!\n"
		sys.exit(1)
	else:
		f = open(fname,'r')
		const = f.readline().strip()
		if(VERBOSE): print '#constraints:',const
		if len(const)!= 3*n+1: 
			print "\nerror: constraint should have length 3n+1 where n is the length of the given peptides\n"
			sys.exit(1)
		for i in const:
			if i not in iupacDict.keys():
				print "\nconstraint file contains character %s which is not a nucleotide iupac code\n"
				sys.exit(1)
	return const
	
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
	replacement=1
	constraint = None
	args     = sys.argv[3:]
	#sequence constraint indicated by -c flag
	for i in range(len(args)):
		if args[i]=='-c':
			constraintfile = args[i+1]
			del[args[i+1]]; del[args[i]]
			break
		else:
			constraintfile = None
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
	#sampling with or w/o replacement with -r 1 / -r 0 flage respectively.
	for i in range(len(args)):
		if args[i]=='-r':
			replacement = int(args[i+1])
			if(replacement!=0 and replacement!=1):
				print "\nerror in -r input!"
				printusage(sys.argv[0])
			del[args[i+1]]; del[args[i]]
			break
	#weighted sampling, with -w flag
	for i in range(len(args)):
		if args[i]=='-w':
			weighted = 1
			break
	#verbose with -v flag
	for i in range(len(args)):
		if args[i]=='-v':
			VERBOSE = True
			break
	n   = len(peptide1)
	if(constraintfile):
		constraint = parseConstraintFile(constraintfile,n)
	else:
		constraint = None
	main(peptide1,peptide2,numSamples,weighted,gc,replacement,filename,threshold,constraint)
	#all = bfs_1(peptide1,peptide2,filename,threshold)
	#print len(all)
	#print len(all),all
