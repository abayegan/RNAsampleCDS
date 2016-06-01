#! /usr/bin/python2.6

# P.Clote - A. Bayegan
"""This program samples RNA sequences that translate into peptide1 and peptide2 in reading frames 0 in the negative and positive strands respectively.
Translation could be based on either exact or Blosum similar amino acids. In addition, RNA sequence constraints and weighted
sampling are supported. In weighted sampling nucleotide frequencies are representative of frequencies in the total partition function."""

import sys,os,copy,random
from math import log,exp
from utility import *

DEBUG  = False
VERBOSE= False
PRINT  = True #set False if use main as function to return list

def createListOfCodonsForAllPairsOfResidues(filename=None,threshold=1):
	if filename!=None:
		Blosum,AAs = readBlosumMatrix(filename)
	D = {}
	for aa1 in singleLetterAAcodes:
		for aa2 in singleLetterAAcodes:
			Tuples = [] 
			for gcode in genCode.keys():
				if (genCode[gcode].upper()!='STOP' and genCode[gcode[::-1]].upper()!='STOP' ):
					aa1bis = aaCode[genCode[gcode].upper()]
					aa2bis = aaCode[genCode[gcode[::-1]].upper()]
					if filename==None:
						if aa1==aa1bis and aa2==aa2bis:
							#print genCode,aa1bis,aa1bis,aa1,aa2
							Tuples.append(gcode)
					else:
						if Blosum[aa1][aa1bis]>=threshold and  Blosum[aa2][aa2bis]>=threshold: Tuples.append(gcode)
			if (Tuples):
				D[aa1+aa2] = Tuples
	#~ for i in sorted(D):
		#~ print i,D[i]
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
	# Function computes the ZF[k] = number of RNA sequences
	#    s = a_{0}...a_{3k-1} where len(s)=3*k  and such that
	# (1) translateMRNA(s) = peptide1[:k]
	# (2) translateMRNA(s[::-1]) = peptide2[n-k:]
	# The final coding sequences are of the form a_{0} ... a_{3n-1}.
	#
	# Parameter explanation:
	# D is dictionary of all compatible codons corresponding to peptide1 and peptide2
	# D is returned by function
	#    createListOfCodonsForAllPairsOfResidues(filename,threshold)
	#
	#WARNING: This assumes that peptide1 is in reading frame 0 of + strand and 
	#         peptide2 is in reading frame 0 of - strand.
	p = peptide1; q = peptide2
	assert (len(p) <= len(q)) 
	n  = min(len(p),len(q)) 
	#------------Initialize ZF---------------------
	ZF = {}
	for k in range(0,n+1):
		ZF[k] = 0.0
	ZF[0] = 1.0
	#~ else: #first position is constrained     
		#~ for ch in constraint[0]:
			#~ ZF[(0,ch)]=1.0
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(1,n+1): #k in [1,...,n] 
			dimer = p[k-1]+q[n-k]
			if dimer in D.keys():
				L     = D[dimer]
				#print k,dimer,D[dimer],L
				for s in L:
					ZF[k] += ZF[k-1]
	else: #must satisfy sequence constraint
		for k in range(1,n+1): #k in [1,...,n] 
			dimer = p[k-1]+q[n-k]
			if dimer in D.keys():
				L0               = D[dimer]
				constraintcodon = constraint[3*(k-1):3*k]
				L                = []
				for tuple in L0:
					ok = 1
					for i in range(len(constraintcodon)):
						ch = constraintcodon[i]
						if tuple[i] not in iupacDict[ch]:
							ok = 0
							#break
					if ok:
						L.append(tuple)
				if DEBUG:
					print "computeZF ",k,L
				for s in L:
					ZF[k] += ZF[k-1]
	return ZF

def computeZF_GC(peptide1,peptide2,D,constraint=None): 
	DEBUG=0
	p = peptide1; q = peptide2
	assert (len(p) <= len(q)) 
	n  = min(len(p),len(q)) 
	#------------Initialize ZF---------------------
	ZF = {}
	for k in range(0,n+1):
			for x in range(0,3*k+4):
				ZF[(k,x)] = 0.0
	if constraint==None or constraint[0]=='N':
		dimer = p[0]+q[k-1]
		L     = D[dimer]
		for s in L:
			print s
			ZF[(1,GCcont(s))] += 1.0
	else: #first tuple is constrained
		dimer = p[0]+q[k-1]
		L0     = D[dimer]
		constraint4tuple = constraint[:3]
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
				ZF[(1,GCcont(s))] += 1.0    
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(2,n+1):
			dimer = p[k-1]+q[n-k]
			if(dimer in D.keys()):
				L     = D[dimer]
				for s in L:
					print k,s
					for x in range(0,3*k+1):
						if(x-GCcont(s)>=0):
							ZF[(k,x)] += ZF[(k-1,x-GCcont(s))]		
	else: #must satisfy sequence constraint
		for k in range(2,n+1): #k in [1,...,n] 
			dimer            = p[k-1]+q[n-k]
			if(dimer in D.keys()):
				L0               = D[dimer]
				constraint4tuple = constraint[3*(k-1):3*k]
				L = []
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
					for x in range(0,3*k+1):
						if(x-GCcont(s)>=0):
							ZF[(k,x)] += ZF[(k-1,x-GCcont(s[1:]))]
	if(DEBUG):
		print "CG-ZF"
		for k in sorted(ZF):
			print k,ZF[k]
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
	assert( len(peptide1)==len(peptide2) )
	SAMPLES = []
	ZF      = computeZF(peptide1,peptide2,D,constraint)
	if weighted==1:
		print "#WARNING:weighted sampling for coding in frames 0 on the opposit strands does not differ from uniform sampling!"
	for i in range(numSamples):
		#initialize empty sequence and set k=n
		k  = seqLen; rna = ""
		#--------construct sample
		while k>0:
			Tuples0 = []
			dimer  = peptide1[k-1]+peptide2[n-k] #peptides are 0-indexed
			if dimer in D.keys():
				for tuple in D[dimer]:
					Tuples0.append(tuple)
			#-----------------------------------------------
			#remove tuples that don't satisfy constraints
			#-----------------------------------------------
			Tuples = []
			if constraint==None: 
				Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
			else:
				constraint4tuple = constraint[3*(k-1):3*k]
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
			tuple = random.choice(Tuples)
			if ZF[k-1]==0: print "no sequence coding in frames 0 on the opposit strands exists!"; return SAMPLES
			#while ZF[k-1]==0: tuple = random.choice(Tuples)
			rna = tuple[0]+tuple[1]+tuple[2]+rna
			k   = k-1
		SAMPLES.append(rna)
	return SAMPLES
			
def sample_GC(peptide1,peptide2,seqLen,numSamples,GC,D,constraint,weighted):
	assert( len(peptide1)==len(peptide2) )
	SAMPLES = []
	ZF      = computeZF_GC(peptide1,peptide2,D,constraint)
	if weighted==1:
		print "#WARNING:weighted sampling for coding in frames 0 on the opposit strands does not differ from uniform sampling!"
	for i in range(numSamples):
		k  =  seqLen;rna="";gc=GC
		while k>0:
			Tuples0 = []
			dimer   = peptide1[k-1]+peptide2[n-k] #peptides are 0-indexed
			if dimer in D.keys():
				for tuple in D[dimer]:
					Tuples0.append(tuple)
			#-----------------------------------------------
			#remove tuples that don't satisfy constraints
			#-----------------------------------------------
			Tuples = []
			if constraint==None: 
				Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
			else:
				constraint4tuple = constraint[3*(k-1):3*k]
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
			random.shuffle(Tuples)
			for tuple in Tuples:
				tupleGC = GCcont(tuple)
				if((gc-tupleGC)>=0):
					if k>1:
						if(ZF[(k-1,gc-tupleGC)]!=0):
							break
					else:
						if(GCcont(tuple)==gc):
							break
			gc -= tupleGC
			rna = tuple[0]+tuple[1]+tuple[2]+rna
			ch4  = tuple[0] #set ch for next round
			k   = k-1
		SAMPLES.append(rna)
	return SAMPLES

def bfs_0(peptide1,peptide2,filename=None,threshold=1):
	lenPeptide = min(len(peptide1),len(peptide2)) 
	D     = createListOfCodonsForAllPairsOfResidues(filename,threshold)
	RNA=['']
	for i in range(0,lenPeptide):
		dimer = peptide1[i]+peptide2[n-i-1]
		if dimer in D.keys():
			L     = D[dimer]
			#print i,L
		else:
			print "no sequence exists that codes the given peptides in reading frames 0 of opposite strands!"
			sys.exit(0)
		RNAbis= []
		for rna in RNA:
			for codon in L:
				RNAbis.append(rna+codon)
		RNA = copy.deepcopy(RNAbis)
	return RNA 

def main(peptide1,peptide2,numSamples,weighted,gc,replacement,filename=None,threshold=1,constraint=None):
	D      = createListOfCodonsForAllPairsOfResidues(filename,threshold)
	seqLen = min(len(peptide1),len(peptide2))
	if (gc==-1):
		ZF     = computeZF(peptide1,peptide2,D,constraint)
		RNA    = sample(peptide1,peptide2,seqLen,numSamples,replacement,D,constraint,weighted)
		ZFval = ZF[seqLen]
	else:
		ZF     = computeZF_GC(peptide1,peptide2,D,constraint)			
		ZFval  = 0.0
		for x in range(3*seqLen+1):ZFval += ZF[(seqLen,x)]	
		RNA    = sample_GC(peptide1,peptide2,seqLen,numSamples,gc,D,constraint,weighted)
	if VERBOSE:
		print "#total number of RNAs:%f" % ZFval
		printRNA(RNA)
		#~ for rna in (RNA):
			#~ print GCcont(rna)
	return ZFval,RNA

def check(peptide1,peptide2,RNA):
  for rna in RNA:
    if translateMRNA(rna) != peptide1:
      return False
    if translateMRNA(rna[::-1]) != peptide2:
      return False
  return True

def parseConstraintFile(fname,n):
	if not (os.path.isfile(fname)):
		print "\nconstraint file does not exist!\n"
		sys.exit(1)
	else:
		f = open(fname,'r')
		const = f.readline().strip()
		if(VERBOSE): print '#constraints:',const
		if len(const)!= 3*n: 
			print "\nerror: constraint should have length 3n where n is the length of the given peptides\n"
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
	#createListOfCodonsForAllPairsOfResidues()
	#~ D      = createListOfCodonsForAllPairsOfResidues(filename,threshold)
	#~ ZF     = computeZF(peptide1,peptide2,D,constraint)
	#~ for k in range(n+1):
		#~ sys.stdout.write( "\tZF[(%s)] = %s\n  " % (k,ZF[k]) )
	#~ print len(bfs_0(peptide1,peptide2,filename,threshold))
