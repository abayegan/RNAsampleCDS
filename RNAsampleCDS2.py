#! /usr/bin/env python

# A.Bayegan - P.Clote 
"""This program samples RNA sequences that translate simultaneously into 6 peptides in reading frames 0,1,2 in both sense and antisence strands.
The given peptides must follow AA Iupac code. In addition to Iupac code it is possible to completely nutralize the coding requirement for codons by 
using letter 'O' in the input peptides. In other words, using 'O' in the input peptide means that the corresponding codon may either code for an amino acid 
or a stop codon.
Translation could be based on either exact or Blosum similar amino acids. In addition this program is able to generate all the possible solutions(BFS),
create the exact PSSM profile, sample based on the given range of GC-content, sample RNA sequences compatible with the given nucleotide iupac constraints,
and sample RNA sequences that nucleotide frequencies are more representative of the total sequences(weighted sampling)."""

import sys,os,copy,random,math,shlex,inspect
from math import log,exp
from utility import *
from subprocess import Popen,PIPE

VERBOSE = 0
DEBUG  = 0
PRINT  = True #set False if use main as function to return list
	
def createListOfOverlappingCodons():
	#create list L of ALL 5-tuples corresponding to overlapping codons
	CODONS = createAllCodons()
	L      = []
	for cod1 in CODONS:
		for cod2 in CODONS:
			if cod1[-1]==cod2[0]:
				L.append(cod1+cod2[1:])
	#print len(L)
	#print 'number of all possible 5-tuples(no coding requirement)', len(L)
	return L

def createListOfCompatible5mers(aalist,filename=None,threshold=1):
	#create list of 5mers which are compatible with the coding requirements in the given reading frames.
	#readingFramesL contains 6 items including the coding requirements in +0,+1,+2,-0,-1,-2 reading frames, respectively. 
	assert(len(aalist)== 6)
	if filename!=None:
		Blosum,AAs = readBlosumMatrix(filename)
		Blosum['O']={}
		for aa in AAs:
			Blosum['O'][aa]= 100
			Blosum[aa]['O']= 100
	#create list L of ALL 5-tuples corresponding to overlapping codons
	allTuples = createListOfOverlappingCodons()
	#D[dimer] is list of codon pairs that agree on overlap
	Tuples = [] 	
	for t in allTuples:
		aa0t = codon2aa[t[:3]]
		aa1t = codon2aa[t[1:4]]
		aa2t = codon2aa[t[2:]]
		aa3t = codon2aa[t[:3][::-1]]
		aa4t = codon2aa[t[1:4][::-1]]
		aa5t = codon2aa[t[2:][::-1]]
		if filename==None:
			if iscompatibleAAWithIupac(aa0t,aalist[0]) and iscompatibleAAWithIupac(aa1t,aalist[1]) and \
			 iscompatibleAAWithIupac(aa2t,aalist[2]) and iscompatibleAAWithIupac(aa3t,aalist[3]) \
			 and iscompatibleAAWithIupac(aa4t,aalist[4]) and iscompatibleAAWithIupac(aa5t,aalist[5]):
				Tuples.append(t)
		else: #use Blosum similarity; recall if aa==aabis then Blosum similar
			if (Blosum[aalist[0]][aa0t]>=threshold and Blosum[aalist[1]][aa1t]>=threshold and Blosum[aalist[2]][aa2t]>=threshold and \
			    Blosum[aalist[3]][aa3t]>=threshold and Blosum[aalist[4]][aa4t]>=threshold and Blosum[aalist[5]][aa5t]>=threshold):
				Tuples.append(t)
	#print Tuples
	return Tuples

def computeZB(peptides,constraint=None,filename=None,threshold=1):
	n  = len(peptides['+0'])
	#------------Initialize ZF---------------------
	ZB = {}
	for k in range(0,n+1):
		for ch1 in NUCL:
			for ch2 in NUCL:
				ZB[(k,ch1,ch2)] = 0.0
	if constraint==None:
		for ch1 in NUCL:
			for ch2 in NUCL:
				ZB[(n,ch1,ch2)] = 1.0
	else:
		for ch2 in iupacNuc[constraint[n-1]]:
			for ch1 in iupacNuc[constraint[n-2]]:
				ZB[(n,ch1,ch2)] = 1.0
	#------------Fill matrix ZF-------------------
	for k in range(n-1,-1,-1): #k in [n-1,...,0] 
		hexa = [peptides['+0'][k], peptides['+1'][k], peptides['+2'][k], peptides['-0'][n-k-1],peptides['-1'][n-k-1],peptides['-2'][n-k-1]]
		L0 = createListOfCompatible5mers(hexa,filename,threshold)
		if constraint==None:
			L = L0
		else:
			L=[];
			tupleconst = constraint[3*k:3*k+5]
			for tuple in L0:
				ok=1
				for i in range(len(tupleconst)):
					if  tuple[i] not in iupacNuc[tupleconst[i]]:
						ok=0
				if ok==1:
					L.append(tuple)
		for s in L:
			ZB[(k,s[0],s[1])] += ZB[(k+1,s[-2],s[-1])]
	
	if (DEBUG):
		Zval  = 0.0
		for ch1 in NUCL:
			for ch2 in NUCL:
				Zval += ZB[(0,ch1,ch2)]			
		print 'Total number of seuqneces from ZB: %f'%Zval
	#~ for key in sorted(ZB):
		#~ print key,ZB[key]
	return ZB

def computeZB_GC(peptides,constraint=None,filename=None,threshold=1):
	n  = len(peptides['+0'])
	#------------Initialize ZB---------------------
	ZB = {}
	for k in range(0,n+1):
		for ch1 in NUCL:
			for ch2 in NUCL:
				for x in range(0,3*(n-k)+8):
					ZB[(k,x,ch1,ch2)] = 0.0
	if constraint==None:
		hexa = [peptides['+0'][n-1], peptides['+1'][n-1], peptides['+2'][n-1], peptides['-0'][0],peptides['-1'][0],peptides['-2'][0]]
		L = createListOfCompatible5mers(hexa,filename,threshold)
		#print hexa , L
		for s in L:
			ZB[(n-1,GCcont(s),s[0],s[1])] += 1.0
	else: 
		hexa = [peptides['+0'][n-1], peptides['+1'][n-1], peptides['+2'][n-1], peptides['-0'][0],peptides['-1'][0],peptides['-2'][0]]
		L0     = createListOfCompatible5mers(hexa,filename,threshold) 
		constraint5tuple = constraint[-5:]
		L=[]
		for tuple in L0:
			ok = 1
			for i in range(len(constraint5tuple)):
				ch = constraint5tuple[i]
				if tuple[i] not in iupacNuc[ch]:
					ok = 0
					#break
			if ok: 
				L.append(tuple)
		for s in L:
				ZB[(n-1,GCcont(s),s[0],s[1])] += 1.0  
	#------------Fill matrix ZB-------------------
	for k in range(n-2,-1,-1): #k in [n-1,...,0] 
		hexa = [peptides['+0'][k], peptides['+1'][k], peptides['+2'][k], peptides['-0'][n-k-1],peptides['-1'][n-k-1],peptides['-2'][n-k-1]]
		L0 = createListOfCompatible5mers(hexa,filename,threshold)
		if constraint==None:
			L = L0
		else:
			L=[];
			tupleconst = constraint[3*k:3*k+5]
			for tuple in L0:
				ok=1
				for i in range(len(tupleconst)):
					if  tuple[i] not in iupacNuc[tupleconst[i]]:
						ok=0
				if ok==1:
					L.append(tuple)
		for s in L:
			for x in range(0,3*(n-k)+3):
				ZB[(k,x+GCcont(s[0:3]),s[0],s[1])] += ZB[(k+1,x,s[-2],s[-1])]
			#print L,k
			#for a in NUCL:
				#for b in NUCL:
					#for x in range(0,3*(n-k)+3):
						#print 'ZB[%d,%d,%s,%s]:%d' % (k,x,a,b,ZB[(k,x,a,b)])
	if (DEBUG):
		Zval  = 0.0
		for ch1 in NUCL:
			for ch2 in NUCL:
				for x in range(3*(n-k)+3):
					Zval += ZB[(0,x,ch1,ch2)]			
		print '#total number of seuqneces from ZB_GC: %f'%Zval
	return ZB
	

def computeZF(peptides,constraint=None,filename=None,threshold=1): 
	#forwards "partition function"
	# Function computes the ZF[(k,ch4,ch5)] = number of RNA sequences
	#    s = a_{0}...a_{3k} where len(s)=3*k+2  and such that
	# (1) translateMRNA(s) = peptide1[:k] 
	# (2) translateMRNA(s[2:]) = peptide2[:k]
	# (3) s[-1] = ch5
	# (4) s[-2] = ch4
	# The final coding sequences are of the form a_{0} ... a_{3n}.
	
	n  = len(peptides['+0'])
	#------------Initialize ZF---------------------
	ZF = {}
	for k in range(0,n+1):
		for ch3 in NUCL:
			for ch4 in NUCL:
				ZF[(k,ch3,ch4)] = 0.0
	if constraint==None:
		for ch3 in NUCL:
			for ch4 in NUCL:
				ZF[(0,ch3,ch4)] = 1.0
	else:
		for ch3 in iupacNuc[constraint[0]]:
			for ch4 in iupacNuc[constraint[1]]:
				ZF[(0,ch3,ch4)] = 1.0
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(1,n+1): #k in [1,...,n] 
			hexa = [peptides['+0'][k-1], peptides['+1'][k-1], peptides['+2'][k-1], peptides['-0'][n-k],peptides['-1'][n-k],peptides['-2'][n-k]]
			L = createListOfCompatible5mers(hexa,filename,threshold)
			for s in L:
				ZF[(k,s[-2],s[-1])] += ZF[(k-1,s[0],s[1])]
	else: #must satisfy sequence constraint
		for k in range(1,n+1): #k in [1,...,n] 
			hexa = [peptides['+0'][k-1], peptides['+1'][k-1], peptides['+2'][k-1], peptides['-0'][n-k],peptides['-1'][n-k],peptides['-2'][n-k]]
			L0 = createListOfCompatible5mers(hexa,filename,threshold)
			constraint5tuple = constraint[3*(k-1):3*k+2]
			L                = []
			for tuple in L0:
				ok = 1
				for i in range(len(constraint5tuple)):
					ch = constraint5tuple[i]
					if tuple[i] not in iupacNuc[ch]:
						ok = 0
						#break
				if ok: 
					L.append(tuple)
			if DEBUG:
				print "computeZF ",k,L
			for s in L:
				ZF[(k,s[-2],s[-1])] += ZF[(k-1,s[0],s[1])]
	if DEBUG:
		print "computing ZF:"
		for k in range(3*n+3):
			for ch3 in NUCL:
				for ch4 in NUCL:
					if((k,ch3,ch4) in ZF.keys()):
						sys.stdout.write( "\tZF[(%s,%s,%s)] = %s\n" % (k,ch3,ch4,ZF[(k,ch3,ch4)]) )
	if(VERBOSE):
		ZFval  = 0.0
		for ch3 in NUCL:
			for ch4 in NUCL:
				ZFval += ZF[(n,ch3,ch4)]
		print "#total number of sequences: %f" %ZFval
	return ZF

def computeZF_GC(peptides,constraint=None,filename=None,threshold=1): 
	#forwards "partition function"
	# Here, peptide1 = p_{0}...p_{n-1} and peptide2 = q_{0}...q_{n-1}
	# Function computes the ZF[(k,x,ch4,ch5)] = number of RNA sequences
	#    s = a_{0}...a_{3k} where len(s)=3*k+2  and such that
	# (1) translateMRNA(s) = peptide1[:k] 
	# (2) translateMRNA(s[2:]) = peptide2[:k]
	# (3) s[-1] = ch5
	# (4) s[-2] = ch4
	# (5) x number of G's and C's in s

	n  = len(peptides['+0'])
	#------------Initialize ZF---------------------
	ZF = {}
	for k in range(0,n+1):
		for ch3 in NUCL:
			for ch4 in NUCL:
				for x in range(0,3*k+8): #Warning: 3k+2 doesn't work.
					ZF[(k,x,ch3,ch4)] = 0.0
	if constraint==None:
		hexa = [peptides['+0'][0], peptides['+1'][0], peptides['+2'][0], peptides['-0'][n-1],peptides['-1'][n-1],peptides['-2'][n-1]]
		L = createListOfCompatible5mers(hexa,filename,threshold)
		#print hexa , L
		for s in L:
			ZF[(1,GCcont(s),s[-2],s[-1])] += 1.0
		if(DEBUG==1):
			print "tuples at level 1:", L
			for a in NUCL:
				for b in NUCL:
					for x in range(0,7):
						print 'ZF[1,%d,%s,%s]:%d' % (x,a,b,ZF[(1,x,a,b)])
	else: #first tuple is constrained
		hexa = [peptides['+0'][0], peptides['+1'][0], peptides['+2'][0], peptides['-0'][n-1],peptides['-1'][n-1],peptides['-2'][n-1]]
		L0     = createListOfCompatible5mers(hexa,filename,threshold) 
		constraint5tuple = constraint[:5]
		L=[]
		for tuple in L0:
			ok = 1
			for i in range(len(constraint5tuple)):
				ch = constraint5tuple[i]
				if tuple[i] not in iupacNuc[ch]:
					ok = 0
					#break
			if ok: 
				L.append(tuple)
		for s in L:
				ZF[(1,GCcont(s),s[-2],s[-1])] += 1.0  
	#------------Fill matrix ZF-------------------
	if constraint==None:
		for k in range(2,n+1):
			hexa = [peptides['+0'][k-1], peptides['+1'][k-1], peptides['+2'][k-1], peptides['-0'][n-k],peptides['-1'][n-k],peptides['-2'][n-k]]
			L = createListOfCompatible5mers(hexa,filename,threshold)
			if(DEBUG==1):
				print "tuples at level %d:" %k , L
			for s in L:
				for x in range(0,3*k+3):
					#print s,x,GCcont(s[1:4])
					if(x-GCcont(s[2:])>=0):
						ZF[(k,x,s[-2],s[-1])] += ZF[(k-1,x-GCcont(s[2:]),s[0],s[1])]
			if(DEBUG==1):
				for a in NUCL:
					for b in NUCL:
						for x in range(0,3*k+3):
							print 'ZF[%d,%d,%s,%s]:%d' % (k,x,a,b,ZF[(k,x,a,b)])
		#ZFval  = 0.0
		#for ch in NUCL: ZFval += ZF1[(n,ch)]
		#print "%f" %ZFval
		
	else: #must satisfy sequence constraint
		for k in range(1,n+1): #k in [1,...,n] 
			hexa = [peptides['+0'][k-1], peptides['+1'][k-1], peptides['+2'][k-1], peptides['-0'][n-k],peptides['-1'][n-k],peptides['-2'][n-k]]
			L0 = createListOfCompatible5mers(hexa,filename,threshold)
			constraint5tuple = constraint[3*(k-1):3*k+2]
			L                = []
			for tuple in L0:
				ok = 1
				for i in range(len(constraint5tuple)):
					ch = constraint5tuple[i]
					if tuple[i] not in iupacNuc[ch]:
						ok = 0
						#break
				if ok: 
					L.append(tuple)
			for s in L:
				for x in range(0,3*k+3):
					if(x-GCcont(s[2:])>=0):
						ZF[(k,x,s[-2],s[-1])] += ZF[(k-1,x-GCcont(s[2:]),s[0],s[1])]
			if(DEBUG==1):
				for a in NUCL:
					for b in NUCL:
						for x in range(0,3*k+3):
							print 'ZF[%d,%d,%s,%s]:%d' % (k,x,a,b,ZF[(k,x,a,b)])
	ZFval  = 0.0
	for ch3 in NUCL:
		for ch4 in NUCL:
			for x in range(3*n+3):
				ZFval += ZF[(n,x,ch3,ch4)]	
	if (VERBOSE) : print "#total number of RNAs:%f" % ZFval
	if (DEBUG) : print "#total number of RNAs from ZF_GC:%f" % ZFval 
	return ZF

def getNucSampleProb(ZF,k,nuclist):
	probdict={}
	sum=0.0
	if len(nuclist[0])==1: #probabity of sampling for a single nucleotide
		for nuc1 in nuclist:
			for nuc2 in nuclist:
				sum += ZF[k,nuc1,nuc2]
		for nuc1 in nuclist:
			for nuc2 in nuclist:
				probdict[(nuc1,nuc2)]= ZF[k,nuc1,nuc2]/sum
	else: #sampling probability of dinucleotide
		for nuc in nuclist:
			sum += ZF[k,nuc[0],nuc[1]]
		for nuc in nuclist:
			probdict[(nuc[0],nuc[1])]= ZF[k,nuc[0],nuc[1]]/sum
	return probdict

def sample(peptides,ZF,n,numSamples,constraint,weighted,filename,threshold):
	SAMPLES = []
	for i in range(numSamples):
		#initialize empty sequence and set k=n
		k  = n; rna = ""
		if DEBUG: print "sampling sequence %d"%i
		#--------choose random nucleotide for which ZF[(k,ch)]>0 
		if weighted==1:
			sampleDict = getNucSampleProb(ZF,k,NUCL)
			ch3,ch4 = rouletteWheel(sampleDict)
			while ZF[(k,ch3,ch4)]==0: 
				ch3,ch4 = rouletteWheel(sampleDict)
			#print 'selected ch4',ch4	
		else: #uniform sampling
			found = 0
			shNUCL1 = NUCL
			shNUCL2 = NUCL
			random.shuffle(shNUCL1)
			random.shuffle(shNUCL2)
			for ch3 in shNUCL1:
				for ch4 in shNUCL2:
					if(ZF[k,ch3,ch4]>0):
						found =1
						break
				if found: break
			if not found:
				print "\nno sequence satisfying the coding requirement exists!"
				sys.exit(0)
		#--------construct sample
		while k>0:
			#create list of 5-tuples having similarity to p[k],q[k] in frame 0,2
			#and whose last characters are ch3,ch4
			Tuples0 = []
			hexa = [peptides['+0'][k-1], peptides['+1'][k-1], peptides['+2'][k-1], peptides['-0'][n-k],peptides['-1'][n-k],peptides['-2'][n-k]]
			L = createListOfCompatible5mers(hexa,filename,threshold)
			for tuple in L:
				if (tuple[-1]==ch4 and tuple[-2]==ch3): Tuples0.append(tuple)
			#-----------------------------------------------
			#remove tuples that don't satisfy constraints
			#-----------------------------------------------
			Tuples = []
			if constraint==None: 
				Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
			else:
				constraint5tuple = constraint[3*(k-1):3*k+2]
				for tuple in Tuples0:
					ok = 1
					for i in range(len(constraint5tuple)):
						ch = constraint5tuple[i]
						if tuple[i] not in iupacNuc[ch]:
							ok = 0
							break
					if ok:
						Tuples.append(tuple)
			#-----------------------------------------------
			#choose random tuple with first char ch1 and second charater ch2 for which ZF[(k-1,ch1,ch2)]>0
			#-----------------------------------------------
			#print 'tuples before filteration:'
			#print k,Tuples
			random.shuffle(Tuples)
			if weighted==1:
				selTuples=[];nuclist=[]
				for tup in Tuples: #make sure the sampled tuples are compatible with k-1 (tup[0]=ch1 and tup[1]=ch2) 
						if tup[0:2] not in nuclist: nuclist.append(tup[0:2])
				#print 'nuclist', nuclist
				sampleDict = getNucSampleProb(ZF,k-1,nuclist)
				#print 'weighted sampling dict:', sampleDict
				ch1,ch2 = rouletteWheel(sampleDict)
				while ZF[(k-1,ch1,ch2)]==0: 
					ch1,ch2 = rouletteWheel(sampleDict)
				#print 'sampled characters:',ch1,ch2
				for tup in Tuples:
					if(tup[0]==ch1 and tup[1]==ch2 and k!=1): 
						selTuples.append(tup)
					if(k==1): #Warning: since ZF[0,nuc]=1 (nuc in {A,C,G,U}) the first nucleotide does not need to be ch1 and can be anything 
						selTuples = Tuples
				#print 'tuples after filteration:'
				#print selTuples
				tuple = random.choice(selTuples)
				#print "selected tuple:",tuple
			else: #uniform sampling
				#print "before filtering" , Tuples
				tuple = random.choice(Tuples)
				while ZF[(k-1,tuple[0],tuple[1])]==0: tuple = random.choice(Tuples)
				#print "after filtering", tuple
			rna = tuple[2]+tuple[3]+tuple[4]+rna
			ch4  = tuple[1] #set ch for next round
			ch3  =  tuple[0]
			k   = k-1
		rna = ch3+ch4+rna #append the first character, which remains to be appended
		#for k in sorted(ZF):
			#print k,ZF[k]
		SAMPLES.append(rna)
	return SAMPLES
			
def sample_GC(peptides,ZF,n,numSamples,GClow,GCup,constraint,weighted,filename,threshold):
	SAMPLES = []
	for i in range(numSamples):
		#initialize empty sequence and set k=n
		k  = n; rna = ""; possible=-1
		gc1= int(GClow) #gc1: lower bound inclusive 
		gc2= int(GCup+1) #gc2: upper bound exclusive
		#----------------------------------------------------------------------------
		#choose random nucleotide and random gc in GC-range for which ZF[k,gc,ch3,ch4]>0
		#----------------------------------------------------------------------------
		found = 0
		shNUCL1 = NUCL
		shNUCL2 = NUCL
		gcRange = range(gc1,gc2)
		random.shuffle(shNUCL1)
		random.shuffle(shNUCL2)
		random.shuffle(gcRange)
		#print gcRange
		for ch3 in shNUCL1:
			for ch4 in shNUCL2:
				for gc in gcRange:
					if(ZF[(n,gc,ch3,ch4)]!=0):
						found=1
						break
				if found: break
			if found: break	
		if not found:
			print "no sequence with GC-content in range [%0d , %d] found!" %(GClow,GCup)
			return SAMPLES
		#gc -= GCcont(ch4)
		#--------construct sample
		while k>0:
			#---------------------------------------------------
			#find all 5-tuples satisfying the coding requiremnt
			#peptides[i][k-1] for i in {+0,+1,+2} and peptides[i][n-k] for i in {-0,-1,-2}
			#---------------------------------------------------
			Tuples0 = []
			hexa = [peptides['+0'][k-1], peptides['+1'][k-1], peptides['+2'][k-1], peptides['-0'][n-k],peptides['-1'][n-k],peptides['-2'][n-k]]
			L = createListOfCompatible5mers(hexa,filename,threshold)
			for tuple in L:
				if (tuple[-1]==ch4 and tuple[-2]==ch3): Tuples0.append(tuple)
			#-----------------------------------------------
			#remove tuples that don't satisfy constraints
			#-----------------------------------------------
			Tuples = []
			if constraint==None: 
				Tuples = Tuples0 #Warning: this sets Tuples to point to Tuples0
			else:
				constraint5tuple = constraint[3*(k-1):3*k+2]
				for tuple in Tuples0:
					ok = 1
					for i in range(len(constraint5tuple)):
						ch = constraint5tuple[i]
						if tuple[i] not in iupacNuc[ch]:
							ok = 0
							break
					if ok: 
						Tuples.append(tuple)
			#-----------------------------------------------
			#choose random tuple with first char ch1 and second char ch2 for which ZF[(k-1,ch1)]>0
			#-----------------------------------------------
			random.shuffle(Tuples) #important
			for tuple in Tuples:
				tupleGC = GCcont(tuple[2:])
				if((gc-tupleGC)>=0):
					if k>1:
						if(ZF[k-1,gc-tupleGC,tuple[0],tuple[1]]!=0):
							break
					else:
						if(GCcont(tuple)==gc):
							break
			gc -= tupleGC
			#print (k-1),tuple,tuple[0],tuple[1]
			rna = tuple[2]+tuple[3]+tuple[4]+rna
			ch4  = tuple[1] #set ch for next round
			ch3  =  tuple[0]
			k   = k-1
		rna = ch3+ch4+rna #append the first character, which remains to be appended
		#for k in sorted(ZF):
			#print k,ZF[k]
		SAMPLES.append(rna)
	return SAMPLES

def bfs(peptides,filename=None,threshold=1,constraint=None):
	n = len(peptides['+0'])
	hexa = [peptides['+0'][0],peptides['+1'][0],peptides['+2'][0],peptides['-0'][n-1],peptides['-1'][n-1],peptides['-2'][n-1]]
	L0 = createListOfCompatible5mers(hexa,filename,threshold)
	RNA=[]
	if constraint==None:
		RNA = L0
	else:
		const = constraint[0:5]
		for t in L0:
			ok=1
			for i in range(5):
				#print t,i,t[i]
				if t[i] not in iupacNuc[const[i]]:
					ok=0
			if ok: RNA.append(t)
	for k in range(1,n):
		hexa = [peptides['+0'][k], peptides['+1'][k], peptides['+2'][k], peptides['-0'][n-k-1],peptides['-1'][n-k-1],peptides['-2'][n-k-1]]
		L0 = createListOfCompatible5mers(hexa,filename,threshold)
		RNAbis = []
		L=[]
		if constraint==None:
			L=L0
		else:
			const = constraint[3*k:3*k+5]
			for t in L0:
				ok=1
				for i in range(len(const)):
					if t[i] not in iupacNuc[const[i]]:
						ok=0
				if ok:
					L.append(t)
				#print L
		for rna in RNA:
			for tup in L:
				if rna[-2]==tup[0] and rna[-1]==tup[1]:
					RNAbis.append(rna+tup[2:])
		RNA = copy.deepcopy(RNAbis)		
	if(VERBOSE): print "#total number of sequences:",len(RNA)
	return RNA

def createPSSM(peptides,ZF,ZB,filename=None,threshold=1):
	# Function computes PSSM[(k,ch) = frequency of ch in position k of the
	# multiple alignment of ALL sequences
	#    s = a_{0}...a_{3n+2}  such that
	# translateMRNA(s[:-2])=peptides[0],translateMRNA(s[1:-1])=peptide[1],translateMRNA(s[2:])=peptide[2]
	# translateMRNA(s[::-1][:-2])=peptides[3],translateMRNA(s[::-1][1:-1])=peptides[4] and translateMRNA(s[::-1][2:])=peptides[3]
	# ZF[(k,ch)] is the number of RNA sequences that satisfy the coding requirements in all reading frames and s[-1] = ch 
	# ZB[(k,ch)] is the number of RNA sequences
	#    s = a_{3k}...a_{3n+2}such that satisfy the coding requirements in peptides[i][k:n] (0<=i<=5) and s[0]=ch
	
	n  = len(peptides['+0'])
	ZFtotal = 0.0; ZBtotal = 0.0
	for ch1 in NUCL:
		for ch2 in NUCL:
			ZFtotal += ZF[(n,ch1,ch2)]
			ZBtotal += ZB[(0,ch1,ch2)]
	assert (ZFtotal == ZBtotal)
	#----------Initialize PSSM ---------------
	PSSM = {}
	for i in range(3*n+2):
		for ch in NUCL:
			PSSM[(i,ch)] = 0.0
	#----------Compute PSSM(3*k+i,ch) for i=0,1,2 ---------------
	for k in range(n):
		hexa = [peptides['+0'][k], peptides['+1'][k], peptides['+2'][k], peptides['-0'][n-k-1],peptides['-1'][n-k-1],peptides['-2'][n-k-1]]
		L = createListOfCompatible5mers(hexa,filename,threshold)
		for s in L:
			for ch in NUCL:
				#~ if(s[0]==ch and k==0):
					#~ PSSM[0,ch] += ZB[(k+1,s[0],s[1])]
				#~ if(s[1]==ch and k==0):				
						#~ PSSM[0,ch] = ZB[(1,s[0],s[1])]
				if(s[0]==ch ):
					PSSM[3*k,ch]+= ZF[(k,s[0],s[1])] * ZB[(k+1,s[3],s[4])]
				if(s[1]==ch):
					PSSM[3*k+1,ch] += ZF[(k,s[0],s[1])]* ZB[(k+1,s[3],s[4])]
				if(s[2]==ch):
					PSSM[3*k+2,ch] += ZF[k,s[0],s[1]]*ZB[(k+1,s[3],s[4])]
	#~ for nuc1 in NUCL:
		#~ for nuc2 in NUCL:
			#~ print "ZF(%d,%s,%s): %f"%(n,nuc1,nuc2,ZF[(n,nuc1,nuc2)])
	for ch1 in NUCL:
		for ch2 in NUCL:
			PSSM[3*n,ch1] += ZF[(n,ch1,ch2)]
			PSSM[3*n+1,ch2] += ZF[(n,ch1,ch2)]
	#----------Normalize PSSM ---------------
	for i in range(3*n+2): #i in [0,...,3n+1]
		for ch in NUCL:
			PSSM[(i,ch)] = PSSM[(i,ch)]
	return PSSM


def logoplot(PSSM,outname):
# PSSM is a dictionary s.t PSSM[(pos,ch)]= frequency of nucleotide 'ch' at position 'pos' of the sequence.
	fout = open('tmp','w')
	header = "ID Matrix\nBF\nPO\tA\tC\tG\tU\n"
	fout.write(header)
	for i in range(len(PSSM.keys())/4):
		As = PSSM[(i,'A')]
		Cs = PSSM[(i,'C')]
		Gs = PSSM[(i,'G')]
		Us = PSSM[(i,'U')]
		fout.write("%d\t%d\t%d\t%d\t%d\n" % (i+1,As,Cs,Gs,Us))
	fout.close()
	pypath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
	cmd = "%s -F pdf -f %s --errorbars NO -o %s" %(pypath+"/weblogo/weblogo",'tmp',outname+'.pdf')
	run = Popen(cmd.split())
	run.communicate()
	#os.remove('tmp')

def main(peptides,numSamples,weighted,gclow,gcup,BFS,PSSMfilename,filename=None,threshold=1,constraint=None):
	seqLen = len(peptides['+0'])
	#convert GC-content to GC-count
	gc1 = int((3*seqLen+2)*gclow/100)
	gc2 = math.ceil((3*seqLen+2)*gcup/100)
	if(VERBOSE):
		#~ for k in sorted(peptides):
			#~ print "#%s %s"%(k,peptides[k])
		print "#length of the input peptides:%d" % seqLen
		#~ if(constraint):
				#~ print "#nucleotide constraint:",constraint		
	if(BFS): #output all of the sequences
		RNA = bfs(peptides,filename,threshold,constraint)
		if(PSSMfilename):
			ZF      = computeZF(peptides,constraint,filename,threshold)
			ZB = computeZB(peptides,constraint,filename,threshold)
			pssm = createPSSM(peptides,ZF,ZB,filename,threshold)
			logoplot(pssm,PSSMfilename)
			
	else:
		if (gclow==-1):
			ZF      = computeZF(peptides,constraint,filename,threshold)
			RNA    = sample(peptides,ZF,seqLen,numSamples,constraint,weighted,filename,threshold)				
			if(PSSMfilename):
				ZB = computeZB(peptides,constraint,filename,threshold)
				pssm = createPSSM(peptides,ZF,ZB,filename,threshold)
				logoplot(pssm,PSSMfilename)
		else:
			if(VERBOSE): print "#GC-Content: %d-%d percent"%(gclow,gcup)
			ZF      = computeZF_GC(peptides,constraint,filename,threshold)
			RNA    = sample_GC(peptides,ZF,seqLen,numSamples,gc1,gc2,constraint,weighted,filename,threshold)
			if(PSSMfilename):
				if(VERBOSE): print "#WARNING: GC-content is not used in creating PSSM."
				ZF      = computeZF(peptides,constraint,filename,threshold)
				ZB = computeZB(peptides,constraint,filename,threshold)
				pssm = createPSSM(peptides,ZF,ZB,filename,threshold)
				logoplot(pssm,PSSMfilename)
	#computeZB_GC(peptides,constraint)
	printRNA(RNA)
	#~ for rna in (RNA):
		#~ print GCcont(rna)	
	return RNA

def parseConstraintFile(fname,n):
	if not (os.path.isfile(fname)):
		print "\nERROR:constraint file does not exist!\n"
		sys.exit(1)
	else:
		f = open(fname,'r')
		const = f.readline().strip()
		if len(const)< 3*n+2: 
			const = const + 'N'*(3*n+2-len(const))
		if len(const) != 3*n+2: 
			print "\nERROR: constraint should have length of maximum 3n+2 where n is the length of the given peptides\n"
			sys.exit(1)
		for i in const:
			if i not in iupacNuc.keys():
				print "ERROR:\nconstraint file contains character %s which is not a nucleotide IUPAC code\n" %i
				sys.exit(1)
	return const

def parsePeptideFile(peptidefile):
	if (not os.path.isfile(peptidefile)):
		print "peptidefile does not exist!\n"
		sys.exit(1)
	d={}
	frameList = ['+0','+1','+2','-0','-1','-2']
	f = open(peptidefile,'r')
	lines = f.readlines()
	for line in lines:
		found=0
		if len(line)<4: 
			print "ERROR: error in peptide file!\n"
			sys.exit(1)
		rf = line.split()[0]
		pep = line.split()[1]
		for i in frameList:
			if rf == i:
				d[i] = pep.upper()
				found=1
				break
		if not found:
			print 'ERROR: error in frames of peptidefile!\n'
			sys.exit(1)
	n = len(pep)
	for rf in d.keys():
		if len(d[rf])!= n:
			print 'ERROR: all the peptides must have equal length!\n'
			sys.exit(1)
		for ch in d[rf]:
			if ch not in IupacAA:
				print 'ERROR: error in peptide sequences!\n'
				sys.exit(1)
	for rf in frameList:
		if rf not in d.keys():
			d[rf] = 'O'*n
	return d

if __name__ == '__main__':
	if len(sys.argv) < 3:
		printusage(sys.argv[0])
	args     = sys.argv[1:]
	
	#------------------default setting------------------
	numSamples=10
	weighted = -1
	gclow=-1;gcup=-1
	BFS=0
	threshold=1
	PSSMfilename=None
	filename=None
	constraintfile=None
	peptidefile=None
	#------------------parse input arguments-------------
	
	#peptide file indicated by '-p' flag
	for i in range(len(args)):
		if args[i]=='-p':
			if args[i+1]:
				if i+1!=len(args) and args[i+1][0]!='-':
					peptidefile = args[i+1]
					del[args[i+1]]; del[args[i]]
					break
				else:
					print "ERROR: peptide file must be secified after -p flag!\n "
					printusage(sys.argv[0])
	#sequence constraint indicated by -c flag
	for i in range(len(args)):
		if args[i]=='-c':
			if i+1!=len(args) and args[i+1][0]!='-':
				constraintfile = args[i+1]
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: constraint file must be secified after -c flag!\n "
				printusage(sys.argv[0])
	#number of samples, with -n flag
	for i in range(len(args)):
		if args[i]=='-n':
			if i+1!=len(args) and args[i+1][0]!='-':
				if args[i+1].isdigit():
					numSamples = int(args[i+1])
				elif args[i+1]=='all':
					BFS = 1
				else: 
					print "ERROR: error in -n flag!\n"
					printusage(sys.argv[0])
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: number of samples must be specified after -n flag!\n"
				printusage(sys.argv[0])
	#similarity matrix filename, with -m flag
	for i in range(len(args)):
		if args[i]=='-m':
			if i+1!=len(args) and args[i+1][0]!='-':
				filename = args[i+1]
				del[args[i+1]]; del[args[i]]
				break
			else: 
				print "ERROR: matrix file must be specified after -m flag!\n"
				printusage(sys.argv[0])
	#blosum threshold, with -t flag
	for i in range(len(args)):
		if args[i]=='-t':
			if i+1!=len(args) and args[i+1][0]!='-':
				threshold = int(args[i+1])
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: similarity threshold must be specified after -t flag!\n"
				printusage(sys.argv[0])
	#GC content with -gc flag
	for i in range(len(args)):
		if args[i]=='-gc':
			if i+1!=len(args) and args[i+1][0]!='-':
				gcinput = args[i+1].split(':')
				if(len(gcinput)==1):
					gclow = int(gcinput[0])
					gcup = gclow
				elif(len(gcinput)==2):
					gclow = int(gcinput[0])
					gcup = int(gcinput[1])
					if gclow>gcup or gclow<0:
						print "ERROR: error in gc content flag!\n"
						printusage(sys.argv[0])
				else:
					print "ERROR: error in gc content flag!\n"
					printusage(sys.argv[0])
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: gc-content must be specified after -gc flag!\n"
				printusage(sys.argv[0])
	#weighted sampling, with -w flag
	for i in range(len(args)):
		if args[i]=='-w':
			weighted = 1
			del[args[i]]
			break
	#compute PSSM with -pssm flag
	for i in range(len(args)):
		if args[i]=='-pssm':
			if i+1!=len(args) and args[i+1][0]!='-':	
				PSSMfilename = args[i+1]
				del[args[i+1]]; del[args[i]]
				break
			else:
				print "ERROR: pssm file name must be specified after -pssm flag!\n"
				printusage(sys.argv[0])
	#verbose output, with -v flag
	for i in range(len(args)):
		if args[i]=='-v':
			VERBOSE = 1
			del[args[i]]
			break
	#chack if -p is defined
	if peptidefile==None:
		print "ERROR: peptide file must be specidfied using -p flag!\n"
		printusage(sys.argv[0])
		
	peptideD = parsePeptideFile(peptidefile)
	n = len(peptideD['+0'])
	if(constraintfile):
		constraint = parseConstraintFile(constraintfile,n)
	else:
		constraint = None
	main(peptideD,numSamples,weighted,gclow,gcup,BFS,PSSMfilename,filename,threshold,constraint)
