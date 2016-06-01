#! /usr/bin/env python
#script to test Overcode.py
#A. Bayegan
from overcode import *
from utility import aa2codon
import sys,os,itertools

def isCompatibleRNAwithIupacConstraint(rna,constraint):
	assert(len(rna)==len(constraint))
	iscompat= True
	for i in range(len(rna)):
		if rna[i] not in iupacNuc(constraint[i]):
			iscompat=False
	return iscompat

def generateAllCodingSeq(pep):
	l=[]
	for i in range(len(pep)):
		l.append(aa2codon[pep[i]])
	allSeq = list(itertools.product(*l))
	return allSeq

def naiveOverlapPf(p1,p2,fr,direc):
	assert len(p1)==len(p2)
	ss=[]
	n = 3*len(p1)
	l1 = generateAllCodingSeq(p1)
	l2 = generateAllCodingSeq(p2)
	for s1 in l1:
		for s2 in l2:
			seq1 = ''.join(s1)
			seq2 = ''.join(s2)
			#print seq1
			#print seq2
			if direc == '+':
				#print seq1[fr:], seq2[:n-fr]
				if seq1[fr:] == seq2[:n-fr]:
					#print "***",seq1,seq2
					ss.append(seq1+seq2[n-fr:])
	print(len(ss))
	#~ for s in ss:
		#~ print s
	return ss

def checkPSSM(PSSM):
	keys = PSSM.keys()
	keys.sort() # keys is list [(0,'A'),...,(3n,'U')] 
							# where n = min(len(peptides),len(peptide2))
	N = len(keys)/4
	for i in range(N):
		As = PSSM[(i,'A')]
		Cs = PSSM[(i,'C')]
		Gs = PSSM[(i,'G')]
		Us = PSSM[(i,'U')]
		print "%d:\t%.4f\t%.4f\t%.4f\t%.4f\n" % (i,As,Cs,Gs,Us)
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

#WARNIGN: This funtion cannot be used if iupac codes used  
def test(peptides,numSamples,weighted,gclow,gcup,BFS,PSSM,filename=None,threshold=1,constraint=None):
	seqLen = len(peptides['+0'])
	gc1 = int((3*seqLen+2)*gclow/100)
	gc2 = math.ceil((3*seqLen+2)*gcup/100)
	error=0
	#------------------check partition function and BFS-----------------
	ZB = computeZB(peptides,constraint,filename,threshold)
	ZB_gc = computeZB_GC(peptides,constraint,filename,threshold)
	ZF = computeZF(peptides,constraint,filename,threshold)
	ZF_gc = computeZF_GC(peptides,constraint,filename,threshold)
	
	ZFval  = 0.0
	for ch3 in NUCL:
		for ch4 in NUCL:
			ZFval += ZF[(n,ch3,ch4)]
	ZFval_gc  = 0.0
	for ch3 in NUCL:
		for ch4 in NUCL:
			for x in range(3*n+3):
				ZFval_gc += ZF_gc[(n,x,ch3,ch4)]	
	ZBval  = 0.0
	for ch1 in NUCL:
		for ch2 in NUCL:
			ZBval += ZB[(0,ch1,ch2)]			
	ZBval_gc  = 0.0
	for ch1 in NUCL:
		for ch2 in NUCL:
			for x in range(3*n+3):
				ZBval_gc += ZB_gc[(0,x,ch1,ch2)]
	
	rnaList = bfs(peptides,filename,threshold,constraint)			
	BFSval = len(rnaList)
	if not (ZFval==ZFval_gc==ZBval_gc==ZBval==BFSval):
		print ZFval,ZFval_gc,ZBval,ZBval_gc,BFSval
		print "something is wrong in partition function computation!"
		error=1
	
	#--------------------check sampling--------------------------------
	RNA    = sample(peptides,ZF,seqLen,numSamples,constraint,weighted,filename,threshold)
	for rna in RNA:
		#~ print rna
		#~ print translateMRNA(rna[:-2]),peptides['+0'],(translateMRNA(rna[:-2])==peptides['+0'] )
		#~ print translateMRNA(rna[1:-1])==peptides['+1']
		#~ print translateMRNA(rna[2:])==peptides['+2']
		#~ print translateMRNA(rna[:-2][::-1])==peptides['-0']
		#~ print translateMRNA(rna[1:-1][::-1])==peptides['-1']
		#~ print translateMRNA(rna[2:][::-1])==peptides['-2']
		if not(peptides['+0'] == translateMRNA(rna[:-2]) and peptides['+1'] == translateMRNA(rna[1:-1]) and 
			peptides['+2'] == translateMRNA(rna[2:]) and peptides['-0'] == translateMRNA(rna[:-2][::-1]) and
			peptides['-1'] == translateMRNA(rna[1:-1][::-1]) and peptides['-2'] == translateMRNA(rna[2:][::-1])):
				print "something is wrong in sampling without GC-content!"
				error=1
		if(constraint): assert(isCompatibleRNAwithIupacConstraint(rna,constraint))
	RNA = sample_GC(peptides,ZF_gc,seqLen,numSamples,gc1,gc2,constraint,weighted,filename,threshold)
	for rna in RNA:
		#~ print rna
		#~ print translateMRNA(rna[:-2]),peptides['+0'],(translateMRNA(rna[:-2])==peptides['+0'] )
		#~ print translateMRNA(rna[1:-1])==peptides['+1']
		#~ print translateMRNA(rna[2:])==peptides['+2']
		#~ print translateMRNA(rna[:-2][::-1])==peptides['-0']
		#~ print translateMRNA(rna[1:-1][::-1])==peptides['-1']
		#~ print translateMRNA(rna[2:][::-1])==peptides['-2']
		if not(peptides['+0'] == translateMRNA(rna[:-2]) and peptides['+1'] == translateMRNA(rna[1:-1]) and 
			peptides['+2'] == translateMRNA(rna[2:]) and peptides['-0'] == translateMRNA(rna[:-2][::-1]) and
			peptides['-1'] == translateMRNA(rna[1:-1][::-1]) and peptides['-2'] == translateMRNA(rna[2:][::-1])):
				print "something is wrong in sampling with GC-content!"
				error=1
		if not (gc1<=GCcont(rna)<= gc2):
			print "something is wrong in sampling with GC-content! GC-content of sampled rnas is not correct"
			error=1
		if(constraint): assert(isCompatibleRNAwithIupacConstraint(rna,constraint))
		
	#--------------------check PSSM ---------------------
	pssm = createPSSM(peptides,ZF,ZB,filename,threshold)
	#checkPSSM(pssm)
	RNA    = sample(peptides,ZF,seqLen,numSamples,constraint,weighted,filename,threshold)
	pssm2= computePSSMfromList(RNA)
	crossCheckPSSMvalues(pssm,pssm2)
	if error==0:
		print 'program seems correct!'

if __name__ == '__main__':
	if len(sys.argv) < 3:
		printusage(sys.argv[0])
	args     = sys.argv[1:]
	
	#------------------default setting------------------
	numSamples=10
	weighted = -1
	gclow=40;gcup=50
	BFS=0
	threshold=1
	PSSM=0
	filename=None
	constraintfile=None
	
	#------------------parse input arguments-------------
	
	#peptide file indicated by '-p' flag
	for i in range(len(args)):
		if args[i]=='-p':
			peptidefile = args[i+1]
			del[args[i+1]]; del[args[i]]
			break
		else:
			print "peptide file must be specified!\n"
			printusage(sys.argv[0])
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
			if args[i+1].isdigit():
				numSamples = int(args[i+1])
			elif args[i+1]=='all':
				BFS = 1
			else: 
				print "error in -n flage!\n"
				printusage(sys.argv[0])
			del[args[i+1]]; del[args[i]]
			break
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
			gcinput = args[i+1].split(':')
			if(len(gcinput)==1):
				gclow = int(gcinput[0])
				gcup = gclow
			elif(len(gcinput)==2):
				gclow = int(gcinput[0])
				gcup = int(gcinput[1])
				if gclow>gcup or gclow<0:
					print "error in gc content flag!"
					printusage(sys.argv[0])
			else:
				print "error in gc content flag!"
				printusage(sys.argv[0])
			del[args[i+1]]; del[args[i]]
			break
	#weighted sampling, with -w flag
	for i in range(len(args)):
		if args[i]=='-w':
			weighted = 1
			break
	#compute PSSM with -pssm flag
	for i in range(len(args)):
		if args[i]=='-pssm':
			PSSM = 1
			break
	#verbose output, with -v flag
	for i in range(len(args)):
		if args[i]=='-v':
			VERBOSE = 1
			break

	peptideD = parsePeptideFile(peptidefile)
	n = len(peptideD['+0'])
	if(constraintfile):
		constraint = parseConstraintFile(constraintfile,n)
	else:
		constraint = None
	test(peptideD,numSamples,weighted,gclow,gcup,BFS,PSSM,filename,threshold,constraint)


