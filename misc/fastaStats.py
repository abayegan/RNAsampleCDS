#! /usr/bin/env python
#A. Bayegan

import sys,os,numpy


def parseFasta(fastafile):
	d={}
	fin = open(fastafile,'r')
	line = fin.readline().strip()
	if not line:
		print "Error in fasta file!"
		sys.exit(1)
	seqid = line[1:]
	line = fin.readline().strip()
	seq=''
	while line:
		if line[0]!='>':
			seq += line
		elif line[0]=='>':
			if seqid in d:
				print 'sequences in the fasta file must have unique ids!'
				sys.exit(1)
			d[seqid]=seq.upper().replace('T','U')
			seq=''
			seqid = line[1:]
		line = fin.readline().strip()
	d[seqid] = seq.upper().replace('T','U')
	return d

if __name__ == '__main__':
	if len(sys.argv)!=2:
		print "usage: %s fastaFile"%sys.argv[0]
	d = parseFasta(sys.argv[1])
	l=[]
	for i in d.values():
		l.append(len(i))
	print "number of sequences: %d" %len(d.keys())
	print "shortest seq length: %f" % numpy.min(l)
	print "longest seq length: %f" % numpy.max(l)
	print "mean seq length: %f" % numpy.mean(l)
	print "seq length standard deviation: %f"% numpy.std(l)
