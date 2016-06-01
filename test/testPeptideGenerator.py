#! /usr/bin/python2.6
#generate random peptides that code in all six reading frames
#The output of this script can be used as input peptide file to overcode.py
import random
from utility import translateMRNA
rnalength = 17
p=[]
while (len(p)!=1):
	s='' 
	for i in range(rnalength):
		s += random.choice('ACUG')
	#print s,s[:-2],s[1:-1],s[2:],s[:-2][::-1],s[1:-1][::-1],s[2:][::-1]
	p0 = translateMRNA(s[:-2])
	p1 = translateMRNA(s[1:-1])
	p2 = translateMRNA(s[2:])
	p_0 = translateMRNA(s[:-2][::-1])
	p_1 = translateMRNA(s[1:-1][::-1])
	p_2 = translateMRNA(s[2:][::-1])
	p = set([len(p0),len(p1),len(p2),len(p_0),len(p_1),len(p_2)])
print s
#print p0,p1,p2,p_0,p_1,p_2
f = open('p-test','w')
f.write("+0\t%s\n+1\t%s\n+2\t%s\n-0\t%s\n-1\t%s\n-2\t%s\n"%(p0,p1,p2,p_0,p_1,p_2))
f.close()
