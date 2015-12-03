#! /usr/bin/env python

# extractHIVgagPolJunction.py
# P.Clote

# Program extracts the gag and overlapping pol genes from HIV-1 genome (GenBank file)
# and translates the Gag and Pol proteins


import sys,os
from aminoAcidAndGeneticCodes import genCode,aaCode

def readingFrame(x):
  x = x%3
  if x==2:
    return -1
  return x

def translateMRNA(mrna):
  #WARNING: mrna is 0-indexed
  aaSeq = ""
  i     = 0
  while i<len(mrna)-1:
    codon  = mrna[i:i+3].replace('T','U')
    if len(codon)<3:
      print aaSeq
      print len(aaSeq)
      print "Untranslated nucleotides: %s" % codon
      sys.exit(1)
    gcode = genCode[codon].upper()
    #check for premature abortion of translation, which occurs in PR55
    #gag polyprotein
    if gcode == 'STOP':
      #print gcode,i,mrna[i:i+3]
      #print aaSeq
      return aaSeq
    aa     = aaCode[gcode]
    aaSeq += aa
    i     += 3
  return aaSeq

def readFastaNuclFile(filename):
  #WARNING: Add dummy character so string is 1-indexed
  file    = open(filename)
  genome  = "$"
  #Get FASTA comment
  line = file.readline().strip()
  if line[0]=='>':
    FASTA = line[1:].strip()
    line  = file.readline()
  else:
    FASTA = "Sequence"
  #Read genome
  while line:
    genome += line.strip().upper()    
    line    = file.readline()
  file.close()
  return genome

def aux(filename):
  genome          = readFastaNuclFile(filename)
  gagStartNt      = 336; gagStopNt       = 1838 
  polStartNt      = 1631; polStopNt       = 4642
  gagMRNA         = genome[gagStartNt:gagStopNt+1]
  polMRNA         = genome[polStartNt:polStopNt+1]
  gagProtein      = translateMRNA(gagMRNA)
  polReadingFrame = readingFrame(polStartNt-gagStartNt)
  polProtein      = translateMRNA(polMRNA)
  #print gagProtein
  #print polProtein
  #print polReadingFrame
  return genome,gagProtein,polProtein

  
def main(filename):
  genome,gagProtein,polProtein = aux(filename)
  ntStart = 1631
  ntStop  = 1682
  lenNuclSeq = ntStop-ntStart+1
  gagPartialProtein = translateMRNA(genome[1632:1839])
  polPartialProtein = translateMRNA(genome[1631:1839])
#
#  gagPartialProtein = translateMRNA(genome[1632:1683])
#  polPartialProtein = translateMRNA(genome[1631:1683])
#  text = "Gag amino acids in region %s..%s: " % (ntStart,ntStop)
  text = "Nucleotides in region %s..%s, a seq of len %s: " 
  print text % (ntStart,ntStop,lenNuclSeq)
  print genome[1631:1683].replace('T','U')
  text = "Gag amino acids in region %s..%s " % (1632,1682)
  print text,gagPartialProtein 
#  text = "Pol amino acids in region %s..%s" % (ntStart,ntStop)
  text = "Pol amino acids in region %s..%s" % (1631,1682)
  print text,polPartialProtein 
  cmd = "echo %s | RNAfold" % genome[1631:1683]
  secStr = os.popen(cmd).readlines()[1].split()[0]
  print len(secStr)
  print secStr  


if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s filename " % sys.argv[0]
    sys.exit(1)
  filename = sys.argv[1]
  main(filename)

