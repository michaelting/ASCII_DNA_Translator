#!/usr/bin/env python

from os import fsync
from subprocess import call
from ASCIIcodons import *

infsta = "infile.fasta"
cfname = "chunker_codon.txt"
mfname = "chunker_rmfront_codon.txt"
trname = "chunker_trans_codon.txt"

print "Running arraychunker.py on %s" % infsta

call(['python','arraychunker.py',infsta,cfname,'76','76','TGAC'])

print "Array chunks placed in %s" % cfname

infile = open(cfname, "r")
newfile = open(mfname, "w")

for line in infile:
	newfile.write("%s\n" % line[22:])
	newfile.flush()
	fsync(newfile.fileno())

infile.close()
newfile.close()

print "Removed front adaptors from DNA, results in %s" % mfname

obj = DNAToText()
obj.translate(mfname, trname)

print "DNA translated in %s" % trname
