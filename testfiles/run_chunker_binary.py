#!/usr/bin/env python

"""
Converter script to produce CustomArray oligos for binary

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

from os import fsync
from subprocess import call
from binaryDNA import *

infsta = "infile.fasta"
cfname = "chunker_bin.txt"
mfname = "chunker_rmfront_bin.txt"
trname = "chunker_trans_bin.txt"

print "Running arraychunker.py on %s" % infsta

call(['python','binarraychunker.py',infsta,cfname,'48','48','ACGACTGT'])

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

obj = DNAToBinaryText()
obj.translate(mfname, trname)

print "DNA translated in %s" % trname
