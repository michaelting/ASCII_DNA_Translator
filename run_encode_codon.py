#!/usr/bin/env python

from ASCIIcodons import *

infile = "infile.txt"
outfile = "dna_codon.txt"
check = "dna_codon_check.txt"

# Codons
t2d = TextToDNA()
d2t = DNAToText()

#========#
# Codons #
#========#

# Text to DNA
t2d.translate(infile, outfile)

# Check encoding
d2t.translate(outfile, check)

