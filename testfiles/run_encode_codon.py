#!/usr/bin/env python

"""
Test script to encode codon DNA

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

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

