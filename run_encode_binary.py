#!/usr/bin/env python

from binaryDNA import *

infile = "infile.txt"
outfile = "dna_binary.txt"
check = "dna_binary_check.txt"

# Binary
bt2d = BinaryTextToDNA()
d2bt = DNAToBinaryText()

# Binary Text to DNA
bt2d.translate(infile, outfile)

# Check
d2bt.translate(outfile, check)
