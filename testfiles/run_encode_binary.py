#!/usr/bin/env python

"""
Test script to encode binary DNA

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

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
