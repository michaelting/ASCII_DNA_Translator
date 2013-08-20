#!/usr/bin/env python

"""
Test script to decode binary DNA

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

from custarr_to_fastq import *
from parse_fastq import *

finput      = "binout_48_48_ACGACTGT_rmfradapt"
fqfile      = "bin_newtestfile.fastq"
trnsfile    = "bin_newtestfile_translated.txt"

print "Now converting CustomArray file to FASTQ format..."

convert(finput,fqfile)

print "FASTQ file written."

print "Now translating FASTQ file to human-readable..."

translate_bin(fqfile,trnsfile)

print "Translated files written to %s" % trnsfile
