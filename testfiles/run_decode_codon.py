#!/usr/bin/env python

"""
Test script to decode codon DNA

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

from custarr_to_fastq import *
from parse_fastq import *

finput      = "miseqtest"
fqfile      = "cod_newtestfile.fastq"
trnsfile    = "cod_newtestfile_translated.txt"

print "Now converting CustomArray file to FASTQ format..."

convert(finput,fqfile)

print "FASTQ file written."

print "Now translating FASTQ file to human-readable..."

translate_codon(fqfile,trnsfile)

print "Translated files written to %s" % trnsfile
