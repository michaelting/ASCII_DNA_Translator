#!/usr/bin/env python

"""
CustomArray to FASTQ converter

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

import os

def convert(infile, outfile):
	with open(infile,"r") as f:
		seqnum = 0
		out = open(outfile,"w")
		for line in f:
			out.write("@seq%d\n" % seqnum)
			out.write("%s" % line)
			out.write("+seq%d\n" % seqnum)
			out.write("%s" % line)
			out.flush()
			os.fsync(out.fileno())
			seqnum += 1
		out.close()
	print "Conversion compelte!"