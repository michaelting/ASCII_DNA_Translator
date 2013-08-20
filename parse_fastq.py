#!/usr/bin/env python

"""
FASTQ file parser

Created 15 July 2013
Updated 20 August 2013

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause
"""

import os, sys, gzip
from itertools import ifilter, islice
from ASCIIcodons import *
from binaryDNA import *

def translate_bin(infile, outfile):
	"""Translates the pure sequence binary encoding information 
	from a FASTQ file into human-readable text"""
	parser = ParseFASTQ(infile)
	bintranslator = DNAToBinaryText()
	out = open(outfile,"w")
	for rec in parser:
		header = rec[0]
		seq = rec[1]
		text = bintranslator.dna_to_text(seq)
		out.write("%s\n" % header)
		out.write("%s\n" % text)
		out.flush()
		os.fsync(out.fileno())
		print header
		print text
	out.close()

def translate_codon(infile, outfile):
	"""Translates the codon encoding sequence from a FASTQ file into
	human-readable text"""
	parser = ParseFASTQ(infile)
	codtranslator = DNAToText()
	out = open(outfile,"w")
	for rec in parser:
		header 	= rec[0]
		seq 	= rec[1]
		text	= codtranslator.dna_to_text(seq)
		out.write("%s\n" % header)
		out.write("%s\n" % text)
		out.flush()
		os.fsync(out.fileno())
		print header
		print text
	out.close()

def readFastq(fastqfile):
    fastqiter = ifilter(lambda l: l, fastqfile)  # skip blank lines
    fastqiter = (l.strip('\n') for l in fastqiter)  # strip trailing newlines
    while True:
        values = list(islice(fastqiter, 4))
        if len(values) == 4:
            header1,seq,header2,qual = values
        elif len(values) == 0:
            raise StopIteration
        else:
            raise EOFError("Failed to parse four lines from fastq file!")

        if header1.startswith('@') and header2.startswith('+'):
            yield header1[1:], seq, qual
        else:
            raise ValueError("Invalid header lines: %s and %s" % (header1, header2))

class ParseFASTQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()
	By Augustine Dunn"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
		