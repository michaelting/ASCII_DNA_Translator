#!/usr/bin/python

"""
DNA representation of 256-ASCII Character Encoding

Michael Ting
15 July 2013
Updated 17 July 2013

DNA is represented alphabetically as base 4:
    A - 0
    C - 1
    G - 2
    T - 3
    
4-base codons will be used to represent the 256 Extended ASCII encodings

"""

import itertools

"""======================================================================="""

"""Variables for language translation"""

# TEMPORARY FIX, SHOULD FOCUS ON UNIFORM MAPPING OF CHARACTERS
"""
basevals = {'A':0,
            'C':1,
            'G':2,
            'T':3}
"""            
basevals = {'T':0,
            'A':1,
            'G':2,
            'C':3}
bases = ['A','C','G','T']
# retrieve all length 4 dna base strings
codons = []
CODONSIZE = 4
for lst in list(itertools.product(bases, repeat=CODONSIZE)):
    final = ''
    # using itertools.product, seqlst looks like 
    # [['A','A','A','A',],['A','A','A','T'],['A','A','A','G'],...] 
    # convert each sublist into a string like
    # ['AAAA','AAAT','AAAG',...]
    for base in lst:
        final += base
    codons.append(final)

"""
ASCII code to ASCII character
"""
# ASCII code to ASCII character
codes = range(0,256)
# map ASCII numerical code to ASCII character
# Ex. {65:'A'}
ind2chr = dict((c, chr(c)) for c in codes)
"""
ASCII character to ASCII code
"""
chr2ind = dict((v,k) for k,v in ind2chr.items())

# DNA to ASCII code in base 4
# ex. {'AAAC':1}

"""
DNA to ASCII code
DNA (Base 4)
"""
# map dna 4-base string to decimal values
dna2ind = dict((cstr,64*basevals[cstr[0]]+16*basevals[cstr[1]]+4*basevals[cstr[2]]+1*basevals[cstr[3]]) for cstr in codons)

"""
ASCII code to DNA
"""
ind2dna = dict((v,k) for k,v in dna2ind.items())

"""
DNA to ASCII character table
"""
dna2chr = {}
# map 4-base dna strings to ASCII character
# Ex. {'CAAC':'A','ATTC':'='}
for c in codons:
    index = dna2ind[c]
    char = ind2chr[index]
    dna2chr[c] = char

"""
ASCII character to DNA
"""
chr2dna = dict((v,k) for k,v in dna2chr.items())

"""========================================================================="""

class TextToDNA:

    def __init__(self):
        """
        Initialize TextToDNA object
        """
        print "Initialized TextToDNA object."

    def translate(self, infile, outfile):
        """
        Translates input text file to output DNA file
        """
        template = open(infile)
        newfile = open(outfile, "w")
        
        print "Translating ASCII text to DNA..."    
        
        for line in self.translate_file_dna(template):
            newfile.write(line+"\n")
            print line
            
        template.close()
        newfile.close()
        
        print "File translation complete."
        print "See " + outfile + " for results."

    def text_to_dna(self, txtstr):
        """
        Converts a string of ASCII characters to DNA
        """
        translated = ''
        for symbol in txtstr:
            translated += chr2dna[symbol]
        return translated
    
    def translate_file_dna(self, infile):
        """
        Generator to yield translated lines from input file
        """
        for line in infile:
            yield self.text_to_dna(line.strip())

class DNAToText:

    def __init__(self):
        """
        Initialize DNAToText object
        """
        print "Initialized DNAToText object."
        
    def translate(self, infile, outfile):
        """
        Translates input DNA file to output text file
        """
        template = open(infile)
        newfile = open(outfile, "w")
        
        print "Translating DNA to ASCII text..."    
        
        for line in self.translate_file_chr(template):
            newfile.write(line+"\n")
            print line
            
        template.close()
        newfile.close()
        
        print "File translation complete."
        print "See " + outfile + " for results."

    def chunkify(self, dna, n):
        """
        Split DNA into n-base codons
        """
        for i in xrange(0, len(dna), n):
            yield dna[i:i+n]
    
    def dna_to_text(self, dnastr):
        """
        Translates a single string of DNA codons into ASCII text
        """
        translated = ''
        for codon in self.chunkify(dnastr, CODONSIZE):
            translated += dna2chr[codon]
        return translated
    
    def translate_file_chr(self, infile):
        """
        Generator to yield lines of a DNA file translated into ASCII text
        """
        for line in infile:
            yield self.dna_to_text(line.strip())