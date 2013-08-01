#!/usr/bin/env python

"""
Alternate DNA representation for 256-ASCII character encoding

Michael Ting
1 August 2013

DNA is a direct conversion of binary code to DNA bases:
    0: A or C
    1: G or T
    
DNA bases are randomly selected based on the binary character.
DNA homopolymers are limited to length 3.

DNA base assignments are chosen such that repetitive strings of
zero or one do not yield excess/deficiencies of:
    - purines/pyrimidines
        - Labeling 0 as A or G can result in excess purines for long strings of 0
    - specific base pairs
        - Labeling 0 as A or T can result in high GC content if more 1's
            are present than 0's
"""

import random

binzero2dna = { 2:'A',
               '2':'A',
                3:'C',
               '3':'C'}
    
binone2dna = { 2:'G',
              '2':'G',
               3:'T',
              '3':'T'}
            
bin2dna = {0:binzero2dna,
           '0':binzero2dna,
           1:binone2dna,
           '1':binone2dna}            
            
dna2bin = {'A':'0',
           'C':'0',
           'G':'1',
           'T':'1'}
           
class BinaryTextToDNA:

    def __init__(self):
        """
        Initialize BinaryTextToDNA object
        """
        print "Initialized BinaryTextToDNA object."

    def translate(self, infile, outfile):
        """
        Translates input text file to output DNA file
        """
        template = open(infile, 'rb')
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
        Converts a string of binary to DNA
        """
        
        binary = ''.join('{:08b}'.format(ord(c)) for c in txtstr)
        translated = ''
        for symbol in binary:
            choice = random.randint(2,3)    # determine which base to use        
            translated += bin2dna[symbol][choice]        
        return translated
    
    def translate_file_dna(self, infile):
        """
        Generator to yield translated lines from input file
        """
        for line in infile:
            yield self.text_to_dna(line.strip())

class DNAToBinaryText:

    def __init__(self):
        """
        Initialize DNAToBinaryText object
        """
        print "Initialized DNAToBinaryText object."
        
    def translate(self, infile, outfile):
        """
        Translates input DNA file to output text file
        """
        template = open(infile)
        newfile = open(outfile, "wb")
        
        print "Translating DNA to ASCII text..."    
        
        for line in self.translate_file_chr(template):
            newfile.write(line+"\n")
            print line
            
        template.close()
        newfile.close()
        
        print "File translation complete."
        print "See " + outfile + " for results."
    
    def dna_to_text(self, dnastr):
        """
        Translates a single string of DNA bases into ASCII text from binary
        """
        
        binary = ''
        for base in dnastr:
            binary += dna2bin[base]     
        txtstr = ''.join(chr(int(binary[i:i+8], 2)) for i in xrange(0, len(binary), 8))     
               
        return txtstr
    
    def translate_file_chr(self, infile):
        """
        Generator to yield lines of a DNA file translated into ASCII text
        """
        for line in infile:
            yield self.dna_to_text(line.strip())