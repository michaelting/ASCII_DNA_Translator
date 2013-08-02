#!/usr/bin/python

"""
Splits a FASTA file of DNA into pieces for array orders

***REQUIRED***

- Run with contig and step size as 48
- TAG at head of dna is formate $ _ _ # _ _ _
  meaning 7 characters, so 56 DNA bases
- To ensure mod 8, pad the ends of the DNA with 2 bases (106 --> 104, padded with AA at the end)
- 104-56 = 48 for the chunk size
- stuffer sequence is TGAC, corresponding to single quote '

Steps to replicate:
$ python arraychunker.py outfile.txt 48 48 ACGACTGT
$ python
>>> infile = open("outfile.txt","r")
>>> newfile = open("framecorrect.txt","w")
>>> for line in infile:
...     newfile.write("%s\n" % line[22:])    # to account for correct reading frame since adaptor is 22 bp
$ ls
$ python
>>> from binaryDNA import *
>>> o = DNAToBinaryText()
>>> o.translate("framecorrect.txt","h_readable.txt")


"""

import os
from argparse import ArgumentParser
 
def process_file(infile, stepsize, chunksize, stuffer):
    """
    Processes a FASTA file by outputting the padded sequences line by line
    Input:
       infile - FASTA file of sequence(s)
    Output:
       List of padded sequences
    """
    seqs2order = []
    pid = 0 # index for person
    for name, seq in process_seq(infile):
        chunklist = get_chunks(seq, stepsize, chunksize, stuffer)        
        stuffed = stuff_ends(chunklist, chunksize, pid)
        for order in stuffed:
            seqs2order.append(order)
        pid += 1
    return seqs2order

def process_seq(infile):
    """
    Generator that finds all sequences in a FASTA file
    Input:
        infile - FASTA file of sequence(s)
    Output:
        name - the name of the sequence following the ">"
        seq - the nucleotide sequence
    """
    name, seq = None, []
    for line in infile:
        if line.startswith(">"):
            if name: 
                yield (name, ''.join(seq))
            name, seq = line.strip(), []
        else:
            seq.append(line.strip())
    if name: 
        yield (name, ''.join(seq))
  
def rev_comp(seq):
    """
    Calculates the reverse complement of a single nucleotide sequence
    Input:
        seq     - (String) nucleotide sequence of the input
    Output:
        rcomp   - (String) The reverse complement of the input seq
    """
    # DNA complementary base pairs
    base_dict = {'A':'T',
                 'a':'t',
                 'T':'A',
                 't':'a',
                 'G':'C',
                 'g':'c',
                 'C':'G',
                 'c':'g'}
    
    rev_seq = seq[::-1]
    
    rcomp = ""
    for base in rev_seq:
        rcomp = rcomp + base_dict.get(base)
    
    return rcomp
   
def get_chunks(seq, stepsize, chunksize, stuffer):
    """
    Splits up DNA sequence into array chunks
    """
    
    #STEPSIZE = 76
    #CHUNKSIZE = 106    
    SEQLEN = len(seq)
    #STUFFER = "TGAC"
    #00100111 == ACGACTGT
    
    clst = []    
    
    index = 0
    while index < SEQLEN:
        # last few chunks are smaller than chunksize
        if index+chunksize > SEQLEN:
            # take the remaining piece at the end
            chunk = seq[index:]
            # add a stuffer sequence
            remaining = chunksize - len(chunk)
            stuffnum = int(remaining / 8)   # the number of stuffer sequences to add
            tail = remaining % 8            # how much of the stuffer to add at the end
            chunk += stuffer*stuffnum + stuffer[:tail]
        else:    
            chunk = seq[index:index+chunksize]

        
        clst.append(chunk)
        index += stepsize
    
    return clst
    
def stuff_ends(clst, chunksize, pid):
    """
    
    CHUNKSIZE = 76    
    
    nth segment = foo[n*76:n*76+106]
    
    universalA = "CTACACGACGCTCTTCCGATCT"
    universalB = "TGCTGAACCGCTCTTCCGATCT"
    universalA + foo[n*76:n*76+106] + rc(universalB)
    
    "CTACACGACGCTCTTCCGATCT" + foo[n*76:n*76+106] + "AGATCGGAAGAGCGGTTCAGCA"    
    """
    
    # Indexing for DNA
    # has the format
    # $ _ _ # _ _ _
    # $ PID # CID    
    
    num2dna =   {'0': 'ACGTACAC',
                 '1': 'ACGTACAG',
                 '2': 'ACGTACGA',
                 '3': 'ACGTACTG',
                 '4': 'ACGTAGAC',
                 '5': 'ACGTATAG',
                 '6': 'ACGTATGA',
                 '7': 'ACGTATGT',
                 '8': 'ACGTTCAC',
                 '9': 'ACGTTACG'}
                 
    pstart = 'ACTACATG' # # (hash) 00100011
    cstart = 'ACTACGCA' # $ (dollar) 00100100
    
    universalA = "CTACACGACGCTCTTCCGATCT"
    universalB = "TGCTGAACCGCTCTTCCGATCT"
    
    RC_universalB = "AGATCGGAAGAGCGGTTCAGCA"
    check = rev_comp(universalB)

    CONTIGSIZE = len(universalA) + len(RC_universalB) + chunksize    

    # pad pid with zeroes if single digit
    if pid < 10:
        pid = '0'+str(pid)

    pidstr = ''        
    for digit in str(pid):
        pidstr += num2dna[digit]

    if RC_universalB != check:
        raise IOError("Reverse complement not calculating correctly")
    
    # add universals to the end of the dna region
    stuffedlist = []
    contigid = 0        # index for contig in assembly
    for chunk in clst:
        contstr = str(contigid)
        if contigid < 10:
            contstr = '00' + contstr
        elif contigid < 100 and contigid >= 10:
            contstr = '0' + contstr

        contdna = ''            
        for digit in contstr:
            contdna += num2dna[digit]
        
        # add index tag to piece
        tag = pstart + pidstr + cstart + contdna
        # create padded dna
        padded = universalA + tag + chunk + RC_universalB + 'AA' # AA needed when using coding length %8 = 0
        #if len(padded) != CONTIGSIZE:
        #    raise IOError("Chunk size incorrect!")
        stuffedlist.append(padded)
        contigid += 1
        
    return stuffedlist
    
def main():
    
    parser = ArgumentParser()
    
    parser.add_argument("infile",metavar="in",help="input FASTA file with sequence(s)")
    parser.add_argument("outfile",metavar="out",help="name of output file")
    parser.add_argument("stepsize",metavar="step",help="stepsize for chunk overlap")
    parser.add_argument("chunksize",metavar="chunk",help="length of chunks in bp")
    parser.add_argument("stuffer",metavar="stuffer",help="stuffer sequence")
    
    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    stepsize = int(args.stepsize)
    chunksize = int(args.chunksize)
    stuffer = str(args.stuffer)

    if stepsize > chunksize:
        raise IOError("Stepsize must be smaller than chunk size to allow for overlap")
        
    template = open(infile, "r")
    newfile = open(outfile, "w")
    
    orderlist = process_file(template, stepsize, chunksize, stuffer)
    
    for dna in orderlist:
        newfile.write("%s\n"% dna)
        newfile.flush()
        os.fsync(newfile.fileno())
        
    template.close()
    newfile.close()
    
    print "Chunker complete!"
    
    
if __name__ == "__main__":
    main()