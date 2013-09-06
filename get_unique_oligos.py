#!/usr/bin/env python

"""
Processes merged sequencing data from Illumina MiSeq

Created 15 July 2013
Updated 5 September 2013

Copyright 2013 Michael Ting
https://github.com/michaelting
Released under the BSD 2-clause license. See LICENSE.
http://opensource.org/licenses/BSD-2-Clause

Merged data was produced using SeqPrep on the forward and reverse
reads from Illumina MiSeq. Reads are merged together into single
reads for greater accuracy and are stripped of adapter sequences.
"""

import parse_fastq
import ASCIIcodons
import os, re, errno, gzip, time

class Counter(dict):
    """
    Counts the number of items placed into dictionary buckets.
    Values are initialized to 0.
    Counter class adapted from Pacman AI projects at UC Berkeley, CS188
    """
    def __getitem__(self, idx):
        self.setdefault(idx, 0)
        return dict.__getitem__(self, idx)
        
    def argMax(self):
        """
        Returns the key with the highest value
        """
        if len(self.keys()) == 0:
            return None
        allitems = self.items()
        values = [v[1] for v in allitems]
        maxIndex = values.index(max(values))
        return allitems[maxIndex][0]
        
    def maxVal(self):
        """
        Returns the maximum value of counts across all keys
        """
        if len(self.keys()) == 0:
            return None
        allitems = self.items()
        values = [v[1] for v in allitems]
        return max(values)
        
    def sortedKeys(self):
        """
        Returns a list of keys sorted by their values. 
        Keys with the highest values appear first.
        """
        sortedItems = self.items()
        compare = lambda x,y: sign(y[1] - x[1])
        sortedItems.sort(cmp=compare)
        return [x[0] for x in sortedItems]
        
    def totalCount(self):
        """
        Returns the sum of the counts for all keys.
        """
        return sum(self.values())
        
    def copy(self):
        """
        Returns a copy of the counter
        """
        return Counter(dict.copy(self))
        
def sign(x):
    """
    Returns 1 or -1 depending on the sign of x
    """
    if(x >= 0):
        return 1
    else:
        return -1
        
def sort_oligos(parser, tfunc):
    """
    Stores oligos by person ID and oligo ID in a nested dictionary structure
    using counts of bases for each base position in each oligo by means of
    a Counter object.
    {{[]}}
    """    
    
    TAGLEN = 28
    
    badcharpattern  = re.compile('[^ACGT]')
    tagpattern      = re.compile('TGTC[ACGT]{8}TGAT[ACGT]{12}')
    tagformat       = re.compile('#\d{2}\$\d{3}')
    
    #initialize RAM storage
    ramdict = {}
    
    # look at each read in the sequencing file
    for rec in parser:
        
        seqhead     = rec[0].strip()    # header
        seqdna      = rec[1].strip()    # dna sequence
        
        # if sequence contains non-ATGC characters, skip it
        findbadchars = badcharpattern.search(seqdna)
        
        if findbadchars:
            continue
        
        # search for starting point of information
        findtag = tagpattern.search(seqdna)

        # exclue bad tags
        if not findtag:
            continue
        
        # check that tag is translated into correct format
        tagcheck = tfunc(findtag.group())
        checkformat = tagformat.match(tagcheck)
        # exclude bad tags found in translation
        if not checkformat:
            continue
        infostart = findtag.start()
        
        # correct the reading frame
        seqdna = seqdna[infostart:]

        # cap maximum length
        if len(seqdna) > 104:
            seqdna = seqdna[:104]
        
        # if length is not divisible by 4, throw it out
        # occurs for sequences shorter than maximum length
        if (len(seqdna) % 4) != 0:
            continue
        
        # check the tag of the sequence       
        tagdna  = seqdna[:TAGLEN]
        msgdna  = seqdna[TAGLEN:]
        person  = tagdna[4:12]
        oligo   = tagdna[16:28]        
        
        pid     = tfunc(person)
        oid     = tfunc(oligo)
   
        # check if person has been listed
        if pid not in ramdict:
            ramdict[pid] = {}
        # check if oligo has been listed for person
        if oid not in ramdict[pid]:
            ramdict[pid][oid] = []
        
        for baseindex in range(len(msgdna)):
            
            base = msgdna[baseindex]
            # counter not yet initialized for that position
            if baseindex >= len(ramdict[pid][oid]):
                ramdict[pid][oid].insert(baseindex, Counter())  # Counter object holds counts of A,C,G,T
                
            ramdict[pid][oid][baseindex][base] += 1

    return ramdict

def get_consensus(ramdict):
    """
    Determines consensus sequences using the base with the highest count
    at each position.
    """

    COUNT_THRESHOLD = 100
    badoligos = []
    
    for pid in ramdict:
        for oid in ramdict[pid]:
            # determine the consensus sequence of a particular pid, oid
            consensus = ""
            for poscounter in ramdict[pid][oid]:
                # Counts below threshold imply erroneous sequences, since correct sequences
                # are copied 100's-1000's of times
                if poscounter.maxVal() < COUNT_THRESHOLD:
                    # throw out the sequence below the threshold
                    badpair = (pid, oid)
                    badoligos.append(badpair)
                    break
                # Grab the base with the most counts
                else:
                    consensus += poscounter.argMax()
            # replace the list of Counters with a consensus string
            ramdict[pid][oid] = consensus
            
    # remove bad oligos from the dictionary
    for pair in badoligos:
        badpers = pair[0]
        badolig = pair[1]
        del ramdict[badpers][badolig]            
            
    return ramdict

def condense(ramdict, tfunc, outdir):
    """
    Condense oligos into a single block of text based on pid,oid order, translates
    DNA to ASCII, and writes output to files in outdir
    """
    
    # create the output file directory
    try:
        os.makedirs(outdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise 
    
    for pid in ramdict:

        fullseq = ""        
        condensedfile = outdir + "/" + pid + "_condensed.txt"
        translatedfile = outdir + "/" + pid + "_translated.txt"
        
        old_oid = "000"        
        
        for oid in sorted(ramdict[pid]): # sorts from 000,001,002,003,...
        
            fullseq += ramdict[pid][oid] # concatenate DNA into one block
    
            curr = int(oid)
            old  = int(old_oid)
            # numerical skips imply erroneous sequences, since oligos have order
            if (curr - old) > 1:
                break
            old_oid = oid
        
        # don't write files with sequences filtered out
        if fullseq == "":
            continue            
            
        # write to condensed file and translate
        with open(condensedfile, "w") as condfile:
            condfile.write(fullseq)
        
        with open(translatedfile, "w") as endfile:
            with open(condensedfile, "r") as basefile:
                for line in basefile:
                    dna = line.strip()
                    text = tfunc(dna)
                    endfile.write(text+"\n")

def main():
    
    fqfile = gzip.open("merged.fastq.gz")
    parser = parse_fastq.readFastq(fqfile)      # faster with generator
    
    translator = ASCIIcodons.DNAToText()
    translate_dna = translator.dna_to_text
    
    treepath = "decodeddna"

    print "Now sorting oligos..."
    sort_start = time.time()
    rd = sort_oligos(parser, translate_dna)
    sort_end = time.time()
    
    print "Retrieving consensus DNA sequences..."
    comb_start = time.time()
    rd = get_consensus(rd)
    comb_end = time.time()
    
    print "Condensing DNA sequences and translating DNA to readable text..."
    cond_start = time.time()
    condense(rd, translate_dna, treepath)
    cond_end = time.time()
    
    print "Run complete. See /condensed.txt and /translated.txt in each subfolder of /%s for results." % treepath
    
    print "Elapsed time to sort was %g seconds" % (sort_end - sort_start)
    print "Elapsed time to combine was %g seconds" % (comb_end - comb_start)
    print "Elapsed time to condense was %g seconds" % (cond_end - cond_start)
            
if __name__ == "__main__":
    main()
