#!/usr/bin/env python3
# BME 160 Final Project
# Group Members: Jeevan Bhullar & Logan Roukis

'''
BLORAT is a program designed to align ORFs from provided DNA sequences to reference genomes, featuring an integrated ORF Finder and multi-genome compatability for ease of use and extensibility.

BLORAT: BLAT-like ORF record alignment tool
BLAT: BLAST-like alignment tool
BLAST: Basic Local Alignment Search Tool

BLAT (BLAST-Like Alignment Tool) is a bioinformatics software tool for comparing nucleotide or protein sequences against large genomic databases. Developed by Jim Kent from UCSC, it is designed as a speed optimized version of BLAT for his use cases, quickly finding sequences of 95% and greater similarity of length 40 bases or more. This makes it highly efficient for mapping DNA sequences to genomes or finding protein homologs. 

Integrating data from BLAT with findings from an Open Reading Frame (ORF) finder can be useful in genomic studies. ORF finders identify regions within a DNA or RNA sequence that could potentially encode proteins, indicating possible gene locations. However, determining the function or relevance of these ORFs and ensuring they map accurately to known genomic locations requires further validation. This is where BLAT comes in; it can align these identified ORFs against a genomic database, confirming their existence in the genome and helping annotate them based on similarity to known genes. This complementary approach enhances the identification and annotation of genes, improving our understanding of genomic structures and functions. Extending the lab 5 ORF finder to incorporate BLAT searches would allow for this extension, creating a final program lovingly named BLORAT. BLAT-like ORF record alignment tool

In conversations with Hiram M Clawson, Sr. Software Developer at the Genomics Institute of UCSC, there are significant numbers of genetic counselors who use the genome browser and its' tools including BLAT in their work to support patient healthcare outcomes. In fact, this usage became quite apparent during a recent outage that resulted in a strong outpouring of emails to his department about the impacts of the outage from genetic councelors. This program builds on that use case by making the local version of BLAT more easily usable and interfacable with common file formats and useful output formats. 

** NOTE: Due to system and timeline limitations this program was only tested on ARM bases macOS systems. While this should work on all systems further validation is needed due to the potential complecations of the BLAT subprocess. **
'''

from sequenceAnalysis import FastAreader, NucParams
from os import path
import subprocess, requests, sys, math
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(prog='BLORAT',
                                              description = '%(prog)s features an integrated ORF Finder and leverages BLAT to find alignments against a reference genome from a provided DNA sequence (or multiple). Tools like this provide value for genomics research and clinical genomics, like genetics counseling. ', 
                                             epilog = '%(prog)s is designed as a student project for BME 160 at UCSC', 
                                             add_help = True, # default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options]'
                                             )
        self.parser.add_argument('-in', '--inputFileName', action = 'store', nargs = '?', const = True, required=True, help = 'Source File Name with Extension')
        self.parser.add_argument('-out', '--outputFileName', action = 'store', nargs = '?', const = True, required=True, help = 'Output File Name without Extension')
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs = '?', const = True, default = False, help = 'longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type = int, choices = (0,100,200,300,500,1000), default = 100, action = 'store', help = 'minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'], nargs = '?', help = 'start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'], nargs = '?', help = 'stop Codon') #allows multiple list options
        self.parser.add_argument('-oS', '--optimizedSearch', action = 'store', nargs = '?', const = True, default = False, help = 'Optimized Search via use of an 11.occ file')
        self.parser.add_argument('-g', '--referenceGenome', action = 'store', nargs = '?', const = True, default = 'hg38', help = 'Reference Genome you want to BLAT against, hg38 is default, file in format genome.2bit must be present in program directory. Will attempt to locate and download file from UCSC''s DB using the standard abbreviation for search (ex. hg38 or mm10)')
        self.parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s 1.0')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class ORF():
    ''' Open Reading Frame (ORF) object that stores metadata about an ORF. Additionally allows for a ORF summaries to be returned with key metadata. '''
    def __init__(self, frame, name, startNuc, codons, seqNucLen):
        ''' Saves metadata to self. '''
        self.name = name # seq name
        self.codons = codons # List of codons in ORF
        self.codonCount = len(codons)
        self.nucleotides = ''.join(codons) # Nucleotide sequence in ORF as string
        self.nucelotideCount = len(self.nucleotides)
        self.frame = frame # +/- 1/2/3
        self.start = startNuc if frame > 0 else seqNucLen - (startNuc + self.nucelotideCount) + 2 # Start Nuc pos in seq
        self.stop = startNuc + self.nucelotideCount - 1 if frame > 0 else seqNucLen - startNuc + 1 # Stop Nuc pos in seq (Last nuc in stop codon)

    def fa(self):
        ''' Returns formatted string in fa format with dynamically generated head/description. '''
        return '>{:+d}.{:05d}..{:05d}\n{}\n'.format(
            self.frame, # +/- 1/2/3
            self.start, # Start Nuc pos in seq
            self.stop, # Stop Nuc pos in seq (Last nuc in stop codon)
            "".join(self.codons) # String of all nucleotides/codons without any spacing/delimiters
        )
    
    def intro(self):
        ''' Returns formatted string of selected metadata. '''
        return '{:+d} {:>5d}..{:>5d} {:>5d}'.format(
            self.frame, # +/- 1/2/3
            self.start, # Start Nuc pos in seq
            self.stop, # Stop Nuc pos in seq (Last nuc in stop codon)
            self.nucelotideCount
        )
    
    def intro_BLAT(self):
        ''' Returns formatted string of selected metadata, with specific reference to BLAT results '''
        # Potential extension to allow user specified print params or additional alaignment filtering or only top 1/n
        return '{:+d} {:>5d}..{:>5d} {:>5d}'.format(
            self.frame, # +/- 1/2/3
            self.start, # Start Nuc pos in seq
            self.stop, # Stop Nuc pos in seq (Last nuc in stop codon)
            self.nucelotideCount
        )

class alignment():
    ''' Open Reading Frame (ORF) object that stores metadata about an ORF. Additionally allows for a ORF summaries to be returned with key metadata. '''
    def __init__(self, colsPSL):
        ''' Saves metadata to self. '''
        self.matches = int(colsPSL[0]) # Number of bases that match that aren't repeats
        self.misMatches = int(colsPSL[1]) # Number of bases that don't match
        self.repMatches = int(colsPSL[2]) # Number of bases that match but are part of repeats
        self.nCount = int(colsPSL[3]) # Number of "N" bases
        self.qNumInsert = int(colsPSL[4]) # Number of inserts in query
        self.qBaseInsert = int(colsPSL[5]) # Number of bases inserted in query
        self.tNumInsert = int(colsPSL[6]) # Number of inserts in target
        self.tBaseInsert = int(colsPSL[7]) # Number of bases inserted in target
        self.strand = colsPSL[8] # "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
        self.qSize = int(colsPSL[10]) # Query sequence size.
        self.qStart = int(colsPSL[11]) # Alignment start position in query
        self.qEnd = int(colsPSL[12]) # Alignment end position in query
        self.tName = colsPSL[13] # Target sequence name
        self.tSize = int(colsPSL[14]) # Target sequence size
        self.tStart = int(colsPSL[15]) # Alignment start position in target
        self.tEnd = int(colsPSL[16]) # Alignment end position in target
        self.blockCount = int(colsPSL[17]) # Number of blocks in the alignment (a block contains no gaps)
        self.blockSizes = colsPSL[18] # Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
        self.qStarts = colsPSL[19] # Comma-separated list of starting positions of each block in query
        self.tStarts = colsPSL[20] # Comma-separated list of starting positions of each block in target
        '''
        Calculations relating to quality of alignment based on Jim Kent's methodology as integrated into webBLAT
        '''
        self.pslScore = int((self.matches + ( self.repMatches / 2) ) - self.misMatches - self.qNumInsert - self.tNumInsert) # A calculated metric for a PSL alignment that rewards matches and penalizes mismatches and insertions, reflecting the overall alignment quality.
        self.percentIdentity = (100.0 - self.calcMilliBad() * 0.1)/100 # The proportion of matched bases in a PSL alignment, expressed as a percentage (FP <= 1), indicating the similarity between the query and target sequences.
        
    def calcMilliBad(self):
        '''
        Adapted from Jim Kent's Perl script, calculate score and percent identity using his logic.
        https://genome-source.gi.ucsc.edu/gitlist/kent.git/raw/master/src/utils/pslScore/pslScore.pl
        
        Novel metric as a measure of alignment quality that quantifies the number of mismatches, insertions, and size differences between the query and target sequences, scaled in parts per thousand.
        Dependancy for percentIdentity
        '''
        milliBad = 0 # Default to 0 errors
        qAliSize = self.qEnd - self.qStart
        tAliSize = self.tEnd - self.tStart
        aliSize = tAliSize if (tAliSize < qAliSize) else qAliSize
        if (aliSize <= 0):
            return milliBad
        sizeDif = abs(qAliSize - tAliSize)
        insertFactor = self.qNumInsert + self.tNumInsert
        total = (self.matches + self.repMatches + self.misMatches)
        if (total != 0):
            roundAwayFromZero = 3*math.log(1+sizeDif)
            if (roundAwayFromZero < 0):
                roundAwayFromZero = int(roundAwayFromZero - 0.5)
            else:
                roundAwayFromZero = int(roundAwayFromZero + 0.5)
            milliBad = (1000 * (self.misMatches + insertFactor + roundAwayFromZero)) / total
        return milliBad
    
    def intro(self):
        ''' Returns formatted string of selected metadata. '''
        return '{:05d}..{:05d} {:2d} {:.2%} {}'.format(
            self.qStart, # Alignment Start Pos by Q
            self.qEnd, # Alignment End Pos by Q
            self.pslScore, # Number of matches less number of deviations, secondary indicator of quality
            self.percentIdentity, # % indicator of alignment accuracy
            self.tName # ChrNum
        )

def generateReverseCompliment(str):
    ''' Given a string of DNA nucleotides return the reverese compliment as a string. ''' 
    return "".join({'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.get(base) for base in str[::-1]) # Reverse string and swap bases

def main(options = None):
    '''
    Scans sequences in a FastA file for ORFs and returns all ORFs.
    Can scan multiple sequences through a single file and return ORFs seperated by head, assuming all heads are unique.

    Min Seq Nuc length, Start/Stop Codons, and Nested ORF search are all CL params.
    '''
    CL = CommandLine(options)
    reader = FastAreader(CL.args.inputFileName)
    
    # Preflight Checks
    # Is BLAT installed
    if not path.isfile('./blat'):
        # No BLAT, therefore alert user and stop program
        print('BLAT is not installed, this program depends on BLAT which may be downloaded from UCSC and must be placed in the same directory as this program')
        print('Go to this URL and find the BLAT program for your OS http://hgdownload.soe.ucsc.edu/admin/exe/')
        sys.exit() # Terminate Run
    # Is Reference Genome 2bit file present
    if not path.isfile('./{}.2bit'.format(CL.args.referenceGenome)):
        try: # Download Genome
            print("Attempting to download requested genome")
            r = requests.get('https://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.2bit'.format(CL.args.referenceGenome), stream = True) 
            with open("{}.2bit".format(CL.args.referenceGenome), "wb") as genome: 
                for chunk in r.iter_content(chunk_size=1024): 
                     # writing one chunk at a time to genome 2bit file 
                     if chunk: 
                         genome.write(chunk)
        except: # If Unable to Download Genome, notify user and stop program
            print('Could not find or download the requested genome, please manually provide it in the format of [genome].2bit for example hg38.2bit in the same directory as this program.')
            sys.exit() # Terminate Run
    # Is optimized search 11.ooc file present if needed
    if CL.args.optimizedSearch and not path.isfile('./hg38.2bit.11.ooc'):
        # Optimized Search 11.ooc file is needed but is not present, therefore generate 11.occ file
        subprocess.run(["./blat {0}.2bit /dev/null /dev/null -makeOoc={0}.2bit.11.ooc -repMatch=2253".format(CL.args.referenceGenome)], shell=True)
    
    # ORF Finder, carry over from Lab 5
    print("ORF Finder is Searching Through Provided Sequence(s)")
    ORFs = {} # Dict keyed by head of seq containing a list of ORF objects
    for head, seq in reader.readFasta(): # For each sequence in FastA file
        ORFs[head] = [] # Creates entry in dict, also clears previous results in event of head conflict
        seqReverseCompliment = generateReverseCompliment(seq) # Reverse Compliment of the seq as a string
        frames = { # Dictionary for 3 frames offset by 1 nuc each time and their 3 reverse compliment counterparts. Integer keys, NucParam object values.
            1: NucParams(seq), 
            3: NucParams('N'+seq), 
            2: NucParams('NN'+seq),
            -1: NucParams(seqReverseCompliment), 
            -3: NucParams('N'+seqReverseCompliment), 
            -2: NucParams('NN'+seqReverseCompliment)
        }
        for frame, sequence in frames.items(): # Iterate through all frames
            codonSequence = sequence.codonSequence # Seq for internal manipulation
            if abs(frame) > 1: # Adjust first codon to remove leading Ns used to offset frames
                codonSequence[0] = codonSequence[0][-abs(frame)+4:]
            seqNucLen = len(''.join(codonSequence)) # seqLen in Nucleotides
            firstStop = False # firstStop is used to identify leading boundary condition
            pos = 0 # Search Position in sequence in codon #
            while pos < len(codonSequence): # Iterate through all codons in seq
                if not firstStop: # firstStop has not yet been found (boundary condition)
                    if codonSequence[pos] in CL.args.stop or pos+1 == len(codonSequence): # firstStop has been found or length has been reached (boundary condition)
                        if pos>=CL.args.minGene/3: # Min Length Check
                            ORFs[head] = ORFs[head] + [ORF(frame, head, 1, codonSequence[0:pos+1], seqNucLen)] # Add to DB
                        firstStop = True
                        if not CL.args.longestGene: # If we are including nested ORFs (longestGene == False) reset search pos to 0
                            pos = 0 # Continue searching from start for start codons
                elif codonSequence[pos] in CL.args.start: # Begin search for stop codon if you reach a start codon
                    pos2 = pos # Begin searching for stop codon at start codon
                    while pos2 <= len(codonSequence): # Continue through rest of sequence
                        if pos2 == len(codonSequence) or codonSequence[pos2] in CL.args.stop: # If end is reached (boundary condition) or if stop codon is reached
                            if (pos2-pos+1)>=CL.args.minGene/3: # Min Length Check
                                ORFs[head] = ORFs[head] + [ORF(frame, head, len(''.join(codonSequence[0:pos]))+1, codonSequence[pos:pos2+1], seqNucLen)] # Add to DB
                            if CL.args.longestGene: # If we are not including nested ORFs (longestGene == True) continue search after stop codon
                                pos = pos2 # Continue searching after stop codon (ORF)
                            break # Exit loop as stop has been found
                        pos2+=1
                pos+=1
            # frame search complete
        # sequence search complete
    # file complete
    print('Found {0} ORFs in {1} sequences'.format(
        sum(len(x) for x in ORFs.values()), # Total ORF count across all sequences
        len(ORFs) # Number of sequences.
    ))
    
    # Write ORFs to file in fastA format
    fileName = CL.args.outputFileName # Root outfile name is saved in varible for readability
    fa = open("{}.fa".format(fileName), "w") # Create/overwrite fa file of ORFs
    for name, sequence in ORFs.items(): # Iterate through all sequences in ORF dict
        sequence.sort(key=lambda x: (x.nucelotideCount, -x.start, x.frame>0, -abs(x.frame)), reverse=True) # Sort ORFs in head by length(h2l), start pos(l2h), frame side (+>-), frame number(l2h)
        for singleORF in sequence: # For each ORF
            fa.write(singleORF.fa()) # Print formatted FA seq
    fa.close() # Close File per best practices
    
    # Run BLAT Search
    '''
    genome 2bit file is used with ORFs.fa to create ORFs.psl file. 
    stepSize, repMatch, minScore, and minIdentity are defaults used in Web BLAT to mimic results.
    noHead creates output without header rows for better machine readability.
    {2} refers to a turnary that may add the occ parameter with generated 11.occ file for search speed based on user provided parameter
    
    Refer to BLAT help for further information on parameters
    https://genome.ucsc.edu/goldenPath/help/blatSpec.html#blatUsage
    '''
    print("Running BLAT Search, this may take a moment. Please Wait.")
    result = subprocess.run(["./blat {0}.2bit {1}.fa {1}.psl -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 -noHead {2}".format(
        CL.args.referenceGenome, # Reference genome without 2bit extension
        fileName, # Outfile name used for intermediate files
        '' if not CL.args.optimizedSearch else '-ooc={}.2bit.11.ooc'.format(CL.args.referenceGenome) # 11.ooc tells the program to load over-occurring 11-mers from an external file. This will increase the speed by a factor of 40 in many cases.
    )], shell=True, capture_output=True, text=True)
    print("Search Completed!")
    print(result.stdout) # Print loading and search details
    
    # Populate Alignments Dictionary
    alignments = {} # Dict keyed by ORF of seq containing a list of alignment objects
    psl = open('{0}.psl'.format(fileName), 'r').read().split('\n')
    for i in range(0, len(psl)-1): # Iterate Through all lines of the BLAT psl output (less the final trailing blank line).
        psl[i] = psl[i].split('\t') # Split line by tab to cols (TSV)
        name = psl[i][9] # FastA Head Name is saved in varible for readability
        if name not in alignments: # If no alignments exist for head
            alignments[name] = [] # Create Blank List
        alignments[name] = alignments[name] + [alignment(psl[i])] # Append to List
    
    # Write output to file in spec format
    aln = open("{}.txt".format(fileName), "w") # Create/overwrite txt file of alignments
    for name, sequence in alignments.items(): # Iterate through all sequences in alignment dict
        aln.write(name+'\n') # Print Head
        sequence.sort(key=lambda x: (x.percentIdentity, x.qStart)) # Sort alignments in head by accuracy(h2l), length(h2l), start pos(l2h)
        for singleAlignment in sequence: # For each alignment
            aln.write(singleAlignment.intro()+'\n') # Print formatted intro
        aln.write('\n') # Extra newline for readability across ORFs
    aln.close() # Close File per best practices
        
    print('Program Completed! Thank You!')

if __name__ == "__main__":
    main() 