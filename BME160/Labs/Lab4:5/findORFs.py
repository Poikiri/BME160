#!/usr/bin/env python3

"""
This script processes sequences from a FastA file to identify open reading frames (ORFs) and outputs the formatted results.

Modules:
    - sequenceAnalysis: Contains the NucParams class for nucleotide parameter analysis.
    - FastAreader: Contains the FastAreader class for reading FastA formatted files.
    - argparse: Standard library for parsing command-line arguments.

Classes:
    - CommandLine: Handles the command line, usage, and help requests using argparse.

Functions:
    - reverseComplement(sequence): Calculates and returns the reverse complement of a DNA sequence.
    - formatOutput(header, orfs): Formats the ORF output based on provided details.
    - main(inFile=None, options=None): Processes sequences from a FastA file, identifies open reading frames (ORFs), and outputs the formatted results.

Usage:
    The script can be run from the command line with the following options:
    - -lG, --longestGene: Only the longest gene.
    - -mG, --minGene: Minimum gene length (choices: 0, 100, 200, 300, 500, 1000).
    - -s, --start: Start codon(s) (default: 'ATG').
    - -t, --stop: Stop codon(s) (default: 'TAG', 'TGA', 'TAA').

Example:
    python script.py -lG -mG 300 -s ATG -t TAG -t TGA -t TAA input.fasta
    Input: fastA file
    Output: [frame] [start location ... stop location] [sequence length]
    
Details:
    The main function reads sequences from a FastA file and identifies ORFs by:
    - Iterating through each sequence.
    - Generating reverse complements for both strands.
    - Iterating through each strand and reading frame to find ORFs.
    - Using start and stop codons specified in command-line arguments.
    - Formatting and printing the identified ORFs.

The CommandLine class initializes an argparse parser to handle command-line options. 
The reverseComplement function computes the reverse complement of a given DNA sequence. 
The formatOutput function formats the identified ORFs for output. 
The main function orchestrates the reading, processing, and outputting of sequences and ORFs.
"""

from sequenceAnalysis import NucParams
from FastAreader import FastAreader
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
        '''Implement a parser to interpret the command line argv string using argparse.'''
        import argparse
        self.parser = argparse.ArgumentParser(description="Find Open Reading Frames (ORFs)")
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False, help='Only the longest gene')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(0, 100, 200, 300, 500, 1000), default=100, action='store', help='Minimum gene length')
        self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?', help='Start codon(s)')
        self.parser.add_argument('-t', '--stop', action='append', default=['TAG', 'TGA', 'TAA'], nargs='?', help='Stop codon(s)')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

def reverseComplement(sequence):
    """Calculate and return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'} # maps each base to its complement
    return ''.join(complement.get(base, base) for base in reversed(sequence)) # iterates through reversed sequence replacing using complement and joining the string together

def formatOutput(header, orfs):
    """Format the ORF output based on provided details."""
    output = [header] # initialize output list with the header
    sortedOrfs = sorted(orfs, key=lambda x: (-x[4], x[2])) # sort ORFs by length (descending) and start position (ascending)
    for orf in sortedOrfs:
        frame, start, end, length = orf[0]*orf[1], orf[2], orf[3], orf[4] # calculate frame and extract start, end, and length
        output.append(f"{frame:+d} {start:>5d}..{end:>5d} {length:>5d}")  #format ORF details into a string and append to output
    return "\n".join(output) # join all output lines into a single string with newline characters

def main(inFile = None, options = None):
    '''Processes sequences from a FastA file, identify open reading frames (ORFs), and output the formatted results.'''
    thisCommandLine = CommandLine(options)
    reader = FastAreader(inFile)
    
    for header, seq in reader.readFasta(): # For each sequence in a FastA file
        orfs = [] # list of ORFs
        seqLen = len(seq)
        reverseComplementSeq = reverseComplement(seq) # Reverse Compliment of sequence
        for strand in [1, -1]:
            sequence = seq if strand == 1 else reverseComplementSeq
            for frame in [0,1,2]:
                start = 0 # Adjust for longest gene
                for i in range(0, len(seq), 3): # Boundary Cond start
                    if sequence[i+frame:i+frame+3] in thisCommandLine.args.stop or i+3 >= len(sequence):
                        if i >= thisCommandLine.args.minGene:
                            orf = (strand, frame+1, 1, i + 3, i+3) if strand == 1 else (strand, frame+1, 1, i + 3, i+3)
                            orfs.append(orf) # add ORF to orfs list
                        if thisCommandLine.args.longestGene: # If we are including nested ORFs (longestGene == False) reset search pos to 0
                            start = i # Continue searching from start for start codons
                        break
                longPause = 0 # Variable to store last stop codon for longest gene
                for i in range(start, len(seq), 3):
                    if i > longPause and sequence[i+frame:i+frame+3] in thisCommandLine.args.start: # Begin search for stop codon if you reach a start codon
                        j = i # Begin searching for stop codon at start codon
                        while j <= len(sequence): # Continue through rest of sequence
                            if j+3 >= len(sequence) or sequence[j+frame:j+frame+3] in thisCommandLine.args.stop: # If end is reached (boundary condition) or if stop codon is reached
                                if (j-i) >= thisCommandLine.args.minGene: # Min Length Check
                                    orf = (strand, frame+1, i+1+frame, j+3+frame, j-i+3) if strand == 1 else (strand, frame+1, seqLen-j-frame-2, seqLen-i-frame, j-i+3)
                                    orfs.append(orf) # add ORF to orfs list
                                if thisCommandLine.args.longestGene: # If we are not including nested ORFs (longestGene == True) continue search after stop codon
                                    longPause = j # Continue searching after stop codon (ORF)
                                break # Exit loop as stop has been found
                            j+=3
        #print
        formatted_output = formatOutput(header, orfs)
        print(formatted_output)

if __name__ == "__main__":
    main()