#!/usr/bin/env python3
# Name: Jeevan Bhullar
# Group Members: Ariana, Adrian, Kelista, Sofia, & Nuohan

import sys, re

class NucParams:
    ''' Nucleotide string that returns metadata about codons and amino acids as relevant. Primarily helpful in working with a codon structured list. '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        ''' Starts building the codon and AA sequence. '''
        self.addSequence(inString)
        
    def addSequence (self, inSeq):
        ''' Given a string of nuc sequence add it to the front of the codon and AA lists. Does not validate but does upcase. '''
        if not hasattr(self, 'codonSequence'): # Create seq list if does not exist
            self.codonSequence = []
        if not hasattr(self, 'aaSequence'): # Create seq list if does not exist
            self.aaSequence = []
        codonTable = self.rnaCodonTable if 'U' in inSeq.upper() else self.dnaCodonTable # If contains U is RNA else DNA/Same
        newCodons = re.findall("(.{1,3})", inSeq.replace(" ", "").upper()) # Regex split into groups of 3 (codons), includes partial codon at end if relevant
        self.codonSequence = newCodons + self.codonSequence # Add New Codons to front of codonSequence list
        newAAs = []
        for codon in newCodons: # Iterate through all new codons and find matching AA
            newAAs += codonTable.get(codon, '')
        self.aaSequence = newAAs + self.aaSequence # Add New AAs to front of aaSequence list

    def getCodonString(self):
        ''' Returns string verson of the codon list (aka Nuc Sequence). '''
        return ''.join(self.codonSequence)

    def getAAString(self):
        ''' Returns string version of the AA list. '''
        return ''.join(self.aaSequence)

    def getValidCodons(self):
        ''' Returns the list of codons less any codons containing invalid nucleotides or shorter than 3 nucleotides. '''
        return [ codon for codon in self.codonSequence if sum([N in ['A', 'T', 'U', 'G', 'C', 'N'] for N in codon]) == 3 ]

    def aaComposition(self):
        ''' Returns dictionary keyed by single digit amino acids with counts of each AA coded by a codon in the sequence. '''
        aaComp = {} # Initialize blank dictionary
        for AA in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'W', '-']: # iterate through all AAs
            aaComp[AA] = self.aaSequence.count(AA) # Set dict value for each AA to count of AA in seq
        return aaComp
        
    def nucComposition(self):
        ''' Returns dictionary keyed by nucleotide with counts of each codon in the sequence. '''
        # This method returns a dictionary of counts of valid nucleotides found in the analysis. (ACGTNU}. If you were given RNA nucleotides, they should be counted as RNA nucleotides. If you were given DNA nucleotides, they should be counted as DNA nucleotides. Any N bases found should be counted also. Invalid bases are to be ignored in this dictionary.
        nucStr = self.getCodonString() # Single string of all nucleotides
        nucComp = {} # Initialize blank dictionary
        for N in ['A', 'T', 'U', 'G', 'C', 'N']: # iterate through all Nucleotide options
            nucComp[N] = nucStr.count(N)
        # Filter if nucComp[N] == 0 ?
        return nucComp
        
    def codonComposition(self):
        ''' Returns a dictionary keyed by rnaCodons with counts of the DNA and RNA nucleotides that match. '''
        codonList = [ codon for codon in list(set(self.getValidCodons())) if 'N' not in codon ] # Get valid (only approved NUCs) codons and filter out codons containing N
        codonComp = {} # Initialize blank dictionary
        for codon in codonList: # iterate through all codons
            codonComp[codon.replace('T', 'U')] = self.codonSequence.count(codon) # Set dict value for each codon (in rnaCodon format) to count of codon in seq
        return codonComp
        
    def nucCount(self):
        ''' Returns a simple sum of all counted nucleotides from the nucleotide composition dictionary. '''
        return sum(self.nucComposition().values())

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class ProteinParam (str):
    '''
    Class of a protein represented as a string of single letter Amino Acids (AAs).
    Contains methods to calculate metadata about the protein.
    Ignores invalid AA characters and normalizes to upper case amino acids.
    '''
    
    def __init__ (self, protein):
        ''' Creates an uppercase version of the protein '''
        self.upper = protein.upper() # Uppercase version of protein

    def aaCount (self): 
        ''' Calculates and returns count of amino acids in the protein '''
        return sum(self.aaComposition().values()) # Sum of counts of all AAs from aaComposition

    def aaComposition (self):
        ''' Counts AA composition, returning a dictionary keyed by single letter AA code with associated values as counts of AA in the protein. '''
        dictAA = {} # Initialize blank dictionary
        for AA in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'W']: # iterate through all AAs
            dictAA[AA] = self.upper.count(AA) # Set dict value for each AA to count of AA in protein
        return dictAA

    def _charge_ (self, pH):
        ''' Calculates and returns net charge on the protein at a specific pH '''
        aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6} # pKa of positively charged Amino Acids
        aaNterm = 9.69 # pKa of Nterm
        chargePos = ((10**aaNterm)/(10**aaNterm+10**pH)) # Charge of Nterm
        for AA, pKa in aa2chargePos.items(): # Factor all positive AA charges
            chargePos += self.aaComposition()[AA]*((10**pKa)/(10**pKa+10**pH))
        
        aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10} # pKa of negatively charged Amino Acids
        aaCterm = 2.34 # pKa of Cterm
        chargeNeg = ((10**pH)/(10**pH+10**aaCterm)) # Charge of Cterm
        for AA, pKa in aa2chargeNeg.items(): # Factor all negative AA charges
            chargeNeg += self.aaComposition()[AA]*((10**pH)/(10**pH+10**pKa))
        
        return chargePos-chargeNeg # Positive - Negative charges = Net Charge

    def pI (self, precision = 2):
        ''' Calculates and returns theoretical isolelectric point '''
        # Initial Values for Binary Search
        min = 0
        mid = 7
        max = 14
        search = True # Used to indicate required precision fidelity has not yet been reached

        while search: # While precision fidelity has not yet been reached, continue to search
            if abs(self._charge_((mid-min)/2+min))<abs(self._charge_((max-mid)/2+mid)): # Binary search branches left
                max = mid
                mid = (mid-min)/2+min
            else: # Binary search branches right
                min = mid
                mid = (max-mid)/2+mid
            if str(round(mid, precision)) == str(round(max, precision)): # If the mid and max match to required precision, the answer has been reached
                search = False
                
        return mid # min = mid = max at requested fidelity
        
        '''
        Old Code for "brute force" method
        pHMin = 999 # Default Start Value for pH with lowest charge
        minCharge = 999 # Default Start Value for lowest charge
        for pH in range(0, 1400): # range of 1400 is needed for 2 digit float fidelity
            charge = self._charge_(pH / 100.00) # test charge
            if abs(charge)<minCharge: # If test charge is lower than current minimum charge, save record
                # Update recorded minimums
                pHMin=pH / 100.00
                minCharge=abs(charge)
        
        # return when search has been completed for all 1400 values
        return pHMin
        '''

    def molarExtinction (self, Cystine=True):
        ''' Calculates and returns molarExtinction coefficent. How much light a protein absorbs at a certain wavelength '''
        aa2abs280 = {'Y':1490, 'W': 5500, 'C': 125} # absorbance at 280 nm
        if not Cystine: # If reducing conditions
            aa2abs280['C']=0 # Cystine is ignored (multipling by 0 ignores Cystine)
        molarExtinction = 0
        for AA, abs in aa2abs280.items():
            molarExtinction += self.aaComposition()[AA]*abs # Add product of AA count and AA abs to running sum of molarExtinction
        return molarExtinction
        pass

    def massExtinction (self, Cystine=True):
        ''' Calculates and returns massExtinction coefficent '''
        myMW =  self.molecularWeight() # Get MW for Extinsion calculation
        return self.molarExtinction(Cystine) / myMW if myMW else 0.0

    def molecularWeight (self):
        ''' Molecular weight (MW) of the protein sequence '''
        aa2mw = { # molecular weights of AAs
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
        mwH2O = 18.015 # mol. weight of H2O
        molecularWeight = mwH2O # Formula given includes one base mwH2O addition
        for AA, count in self.aaComposition().items(): # Iterate through all Amino Acids
            molecularWeight += (aa2mw[AA]-mwH2O)*count # Add the product of count and (AA mass less water mass) to running MW
        return molecularWeight