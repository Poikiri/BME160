class NucParams:
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
        self.nucComp = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0} # creates a dict of all nucs with values of 0
        self.codonComp = {codon: 0 for codon in self.rnaCodonTable.keys()} # creates a dict of all the codons used with their values set to 0
        self.aaComp = {aa: 0 for aa in set(self.rnaCodonTable.values())} # creates a dict of all the AAs used with all their values set to 0
        
    def addSequence (self, inSeq):
        ''' Adds new sequences, assuming they start in frame 1 '''
        cleanedSeq = inSeq.upper().replace(' ', '') # parses given sequence
        cleanedSeq = cleanedSeq.replace('T', 'U') # Convert DNA to RNA for uniform processing

        for nuc in cleanedSeq:
            if nuc in self.nucComp:
                self.nucComp[nuc] += 1 # counts amount of each nucleotide

        for i in range(0, len(cleanedSeq) - 2, 3): # loop always in multiples of 3
            codon = cleanedSeq[i:i+3] # gets a 3 char substring for a codon
            if 'N' not in codon and codon in self.codonComp: # makes sure codon is valid
                self.codonComp[codon] += 1 # counts amount of each codon
                aa = self.rnaCodonTable[codon] # translates codons to AAs
                self.aaComp[aa] += 1 # counts amount of each AA
                
    def aaComposition(self):
        ''' Returns the amino acid composition '''
        return self.aaComp
        
    def nucComposition(self):
        ''' Returns the nucleotide composition'''
        return self.nucComp
        
    def codonComposition(self):
        ''' Returns the codon composition'''
        return self.codonComp
        
    def nucCount(self):
        ''' Returns the total count of valid nucleotides '''
        return sum(self.nucComp[nuc] for nuc in 'ACGUN') # adds together all nucleotide amounts in dict
    


    
    def reverseComplement(sequence):
        """Calculate and return the reverse complement of a DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'} # maps each bas to its complement
        return ''.join(complement.get(base, base) for base in reversed(sequence)) # iterates through reversed sequence replacing using complement and joining the string together

    def findOrfs(self, sequence, startCodons, stopCodons, minLength, longestGene):
        """Identify and return a list of open reading frames (ORFs) in the given DNA sequence."""
        orfs = [] # empty list for storing ORFs

        def scanFrame(frame, strand, s):
            """Scan a specific reading frame on a given strand for ORFs."""
            startPositions = [] #empty list for storing start positions
            frameOffset = frame - 1
            for i in range(frameOffset, len(s) - 2, 3): 
                codon = s[i:i + 3] # ensures codons are 3 nucs long
                if codon in stopCodons: # finds stop codon
                    if startPositions:
                        start = startPositions[0]
                        orfLength = i + 3 - start
                        if orfLength >= minLength: # ensures ORF is longer than min length
                            orfs.append((strand, frame, start + 1, i + 3, orfLength)) # add ORF to orfs list
                        startPositions = [] if longestGene else startPositions[1:]
                elif codon in startCodons:
                    startPositions.append(i)

            if startPositions:  # check for open ended ORFs
                start = startPositions[0]
                orfLength = len(sequence) - start
                if orfLength >= minLength:
                    orfs.append((strand, frame, start + 1, len(sequence), orfLength))

        for frame in range(1, 3 + 1): # scans forward strand frames
            scanFrame(frame, '+', s=sequence)

        revSequence = NucParams.reverseComplement(sequence) # gets reverse of original sequence
        for frame in range(1, 3 + 1): # scans reverse strand frames
            scanFrame(frame, '-', s=revSequence)
        print(frame[0:9])
        print(revSequence[0:9])
        print(revSequence[0:9:-1])
        return orfs # returns complete list of all ORFs across all 6 frames