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