from sequenceAnalysis import *
from FastAreader import *

def subMain (fileName=None):
    
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq) # calls addSequence method and has seq as an attribute

    totalNucleotides = myNuc.nucCount() # Calculate the sequence length
    sequenceLengthMb = totalNucleotides / 1_000_000  # Convert base count to megabases
    # calculates the GC content
    # has default val of 0 for G and C
    # if  no nucleotides then default val of 0 so no division by 0 error
    gcContent = ((myNuc.nucComposition().get('G', 0) + myNuc.nucComposition().get('C', 0)) / totalNucleotides) * 100 if totalNucleotides > 0 else 0 

    print(f'sequence length = {sequenceLengthMb:.2f} Mb\n')
    print(f'GC content = {gcContent:.1f}%\n')
    
    codonUsage = myNuc.codonComposition()  # Get codon counts
    aaComposition = myNuc.aaComposition()  # Get amino acid counts from codons

    sortedCodons = sorted(codonUsage.keys(), key=lambda x: (NucParams.rnaCodonTable[x], x)) # Sort codons by amino acid and then by codon (within the same AA group)

    # Print relative codon usage for each codon
    for codon in sortedCodons:
        aa = NucParams.rnaCodonTable[codon]  # Get amino acid for the codon
        count = codonUsage[codon]  # Count of this specific codon
        total_aaCount = aaComposition.get(aa, 0)  # Total counts of this amino acid, default val of 0
        frequency = (count / total_aaCount) * 100 if total_aaCount > 0 else 0 # if none of this AA default val of 0 so no division by 0 error
        print(f'{codon} : {aa} {frequency:5.1f} ({count:6d})')

    print('\n')

def main(): 
    
    genome1 = input() # get inputs
    genome2 = input()

    subMain(genome1) # call old main but now called submain
    subMain(genome2)

if __name__ == "__main__":
    main()