#!/usr/bin/env python3
# Name: Logan Roukis (lroukis)
# Group Members: Mikayla Lorico, Judah Tapert, Pietro Gazzola, Sibhi Sakthivel

''' Analyze a protein string and return needed information.
    No assumptions were made beyond the design parameter's scope.'''
  
class ProteinParam :
    ''' The protein string can have other characters within without affecting the output.
    aaCount() will count the valid AAs in the seq.
    aaComposition() will create a dict with the amount each AA.
    molecularWeight() adds up the molecular weights of all the AAs.
    _charge_() calculates the charge of the AAs for the pI method.
    pI() takes the charge and finds the theoretical pI.
    molarExtinction() finds the molar extinction coefficient.
    massExtinction() finds the mass extinction coefficient.
    
    input: VLSPADKTNVKAAW
    output: protein sequence? VLSPADKTNVKAAW
    Number of Amino Acids: 14
    Molecular Weight: 1499.7
    molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67
    Theoretical pI: 9.88
    Amino acid composition:
	A = 21.43%
	C = 0.00%
	D = 7.14%
	E = 0.00%
	F = 0.00%
	G = 0.00%
	H = 0.00%
	I = 0.00%
	K = 14.29%
	L = 7.14%
	M = 0.00%
	N = 7.14%
	P = 7.14%
	Q = 0.00%
	R = 0.00%
	S = 7.14%
	T = 7.14%
	V = 14.29%
	W = 7.14%
	Y = 0.00% '''
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }
    mwH2O = 18.015
    
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        self.protein = protein.upper()
        self.aaComposition = self.aaComposition() # bad for looking through later oops

    
    def aaCount (self):
        ''' Counts how many Amino Acids are in the protein sequence given through the input'''
        count = 0 # set a value for the count of AAs
        for aa in self.protein: # iterates through the protein sequence
            if aa in ProteinParam.aa2mw: # make sure to call the class so that it registers aa2mw
                count += 1 # if the character correlates to a valid AA in aa2mw, then we include it in the length
        return count # return the valid length of AAs

    
    def pI (self): 
        ''' Uses the charge to find the theoretical PI/ Isoelectric point'''
        pHVals = [i * 0.01 for i in range(0, 1401)]  # creates a list from 0 to 14.00 in increments of 0.01
        charges = [self._charge_(pH) for pH in pHVals] #calls the charge method and uses the charge to search through the list we just made
        # then finds the pH value where the net charge is closest to zero
        minCharge = min(charges, key=abs)  # find the minimum charge by absolute value, iterates the abs operation through the list
        closestIndex = charges.index(minCharge)  # find the index of this minimum charge
        return pHVals[closestIndex]

    
    def aaComposition (self):
        ''' Creates a dict containing all the AAs used and the amount of times they are used'''
        aaComp = {} # creates an empty dict for amino acids and their freq
        for aa in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'K', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'Y', 'W']:
            aaComp[aa] = self.protein.count(aa)  # goes through all the amino acids, if in sequence then .count will count them
        return aaComp

    
    def _charge_ (self, pH):
        ''' Calculates the total charge of the Amino Acids given the protein sequence'''
        charge = 0
        for aa, pKa in self.aa2chargePos.items(): # gets positive charge
            pKaPosVal = self.aaComposition.get(aa, 0)  # gets pKa value for the AA key, has default of 0
            charge += pKaPosVal * ((10 ** pKa) / (10 ** pKa + 10 ** pH)) # adds pos charge to netcharge
            
        for aa, pKa in self.aa2chargeNeg.items(): # gets negative charge
            pKaNegVal = self.aaComposition.get(aa, 0) # gets pKa value for the AA key, has default of 0
            charge -= pKaNegVal * ((10 ** pH) / (10 ** pH + 10 ** pKa)) # adds neg charge to netcharge
            
        charge += ((10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)) # get n terminus charge
        charge -= ((10 ** pH) / (10 ** pH + 10 ** self.aaCterm)) # get c terminus charge
        return charge

    
    def molarExtinction (self, Cystine = True):
        ''' Calculates molar extinction using given coefficient and given protein sequence'''
        mExtinction = 0
        for aa in self.aaComposition: 
            if aa in self.aa2abs280: # if the key/AA is in both dicts then 
                if aa == 'C' and not Cystine: # Skip cystine calculation under reducing conditions
                    continue  
                product = self.aaComposition[aa] * self.aa2abs280[aa] # amount of aa multiplied by coeeficients
                mExtinction += product # all products added together
        return mExtinction

    
    def massExtinction (self, Cystine = True):
        ''' Calculates mass extinction using molar extinction and molecular weight'''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    
    def molecularWeight (self):
        ''' Calculates the molecular weight of the given Amino Acids in the sequence'''
        mWeight = 0
        for aa in self.aaComposition: 
            if aa in self.aa2mw: # if the key/AA is in both dicts then 
                product = self.aaComposition[aa] * self.aa2mw[aa] # both of its values are multiplied
                mWeight += product # all of the AAs weights are added together
                mWeight -= (self.aaComposition[aa] * self.mwH2O) # subtracts the weight of one molecule of water per amino acid 
        mWeight = mWeight + self.mwH2O # there's one less peptide bond than amino acid in a sequence so add back one water molecule
        return mWeight # need to subtract H2O weight
        

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition.items(), # took out () after aaComposition as error said that 'dict' object was not callable
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()