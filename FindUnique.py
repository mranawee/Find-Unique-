##!/usr/bin/env python3
# Name: Mano Ranaweear
# Group Members: Ikenna Anigbogu (ianigbog), Andrew Gjelseen(agjelste)
'''
This is a python command line program that reads a FASTA file of 22 tRNA sequences.
The program must make a powerset of each tRNA and then be filtered to unique powersets 
in one tRNA that will not be found in another. Once all tRNA sets are completly unique,
the proram will filter down to what is essential.

Input: python findUnique.py < bos-tRNA.fa
'''

#from sequenceAnalysis import FastAreader
import sys
class FastAreader :

    

    def __init__ (self, fname=''):

        '''contructor: saves attribute fname '''

        self.fname = fname

            

    def doOpen (self):

        if self.fname is '':

            return sys.stdin

        else:

            return open(self.fname)

 

    def readFasta (self):

        

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

class FindUnique():
	'''
	Class contains a constructor to instantiate important objects, 
	find the unique subsets, and using the FASTAReader to print in the proper format.
	'''	
	powerSetList = []  #adding all powersets to this list
	def __init__(self, head, seq):
		'''
		Constructor method initiates the header, sequence, list of powerSets, and an
		iterator that ends up adding all 22 powerSets to a list.
		'''
		
		self.seq = seq.replace(".","").replace("_","").replace("-","")
		self.head = head

		self.powerSet = set()
		
		for start in range(len(self.seq)):
			for stop in range(start+1, len(self.seq)+1):
				subSeq = self.seq[start:stop]
				self.powerSet.add(subSeq) #powerset made(help by Mohammad Abdulqader)

		self.powerSetList.append(self.powerSet) #powersets added to list

		self.unique = set()

	def uniquesAndEssentials(self):

		'''
		For finding only the unique power sets and trimming down to the essential subset needed
		for the tRNA.  
		'''
		temp = self.powerSetList.copy()
		temp.remove(self.powerSet)
		unique = self.powerSet - set.union(*temp)
		'''
		For finding what's essential
		'''
		nonesntls = set()
		'''
		Algorithm adds base to the left and to the right,
		then adds to the nonesntls set.(Help By Dennis Mulligan)
		'''
		for sub in unique:
			rSub = sub[:-1]
			lSub = sub[1:]
			if rSub in unique:
				nonesntls.add(sub)
			elif lSub in unique:
				nonesntls.add(sub)

		esntls = unique - nonesntls
		esntlPosList = []
		'''
		For possible duplicate essentials
		'''
		for start in range(len(self.seq)):
			for stop in range(start+1, len(self.seq)+1):
				subSeq = self.seq[start:stop]
				if subSeq in esntls:
					esntlPosList.append((start, subSeq))
			
		return (esntlPosList)

def main():
	'''
	Main methods structures the output format, printing 
	the needed analysis for each tRNA.  Order is sorted
	by header alphabetically.
	'''
	
	'''
	This for loop is for feeding the FastAReader, and sorting by
	the header alphabetically.
	'''
	listOftRNAs = []
	for head, seq in FastAreader().readFasta():
		head = head.replace(' ','') 
		tRNA = FindUnique(head, seq)
		listOftRNAs.append(tRNA)
		
		finalListOftRNAs = sorted(listOftRNAs, key = lambda x:x.head)

	'''
	For each tRNA, the print format is arrnanged in this for loop and 
	sorted.
	'''
	for tRNA in finalListOftRNAs:
		listOfUniques = tRNA.uniquesAndEssentials()
		print(tRNA.head)
		print(tRNA.seq)
		
		final = sorted(listOfUniques, key = lambda x:x[0])
	
		'''
		The dots are multiplied by the index of what is essential 
		to output in this desired format. 
		'''
		for ix, esntl in final:
			final = ix * '.' + esntl 
			print(final)

if __name__ == "__main__":
	
	main()


