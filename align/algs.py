import numpy as np
import pandas as pd

class PairwiseAligner:
	def __init__(self, smat)
		self.smat = read_score_matrix(smat)
		self.scores = None
		self.gapX = None
		self.gapY = None
		self.traceback = None
		self.score = 0

	'''
	Method that reads in a protein sequence from a fasta file
	Parameters:
		file: Path to fasta file
	Returns: list containing the sequence from the file
	'''
	def read_sequence(file):
	    f = open(file,'r')
	    lines = f.readlines()[1:] # skip header
	    seq = []
	    for line in lines:
	        seq += line.strip()
	    return seq

	'''
	 Method that reads in score matrix 
	 Parameters:
		smat: path to score matrix
	 Returns: score matrix as a pandas df, where score lookup can be done
	 	as mat.loc["A","N"]
	'''
	def read_score_matrix(smat):
	    skiplines = 0
	    with open(smat,'r') as fh:
	        for line in fh:
	            if line.startswith("#") | line.startswith("."):
	                skiplines += 1
	            else:
	                break
	    scores = pd.read_csv(smat, delim_whitespace=True, skiprows=skiplines)
	    scores.set_index(scores.columns, inplace=True)
	    return scores

	'''
	Initialize matrices needed for alignment
	Parameters:
	Returns: score, gap X, gap Y, and traceback matrices
	'''
	def init_mats(seq1, seq2):
		self.scores = np.empty(shape=(len(seq1),len(seq2)))
		self.traceback = np.empty(shape=(len(seq1),len(seq2)))
		self.gapX = np.empty(shape=(len(seq1),len(seq2)))
		self.gapY = np.empty(shape=(len(seq1),len(seq2)))
		# TODO: initialize()?

	'''
	 Align method: 
	 Parameters:
	 	seq1: path to first sequence
	 	seq2: path to second sequence
		with_score: if true, returns the score of the alignment
			in addition to the alignments
		readable: if true, returns the two sequences with lines
			for readability
	 Returns: returns seq1 and seq2 with best alignment
	'''
	def align(self, seq1, seq2, with_score=False, readable=False)
		seq1 = read_sequence(seq1)
		seq2 = read_sequence(seq2)

		# initialize score matrix with seq1/seq2 
		init_mats(seq1, seq2)

		num_rows, num_cols = scores.shape
		for i in range(num_rows):
			for j in range(num_cols):


		# for 
class SmithWaterman(PairwiseAligner):
	self.max = tuple(0,0)

	'''
	Method to initialize the score-keeping and traceback matrices.
	seq1 is initialized on the rows, while seq2 on the columns
	Returns: two mat
	'''
	def initialize():
		# initialize first row and column with zeroes
		scores[0,:] = 0
		scores[:,0] = 0

		

	def get_score():

		# Fill in a zero if best score is negative
		# update max if needed

class NeedlemanWunsch(PairwiseAligner):
	def initialize()
		# initialize first row and column with gap cost

	def get_score():
