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
	def init_mats():
		self.scores = np.empty(shape=(len(self.seq1),len(self.seq2)))
		self.traceback = np.empty(shape=(len(self.seq1),len(self.seq2)))
		self.gapX = np.empty(shape=(len(self.seq1),len(self.seq2)))
		self.gapY = np.empty(shape=(len(self.seq1),len(self.seq2)))
		# TODO: initialize()?

	def set_seqs(seq1, seq2):
		self.seq1 = read_sequence(seq1)
		self.seq2 = read_sequence(seq2)

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
		set_seqs(seq1, seq2)

		# initialize score matrix with seq1/seq2 
		init_mats()

		# loop through matrix, skipping boundary row/col
		for i in range(1,len(self.seq1)):
			for j in range(1,len(self.seq2)):
				fill_in(i,j)


		# traceback and return the score and alignment
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

		
	'''
	Fills in score in the scores matrix and arrow in the traceback
	for the current cell
	Parameters:
		i: row index
		j: col index
	Returns: score/arrows are filled in
	'''
	def fill_in(i, j, in_gap_x, in_gap_y):

		match = self.scores[i-1,j-1] + self.smat.loc[self.seq1[i], self.seq2[j]]

		gap_x = max(self.scores[i,j-1] + gap_open + gap_extend,
			self.gapY[i,j-1] + gap_extend)
		gapX[i,j] = gap_x

		gap_y = max(self.scores[i-1,j] + gap_open + gap_extend,
			self.gapX[i-1,j] + gap_extend)
		gapY[i,j] = gap_y

		score = max(match, gap_x, gap_y)

		# if there's a tie in gaps, just go with a gap in x
		if score == gap_x == gap_y:
			in_gap_x = True
		else if score == gap_x:
			in_gap_x = True
		else if score == gap_y:
			in_gap_y = True

		# reset in_gap_x or in_gap_y if gap is closed

		# Fill in a zero if best score is negative
		if score < 0:
			score = 0
		# update max if needed
		if score >= self.scores[self.max]:
			self.max = (i,j)

		# fill in arrow
		if score != 0:
			self.traceback[i,j] = 


class NeedlemanWunsch(PairwiseAligner):
	def initialize():
		# initialize first row and column with gap cost

	def fill_in(i,j):
