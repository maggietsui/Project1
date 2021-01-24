import numpy as np
import pandas as pd

class PairwiseAligner:
    def __init__(self, smat):
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
    def align(self, seq1, seq2, with_score=False, readable=False):
        if len(seq1) == 0 | len(seq2) == 0:
            raise ValueError("sequences cannot be empty")

        set_seqs(seq1, seq2)

        # initialize score matrix with seq1/seq2 
        init_mats()

        # loop through matrix, skipping boundary row/col
        for i in range(1,len(self.seq1)):
            for j in range(1,len(self.seq2)):
                fill_in(i,j)

        # traceback and return the score and alignment
        return traceback(with_score, readable)

class SmithWaterman(PairwiseAligner):
    def __init__(self, smat):
        super().__init__(self, smat)
        self.max = (0,0)

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
    def fill_in(i, j):
        gap_open = -5
        gap_extend = -1

        # calculate match score
        residue_score = self.smat.loc[self.seq1[i], self.seq2[j]]
        match = self.scores[i-1,j-1] + residue_score
        end_gapx = self.gapX[i-1,j-1] + residue_score
        end_gapy = self.gapY[i-1,j-1] + residue_score

        match_score = max(match, end_gapx, end_gapy)

        # calculate gap in X score
        gapx_open = self.scores[i,j-1] + gap_open + gap_extend
        gapx_extend = self.gapX[i,j-1] + gap_extend
        gapx_score = max(gapx_open, gapx_extend)

        # calculate gap in Y score
        gapy_open = self.scores[i-1,j] + gap_open + gap_extend
        gapy_extend = self.gapY[i-1,j] + gap_extend
        gapy_score = max(gapy_open, gapy_extend)

        score_dict = {"M":match_score, "Ix":gapx_score, "Iy":gapy_score}

        # set negative scores to 0
        for key in score_dict.keys():
            if score_dict[key] < 0:
                score_dict[key] = 0
        self.scores[i,j] = score_dict["M"]
        self.gapX[i,j] = score_dict["Ix"]
        self.gapY[i,j] = score_dict["Iy"]

        max_score = max(score_dict.values())

        # assign arrows (nonzero scores only)
        if max_score != 0:
            arrows = {}
            # for max score (including ties), get the cell that it came from
            max_mats = [k for k,v in score_dict.items() if v == max_score]
            for mat in max_mats:
                if mat == "M":
                    if max_score == match:
                        arrows[(i-1,j-1)] = "M"
                    elif max_score == end_gapx:
                        arrows[(i-1,j-1)] = "Ix"
                    elif max_score == end_gapy:
                        arrows[(i-1,j-1)] = "Iy"
                if mat == "Ix":
                    if max_score == gapx_open:
                        arrows[(i,j-1)] = "M"
                    else:
                        arrows[(i,j-1)] = "Ix"
                if mat == "Iy":
                    if max_score == gapy_open:
                        arrows[(i-1,j)] = "M"
                    else:
                        arrows[(i-1,j)] = "Iy"



            self.traceback[i,j] = arrows

        # update max if needed
        if max_score >= self.scores[self.max]:
            self.max = (i,j)


    def traceback():
        pass
        #score = 0
        #return score


class NeedlemanWunsch(PairwiseAligner):
    def __init__(self, smat):
        super().__init__(self, smat)
        
    def initialize():
        pass
        # initialize first row and column with gap cost

    def fill_in(i,j):
        pass