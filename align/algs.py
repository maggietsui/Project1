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


    def set_seqs(self, seq1, seq2):
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
        #return traceback(with_score, readable)

class SmithWaterman(PairwiseAligner):
    def __init__(self, smat):
        super().__init__(smat)
        self.max = (0,0)

    '''
    Method to initialize the score-keeping and traceback matrices.
    seq1 is initialized on the rows, while seq2 on the columns
    Parameters:
    Returns: score, gap X, gap Y, and traceback matrices
    '''
    def init_mats(self):
        self.scores = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.traceback = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.gapX = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.gapY = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        
        self.scores[:,0] = 0
        self.scores[0,:] = 0
        
        self.gapX[:,0] = 0
        self.gapX[0,:] = 0
        
        self.gapY[:,0] = 0
        self.gapY[0,:] = 0
        # TODO: initialize()?


    '''
    Fills in score in the scores matrix and arrow in the traceback
    for the current cell
    Parameters:
        i: row index
        j: col index
    Returns: score/arrows are filled in
    '''
    def fill_in(self, i, j):
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
            # for max score, get the cell that it came from
            # to resolve ties, a match is preferred, then X gap, then Y gap
            max_mats = [k for k,v in score_dict.items() if v == max_score]
            if "M" in max_mats:
                if max_score == match:
                    arrow = (i-1,j-1,"M")
                elif max_score == end_gapx:
                    arrow = (i-1,j-1,"Ix")
                elif max_score == end_gapy:
                    arrow = (i-1,j-1,"Iy")
            elif "Ix" in max_mats:
                if max_score == gapx_open:
                    arrow = (i,j-1,"M")
                else:
                    arrow = (i,j-1,"Ix")
            elif "Iy" in max_mats:
                if max_score == gapy_open:
                    arrow = (i-1,j,"M")
                else:
                    arrow = (i-1,j,"Iy")



            self.traceback[i,j] = arrow

        # update max if needed
        if max_score >= self.scores[self.max]:
            self.max = (i,j)


    def traceback():
        # a match
        if self.traceback[i,j] == "match":
            self.scores[i-1,j-1]
        if self.traceback[i,j] == "endgapx":
            self.gapX[i-1,j-1]
        if self.traceback[i,j] == "endgapy":
            self.gapY[i-1,j-1]
        if self.traceback[i,j] == "opengapx":
            self.scores[i,j-1]
        if self.traceback[i,j] == "extendgapx":
            self.gapX[i,j-1]
        if self.traceback[i,j] == "opengapy":
            self.scores[i-1,j]
        if self.traceback[i,j] == "extendgapy":
            self.gapY[i-1,j]
            
        #score = 0
        self.score = score
        #return score


class NeedlemanWunsch(PairwiseAligner):
    def __init__(self, smat):
        super().__init__(smat)
        
    def initialize():
        pass
        # initialize first row and column with gap cost

    def fill_in(i,j):
        pass