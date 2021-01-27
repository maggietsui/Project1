import numpy as np
import pandas as pd

class PairwiseAligner:
    def __init__(self, smat):
        self.smat = self.read_score_matrix(smat)
        self.scores = None
        self.gapX = None
        self.gapY = None
        self.traceback = None
        self.score = 0
        self.gap_open = -5
        self.gap_extend = -1
        self.seq1 = None
        self.seq2 = None

    '''
    Method that reads in a protein sequence from a fasta file
    Parameters:
        file: Path to fasta file
    Returns: list containing the sequence from the file, in uppercase
    '''
    def read_sequence(self,file):
        f = open(file,'r')
        lines = f.readlines()[1:] # skip header
        seq = ""
        for line in lines:
            seq += line.strip()
        return seq.upper()

    '''
     Method that reads in score matrix 
     Parameters:
        smat: path to score matrix
     Returns: score matrix as a pandas df, where score lookup can be done
        as mat.loc["A","N"]
    '''
    def read_score_matrix(self,smat):
        skiplines = 0
        with open(smat,'r') as fh:
            for line in fh:
                # Skip header lines
                if line.startswith("#") | line.startswith("."):
                    skiplines += 1
                else:
                    break
        scores = pd.read_csv(smat, delim_whitespace=True, skiprows=skiplines)
        scores.set_index(scores.columns, inplace=True)
        return scores

    def get_seq1(self):
        return self.seq1
    
    def get_seq2(self):
        return self.seq2
    
    def get_score_mat(self):
        return self.smat
    
    def set_score_mat(self, smat):
        self.smat = self.read_score_matrix(smat)

    def set_gap_open(self, gap_open):
        self.gap_open = gap_open
    
    def set_gap_extend(self, gap_extend):
        self.gap_extend = gap_extend
        
    def set_seqs(self, seq1, seq2):
        self.seq1 = self.read_sequence(seq1)
        self.seq2 = self.read_sequence(seq2)

    '''
     Method that aligns two sequences
     Parameters:
        seq1: path to first sequence
        seq2: path to second sequence
        return_alignment: if true, returns the two aligned sequences with 
        gaps if necessary
     Returns: Best alignment score and resulting aligned sequences
    '''
    def align(self, seq1, seq2, return_alignment = True):
        # Exit if empty sequence(s)
        if len(seq1) == 0 | len(seq2) == 0:
            raise ValueError("sequences cannot be empty")

        self.set_seqs(seq1, seq2)

        # initialize matrices needed for alignment
        self.init_mats()

        # loop through matrix, skipping boundary row/col
        # and fill in each cell
        for i in range(1,len(self.seq1)+1):
            for j in range(1,len(self.seq2)+1):
                self.fill_in(i,j)
        
        # traceback and return the best score and alignment
        if return_alignment == True:
            score, alignment = self.do_traceback(return_alignment)
            return score, alignment
        else:
            return self.do_traceback(return_alignment)

class SmithWaterman(PairwiseAligner):
    def __init__(self, smat):
        super().__init__(smat)
        self.max = (0,0)

    '''
    Method to initialize the score-keeping and traceback matrices.
    seq1 is initialized on the rows, while seq2 on the columns.
    Matrices are of size M+1 x N+1, where M and N are the lengths
    of seq1 and seq2
    Returns: None, sets the instance's score, gap X, gap Y, and 
    traceback matrices
    '''
    def init_mats(self):
        # Create empty numpy arrays
        self.scores = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.traceback = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1),dtype=object)
        self.gapX = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.gapY = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        
        # Set the first row and column
        self.scores[:,0] = 0
        self.scores[0,:] = 0
        
        self.gapX[:,0] = 0
        self.gapX[0,:] = 0
        
        self.gapY[:,0] = 0
        self.gapY[0,:] = 0


    '''
    Method that fills in score in the scores matrix and arrow in the traceback
    for the current cell
    Parameters:
        i: current row index
        j: current col index
    Returns: None, score/arrows are filled in
    '''
    def fill_in(self, i, j):
        # calculate match score
        residue_score = self.smat.loc[self.seq1[i-1], self.seq2[j-1]]
        match = self.scores[i-1,j-1] + residue_score
        end_gapx = self.gapX[i-1,j-1] + residue_score
        end_gapy = self.gapY[i-1,j-1] + residue_score

        match_score = max(match, end_gapx, end_gapy)

        # calculate gap in X score
        gapx_open = self.scores[i-1,j] + self.gap_open + self.gap_extend
        gapx_extend = self.gapX[i-1,j] + self.gap_extend
        gapx_score = max(gapx_open, gapx_extend)

        # calculate gap in Y score
        gapy_open = self.scores[i,j-1] + self.gap_open + self.gap_extend
        gapy_extend = self.gapY[i,j-1] + self.gap_extend
        gapy_score = max(gapy_open, gapy_extend)

        score_dict = {"M":match_score, "Ix":gapx_score, "Iy":gapy_score}

        # set negative scores to 0
        for key in score_dict.keys():
            if score_dict[key] < 0:
                score_dict[key] = 0

        # fill in the matrices with the corresponding value
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
                    self.traceback[i,j] = "match"
                elif max_score == end_gapx:
                    self.traceback[i,j] = "endgapx"
                elif max_score == end_gapy:
                    self.traceback[i,j] = "endgapy"
            elif "Ix" in max_mats:
                if max_score == gapx_open:
                    self.traceback[i,j] = "opengapx"
                else:
                    self.traceback[i,j] = "extendgapx"
            elif "Iy" in max_mats:
                if max_score == gapy_open:
                    self.traceback[i,j] = "opengapy"
                else:
                    self.traceback[i,j] = "extendgapy"

        # update max if needed
        if max_score >= self.scores[self.max]:
            self.max = (i,j)

    '''
    Method that uses the traceback matrix to get the optimal 
    alignment score and reconstruct the alignment
    Parameters:
        return_alignment: whether or not to return the two
        aligned sequences
    Returns: The max score in the matrix and the aligned sequences
    '''
    def do_traceback(self, return_alignment):
        # Start at the cell with the max score
        i = self.max[0]
        j = self.max[1]

        # reset max
        self.max = (0,0)

        fin_score = self.scores[i, j]
        
        if return_alignment == True:
            # start traceback at max value
            #current_val = self.scores[self.max[0], self.max[1]]
            str1 = ""
            str2 = ""
            # Stop when the current cell is None, which means the score was 0
            while self.traceback[i,j] != None:
                if self.traceback[i,j] == "match" or self.traceback[i,j] == "endgapx" or self.traceback[i,j] == "endgapy":
                    str1 = self.seq1[i-1] + str1
                    str2 = self.seq2[j-1] + str2
                    i = i-1
                    j = j-1
                # Add a gap in x direction (seq2)
                elif self.traceback[i,j] == "opengapx" or self.traceback[i,j] == "extendgapx":
                    str2 = "-" + str2
                    str1 = self.seq1[i-1] + str1
                    i = i-1
                # Add a gap in y direction (seq1)
                elif self.traceback[i,j] == "opengapy" or self.traceback[i,j] == "extendgapy":
                    str2 = self.seq2[j-1] + str1
                    str1 = "-" + str1
                    j = j-1

            alignment = str1 + '\n' + str2
            return fin_score, alignment
        return fin_score

class NeedlemanWunsch(PairwiseAligner):
    def __init__(self, smat):
        super().__init__(smat)
        
    def init_mats(self):
        self.scores = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.traceback = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1),dtype=object)
        self.gapX = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        self.gapY = np.empty(shape=(len(self.seq1)+1,len(self.seq2)+1))
        
        # initialize the first rows/cols of the matrices to boundary values
        self.scores[0,0] = 0
        self.scores[0,1:] = float("-inf")
        self.scores[1:,0] = float("-inf")
        
        self.gapX[:,0] = list(range(len(self.gapX[:,0])))
        self.gapX[:,0] = [idx*self.gap_extend+self.gap_open for idx in self.gapX[:,0]]
        self.gapX[0,1:] = float("-inf")
        
        self.gapY[1:,0] = float("-inf")
        self.gapY[0,:] = list(range(len(self.gapY[0,:])))
        self.gapY[0,:] = [idx*self.gap_extend+self.gap_open for idx in self.gapY[0,:]]
        
        self.traceback[0,1:] = "extendgapy"
        self.traceback[1:,0] = "extendgapx"

    '''
    Fills in score in the scores matrix and arrow in the traceback
    for the current cell
    Parameters:
        i: row index
        j: col index
    Returns: score/arrows are filled in
    '''
    def fill_in(self, i, j):
        # calculate match score
        residue_score = self.smat.loc[self.seq1[i-1], self.seq2[j-1]]
        match = self.scores[i-1,j-1] + residue_score
        end_gapx = self.gapX[i-1,j-1] + residue_score
        end_gapy = self.gapY[i-1,j-1] + residue_score

        match_score = max(match, end_gapx, end_gapy)

        # calculate gap in X score
        gapx_open = self.scores[i-1,j] + self.gap_open + self.gap_extend
        gapx_extend = self.gapX[i-1,j] + self.gap_extend
        gapx_score = max(gapx_open, gapx_extend)

        # calculate gap in Y score
        gapy_open = self.scores[i,j-1] + self.gap_open + self.gap_extend
        gapy_extend = self.gapY[i,j-1] + self.gap_extend
        gapy_score = max(gapy_open, gapy_extend)

        score_dict = {"M":match_score, "Ix":gapx_score, "Iy":gapy_score}


        self.scores[i,j] = score_dict["M"]
        self.gapX[i,j] = score_dict["Ix"]
        self.gapY[i,j] = score_dict["Iy"]

        max_score = max(score_dict.values())

        # assign arrows
        # for max score, get the cell that it came from
        # to resolve ties, a match is preferred, then X gap, then Y gap
        max_mats = [k for k,v in score_dict.items() if v == max_score]
        if "M" in max_mats:
            if max_score == match:
                self.traceback[i,j] = "match"
            elif max_score == end_gapx:
                self.traceback[i,j] = "endgapx"
            elif max_score == end_gapy:
                self.traceback[i,j] = "endgapy"
        elif "Ix" in max_mats:
            if max_score == gapx_open:
                self.traceback[i,j] = "opengapx"
            else:
                self.traceback[i,j] = "extendgapx"
        elif "Iy" in max_mats:
            if max_score == gapy_open:
                self.traceback[i,j] = "opengapy"
            else:
                self.traceback[i,j] = "extendgapy"

    def do_traceback(self, return_alignment):
        # start in the bottom right cell
        i = len(self.scores[:,0]) - 1
        j = len(self.scores[0,:]) - 1
        
        # score is the biggest value 
        fin_score = max(self.scores[i,j], self.gapX[i,j], self.gapY[i,j])

        if return_alignment == True:
            str1 = ""
            str2 = ""

            # traceback until you reach (0,0)
            while i > 0 or j > 0:
                if self.traceback[i,j] == "match" or self.traceback[i,j] == "endgapx" or self.traceback[i,j] == "endgapy":
                    str1 = self.seq1[i-1] + str1
                    str2 = self.seq2[j-1] + str2
                    i = i-1
                    j = j-1
                elif self.traceback[i,j] == "opengapx" or self.traceback[i,j] == "extendgapx":
                    str2 = "-" + str2
                    str1 = self.seq1[i-1] + str1
                    i = i-1
                elif self.traceback[i,j] == "opengapy" or self.traceback[i,j] == "extendgapy":
                    str2 = self.seq2[j-1] + str1
                    str1 = "-" + str1
                    j = j-1
        
            alignment = str1 + '\n' + str2
            return fin_score, alignment
        return fin_score