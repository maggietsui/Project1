# Project 1 - Sequence Alignment

![BuildStatus](https://github.com/maggietsui/Project1/workflows/HW1/badge.svg?event=push)

## Documentation

### PairwiseAligner Class
A parent class used to align two sequences

#### read_sequence(self,file)
Method that reads in a protein sequence from a fasta file  
Parameters:  
	* file: Path to fasta file  
Returns: list containing the sequence from the file, in uppercase

#### read_score_matrix(self,smat)
Method that reads in score matrix  
Parameters:  
	* smat: path to score matrix  
Returns: score matrix as a pandas df, where score lookup can be done as mat.loc["A","N"]

#### get_seq1(self)
Getter method for seq1

#### get_seq2(self)
Getter method for seq2

#### get_score_mat(self)
Getter method for the current reference scoring matrix (BLOSUM50, etc)

#### set_score_mat(self, smat)
Setter method for the reference scoring matrix (BLOSUM50, etc)  
Parameters:  
	* smat: path to score matrix

#### set_gap_open(self, gap_open)
Setter method for gap opening cost. Takes in a negative value  
Parameters:  
	* gap_open: gap opening cost 

#### set_gap_extend(self, gap_extend)
Setter method for gap extend cost. Takes in a negative value  
Parameters:  
	* gap_extend: gap extend cost

#### set_seqs(self, seq1, seq2)
Setter method for the two sequences to be aligned  
Parameters:  
	* seq1: path to seq1  
	* seq2: path to seq2

#### align(self, seq1, seq2, return_alignment = True)
Method that aligns two sequences  
Parameters:  
* seq1: path to first sequence  
* seq2: path to second sequence  
* return_alignment: if true, returns the two aligned sequences with gaps if necessary  
Returns: Best alignment score and resulting aligned sequences



### SmithWaterman class
Class that inherits from PairwiseAligner, implementing the Smith-Waterman alignment method

#### init_mats(self)
Method to initialize the score-keeping and traceback matrices.  
seq1 is initialized on the rows, while seq2 on the columns.  
Matrices are of size M+1 x N+1, where M and N are the lengths  
of seq1 and seq2  
Returns: None, sets the instance's score, gap X, gap Y, and traceback matrices

#### fill_in(self, i, j)
Method that fills in score in the scores matrix and arrow in the traceback for the current cell  
Parameters:  
* i: current row index  
* j: current col index  
Returns: None, score/arrows are filled in

#### do_traceback(self, return_alignment)
Method that uses the traceback matrix to get the optimal alignment score and reconstruct the alignment  
Parameters:  
* return_alignment: whether or not to return the two  
* aligned sequences  
Returns: The max score in the matrix and the aligned sequences


### NeedlemanWunsch class
Class that inherits from PairwiseAligner, implementing the Needleman-Wunsch alignment method

#### init_mats(self)
Method to initialize the score-keeping and traceback matrices.  
seq1 is initialized on the rows, while seq2 on the columns.  
Matrices are of size M+1 x N+1, where M and N are the lengths  
of seq1 and seq2  
Returns: None, sets the instance's score, gap X, gap Y, and traceback matrices

#### fill_in(self, i, j)
Method that fills in score in the scores matrix and arrow in the traceback for the current cell  
Parameters:  
* i: current row index  
* j: current col index  
Returns: None, score/arrows are filled in

#### do_traceback(self, return_alignment)
Method that uses the traceback matrix to get the score and reconstruct the alignment  
Parameters:  
* return_alignment: whether or not to return the two  
* aligned sequences  
Returns: The max score in the bottom right cell of the matrix and the aligned sequences
