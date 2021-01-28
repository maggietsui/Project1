# Project 1 - Sequence Alignment
## Due 01/27/2021

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