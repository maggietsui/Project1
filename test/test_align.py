import pytest
from align import algs

@pytest.fixture
def some_relevant_data():
	return np.ones(10)

def test_fasta_io():
    file = "/Users/mtsui1/Documents/Classes/Algs/Project1/sequences/prot-0026.fa"
    seq = read_sequence(file)
    assert seq

def test_scoring_matrix_io():
	m = read_score_matrix("/Users/mtsui1/Documents/Classes/Algs/Project1/scoring_matrices/BLOSUM50.mat")
    assert m.shape == (24, 24)
    assert not m.empty

def test_identical():
	assert True

def test_alignment_score():
	assert True
