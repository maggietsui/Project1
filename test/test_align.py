import pytest
from align.algs import *

@pytest.fixture
def test_fasta_io():
    pa = PairwiseAligner("scoring_matrices/BLOSUM50.mat")
    file1 = "test_data/prot-0088.fa"
    file2 = "test_data/prot-0004.fa"
    pa.set_seqs(file1, file2)

    # check uppercase
    assert pa.get_seq1().isupper()
    assert pa.get_seq2().isupper()
    assert pa.get_seq1() == "YGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ"
    assert pa.get_seq2() == "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"

def test_scoring_matrix_io():
    pa = PairwiseAligner("scoring_matrices/BLOSUM50.mat")
    m = pa.get_score_mat()
    assert m.shape == (24, 24)
    assert not m.empty
    # Check that the row and column indices are correct
    assert (m.index == list("ARNDCQEGHILKMFPSTWYVBZX*")).all()
    assert (m.columns == list("ARNDCQEGHILKMFPSTWYVBZX*")).all()

def test_identical():
    # Using gap opening -5 and gap extension -1 (default)
    # Checked score using https://www.ebi.ac.uk/Tools/psa/emboss_needle/
    nw = NeedlemanWunsch("scoring_matrices/BLOSUM50.mat")
    seq = "test_data/prot-0088.fa"
    score, alignment = nw.align(seq, seq)
    assert alignment == "YGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ\nYGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ"
    
    sw = SmithWaterman("scoring_matrices/BLOSUM50.mat")
    seq = "test_data/prot-0088.fa"
    score, alignment = sw.align(seq, seq)
    assert alignment == "YGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ\nYGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ"
    

def test_alignment_score():
    nw = NeedlemanWunsch("scoring_matrices/BLOSUM50.mat")
    score, alignment = nw.align("test_data/test1.fa", "test_data/test2.fa")
    assert score == 4
    assert alignment == "BTN\nBT-"

    
    sw = SmithWaterman("scoring_matrices/BLOSUM50.mat")
    score, alignment = sw.align("test_data/test1.fa", "test_data/test2.fa")
    assert score == 10
    assert alignment == "BT\nBT"


