import pytest
from align import algs

@pytest.fixture
def create_pa_instance():
    return PairwiseAligner("/Users/mtsui1/Documents/Classes/Algs/Project1/scoring_matrices/BLOSUM50.mat")

def test_fasta_io():
    pa = create_pa_instance()
    file1 = "/Users/mtsui1/Documents/Classes/Algs/Project1/sequences/prot-0088.fa"
    file2 = "/Users/mtsui1/Documents/Classes/Algs/Project1/sequences/prot-0004.fa"
    pa.set_seqs(file1, file2)
    assert pa.get_seq1().isupper()
    assert pa.get_seq2().isupper()
    assert pa.get_seq1() == "YGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ"
    assert pa.get_seq2() == "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"

def test_scoring_matrix_io():
    pa = create_pa_instance()
    m = pa.get_score_mat()
    assert m.shape == (24, 24)
    assert not m.empty
    assert (m.index == list("ARNDCQEGHILKMFPSTWYVBZX*")).all()
    assert (m.columns == list("ARNDCQEGHILKMFPSTWYVBZX*")).all()

def test_identical():
    # Using gap opening -5 and gap extension -1 (default)
    # Checked score using https://www.ebi.ac.uk/Tools/psa/emboss_needle/

    nw = NeedlemanWunsch("/Users/mtsui1/Documents/Classes/Algs/Project1/scoring_matrices/BLOSUM50.mat")
    seq = "/Users/mtsui1/Documents/Classes/Algs/Project1/sequences/prot-0088.fa"
    score, alignment = nw.align(seq, seq)
    assert score == 870
    assert alignment == "YGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ\nYGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ"
    
    sw = SmithWaterman("/Users/mtsui1/Documents/Classes/Algs/Project1/scoring_matrices/BLOSUM50.mat")
    seq = "/Users/mtsui1/Documents/Classes/Algs/Project1/sequences/prot-0088.fa"
    score, alignment = sw.align(seq, seq)
    assert score == 870
    assert alignment == "YGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ\nYGKNQREAAQMDMVNDGVEDLRGKYVTLIYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLLLIHQVLAPGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRPINGNGKQ"
    

def test_alignment_score():
    assert True

