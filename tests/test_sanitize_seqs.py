from Bio.Seq import Seq
from numpy.testing import assert_equal
from src.sanitize import remove_unrelated_sequences


def test_should_pass_dna_sequences(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    test_string = ">testing\nATGCTGA\n"
    with open(fasta_file, "w") as f:
        f.write(test_string)
    remove_unrelated_sequences(str(fasta_file), "N")
    with open(result_file) as f:
        assert f.read() == test_string

def test_should_pass_ptn_sequences(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    test_string = ">testing\nAYWRSQL\n"
    with open(fasta_file, "w") as f:
        f.write(test_string)
    remove_unrelated_sequences(str(fasta_file), "P")
    with open(result_file) as f:
        assert f.read() == test_string

def test_should_remove_dna_sequences_with_unknown(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    with open(fasta_file, "w") as f:
        f.write(">testing\nATGCTNA\n")
    remove_unrelated_sequences(str(fasta_file), "N")
    with open(result_file) as f:
        assert f.read() == ""

def test_should_remove_ptn_sequences_with_unknown(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    with open(fasta_file, "w") as f:
        f.write(">testing\nAXWRSQL\n")
    remove_unrelated_sequences(str(fasta_file), "P")
    with open(result_file) as f:
        assert f.read() == ""

def test_should_remove_dna_sequences_with_doubt(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    with open(fasta_file, "w") as f:
        f.write(">testing\nATGCTBA\n")
    remove_unrelated_sequences(str(fasta_file), "N")
    with open(result_file) as f:
        assert f.read() == ""

def test_should_remove_ptn_sequences_with_doubt(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    with open(fasta_file, "w") as f:
        f.write(">testing\nABWRSQL\n")
    remove_unrelated_sequences(str(fasta_file), "P")
    with open(result_file) as f:
        assert f.read() == ""

def test_should_only_remove_dna_sequences_not_supported(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    test_string = ">testing\nATGCTGA\n"
    with open(fasta_file, "w") as f:
        f.write(test_string + ">testing2\nATGCTBA\n")
    remove_unrelated_sequences(str(fasta_file), "N")
    with open(result_file) as f:
        assert f.read() == test_string

def test_should_only_remove_ptn_sequences_with_not_supported(tmp_path):
    fasta_file = tmp_path / "fasta.fasta"
    result_file = fasta_file.with_suffix(fasta_file.suffix + ".sanitized")
    test_string = ">testing\nAYWRSQL\n"
    with open(fasta_file, "w") as f:
        f.write(test_string + ">testing2\nABWRSQL\n")
    remove_unrelated_sequences(str(fasta_file), "P")
    with open(result_file) as f:
        assert f.read() == test_string