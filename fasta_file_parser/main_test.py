import re
import unittest
from unittest.mock import mock_open, patch
from collections import Counter
from main import parse_fasta, getSequenceType, find_ORFs, getORFs, computeCounts, printSequenceDetails, printSummary

class TestFastaProcessing(unittest.TestCase):
    def test_parse_fasta(self):
        fasta_content = """>seq1\nCGAATATTAGA\n>seq2\nGUGAAUGUAG\n"""
        expected_output = [(">seq1", "CGAATATTAGA"), (">seq2", "GUGAAUGUAG")]
        with patch("builtins.open", mock_open(read_data=fasta_content)):
            result = list(parse_fasta("fake_path"))
        self.assertEqual(result, expected_output)

    def test_getSequenceType(self):
        self.assertEqual(getSequenceType("ATTAG"), "DNA")
        self.assertEqual(getSequenceType("AUGAUUCGAA"), "RNA")
        self.assertEqual(getSequenceType("TATATUUUUUUGCG"), "Invalid")

    def test_find_ORFs(self):
        sequence = "ATGAAACCCGGGTTTTAA"
        expected_orfs = [(0, 15, 'ATGAAACCCGGGTTTTAA')]
        result = find_ORFs(sequence, re.compile(r"ATG"), re.compile(r"TAG|TAA|TGA"))
        self.assertEqual(result, expected_orfs)

        sequence = "AUGAAACCCGGGUUUUAA"
        expected_orfs = [(0, 15, 'AUGAAACCCGGGUUUUAA')]
        result = find_ORFs(sequence, re.compile(r"AUG"), re.compile(r"UAA|UAG|UGA"))
        self.assertEqual(result, expected_orfs)

    def test_getORFs(self):
        dna_sequence = "ATGAAACCCGGGTTTTAA"
        rna_sequence = "AUGAAACCCGGGUUUUAA"
        invalid_sequence = "ZEBRA"

        dna_expected_orfs = [(0, 15, 'ATGAAACCCGGGTTTTAA')]
        rna_expected_orfs = [(0, 15, 'AUGAAACCCGGGUUUUAA')]
        result = getORFs(dna_sequence, "DNA")
        self.assertEqual(result, dna_expected_orfs)
        result = getORFs(rna_sequence, "RNA")
        self.assertEqual(result, rna_expected_orfs)
        self.assertEqual(getORFs(invalid_sequence, "Invalid"), [])

    def test_computeCounts(self):
        dna_sequence = "ATGAAACCCGGGTTTTAA"
        rna_sequence = "AUGAAACCCGGGUUUUAA"
        ambiguous_sequence = "ZEBRACROSSING"
        self.assertEqual(computeCounts(dna_sequence), "A=6, C=3, G=4, T=5")
        self.assertEqual(computeCounts(rna_sequence), "A=6, C=3, G=4, U=5")
        self.assertEqual(computeCounts(ambiguous_sequence), "A=1, C=1, G=1, Ambiguous=10")

    @patch("builtins.print")
    def test_printSequenceDetails(self, mock_print):
        sequence_metadata = {"DNA": {"lenSum": 0, "count": 0}, "RNA": {"lenSum": 0, "count": 0}, "Invalid": {"lenSum": 0, "count": 0}}
        printSequenceDetails(">seq1", "ATGAAACCCGGGTTTTAA", sequence_metadata)

        self.assertEqual(sequence_metadata["DNA"], {"lenSum": 18, "count": 1})
        mock_print.assert_any_call(">seq1")
        mock_print.assert_any_call("Type: DNA")

    @patch("builtins.print")
    def test_printSummary(self, mock_print):
        sequence_metadata = {
            "DNA": {"lenSum": 46, "count": 2},
            "RNA": {"lenSum": 23, "count": 1},
            "Invalid": {"lenSum": 0, "count": 0},
        }
        printSummary(sequence_metadata)

        mock_print.assert_any_call("--SUMMARY--")
        mock_print.assert_any_call("Valid DNA sequences: 2")
        mock_print.assert_any_call("Mean DNA length: 23.00")
        mock_print.assert_any_call("Valid RNA sequences: 1")
        mock_print.assert_any_call("Mean RNA length: 23.00")
        mock_print.assert_any_call("Valid Invalid sequences: 0")

if __name__ == "__main__":
    unittest.main()