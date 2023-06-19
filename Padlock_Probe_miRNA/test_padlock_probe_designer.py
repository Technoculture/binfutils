import unittest
from Bio.Seq import Seq
from io import StringIO
from contextlib import redirect_stdout
from Padlock_Probe_Designer import short_formula, long_formula, melting_temp, annealing_temp, index_with_lowest_at, get_padlock_arms

class TestMeltingTemp(unittest.TestCase):

    def test_short_formula(self):
        self.assertEqual(short_formula("ACGT"), 16)

    def test_long_formula(self):
        self.assertEqual(long_formula("ACGT" * 14), 80.94)

    def test_melting_temp_short_seq(self):
        self.assertEqual(melting_temp("ACGT"), 16)

    def test_melting_temp_long_seq(self):
        self.assertEqual(melting_temp("ACGT" * 14), 80.94)

    def test_annealing_temp(self):
        self.assertAlmostEqual(annealing_temp("ACGT", "TGCA"), 7.24, places=2)

    def test_index_with_lowest_at(self):
        self.assertEqual(index_with_lowest_at("ACGTTGCA"), 4)

    def test_get_padlock_arms(self):
        miRNA = Seq("ACGTTGCA")
        arm1, arm2 = get_padlock_arms(miRNA)
        self.assertEqual(arm1, Seq("ACGTT"))
        self.assertEqual(arm2, Seq("GCA"))

    def test_design_padlock_probe(self):
        # Test the design_padlock_probe function by capturing its output
        output = StringIO()
        with redirect_stdout(output):
            design_padlock_probe()
        result = output.getvalue().strip()
        self.assertEqual(result, "Expected output string")

if __name__ == '__main__':
    unittest.main()
