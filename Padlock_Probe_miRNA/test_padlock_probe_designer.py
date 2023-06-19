from Bio.Seq import Seq
import unittest
from Padlock_Probe_Designer import short_formula, long_formula, melting_temp, annealing_temp, index_with_lowest_at, get_padlock_arms


class TestMeltingTemp(unittest.TestCase):
    def test_annealing_temp(self):
        aseq = "ATCGGCTA"
        bseq = "CGGCTA"
        expected_temp = 4.58  # Update the expected temperature based on your calculations
        result = annealing_temp(aseq, bseq)
        self.assertAlmostEqual(result, expected_temp, places=2)

    def test_get_padlock_arms(self):
        miRNA = Seq("ATCGGCTA")
        expected_arm_a = Seq("ATCGGCTA")  # Update the expected arm_a value to be a Seq object
        expected_arm_b = Seq("")  # Update the expected arm_b value to be a Seq object
        result = get_padlock_arms(miRNA)
        self.assertEqual(result[0], expected_arm_a)

    def test_index_with_lowest_at(self):
        cdna = Seq("ATCGGCTATCGGCTA")  # Convert the input to a Seq object
        expected_index = 4
        result = index_with_lowest_at(cdna)
        self.assertEqual(result, expected_index)

    def test_long_formula(self):
        dna_seq = Seq("ATCGGCTATCGGCTA")  # Convert the input to a Seq object
        expected_tm = 41.94  # Update the expected melting temperature based on your calculations
        result = long_formula(dna_seq)
        self.assertAlmostEqual(result, expected_tm, places=2)

    def test_melting_temp_long_sequence(self):
        dna_seq = Seq("ATCGGCTATCGGCTA")  # Convert the input to a Seq object
        expected_tm = 41.94  # Update the expected melting temperature based on your calculations
        result = melting_temp(dna_seq)
        self.assertAlmostEqual(result, expected_tm, places=2)


if __name__ == '__main__':
    unittest.main()
