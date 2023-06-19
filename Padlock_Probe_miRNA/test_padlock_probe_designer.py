import unittest
from Bio.Seq import Seq
from your_module_name import short_formula, long_formula, melting_temp, annealing_temp, index_with_lowest_at, get_padlock_arms

class TestMeltingTemp(unittest.TestCase):
    def test_short_formula(self):
        dna_seq = "ATCGGCTA"
        expected_tm = 24.0
        result = short_formula(dna_seq)
        self.assertEqual(result, expected_tm)

    def test_long_formula(self):
        dna_seq = "ATCGGCTATCGGCTA"
        expected_tm = 66.47
        result = long_formula(dna_seq)
        self.assertAlmostEqual(result, expected_tm, places=2)

    def test_melting_temp_short_sequence(self):
        dna_seq = "ATCGGCTA"
        expected_tm = 24.0
        result = melting_temp(dna_seq)
        self.assertEqual(result, expected_tm)

    def test_melting_temp_long_sequence(self):
        dna_seq = "ATCGGCTATCGGCTA"
        expected_tm = 66.47
        result = melting_temp(dna_seq)
        self.assertAlmostEqual(result, expected_tm, places=2)

    def test_annealing_temp(self):
        aseq = "ATCGGCTA"
        bseq = "CGGCTA"
        expected_temp = 11.17
        result = annealing_temp(aseq, bseq)
        self.assertAlmostEqual(result, expected_temp, places=2)

    def test_index_with_lowest_at(self):
        cdna = "ATCGGCTATCGGCTA"
        expected_index = 4
        result = index_with_lowest_at(cdna)
        self.assertEqual(result, expected_index)

    def test_get_padlock_arms(self):
        miRNA = Seq("ATCGGCTA")
        expected_arm_a = "ATCGGCTA"
        expected_arm_b = ""
        result = get_padlock_arms(miRNA)
        self.assertEqual(result[0], expected_arm_a)
        self.assertEqual(result[1], expected_arm_b)

if __name__ == '__main__':
    unittest.main()
