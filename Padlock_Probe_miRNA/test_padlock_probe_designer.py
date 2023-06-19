import unittest
from Padlock_Probe_Designer import melting_temp, annealing_temp, index_with_lowest_at, get_padlock_arms

class TestPadlockProbeDesigner(unittest.TestCase):
    def test_melting_temp(self):
        self.assertEqual(melting_temp('ATCG'), 58.0)
        self.assertEqual(melting_temp('GCGCTATA'), 44.0)
        self.assertEqual(melting_temp('TACGATCGATCG'), 48.0)
        self.assertEqual(melting_temp('CCCC'), 56.0)

    def test_annealing_temp(self):
        self.assertAlmostEqual(annealing_temp('ATCG', 'CGAT'), 39.05, places=2)
        self.assertAlmostEqual(annealing_temp('GCGCTATA', 'TATAGCGC'), 34.45, places=2)
        self.assertAlmostEqual(annealing_temp('TACGATCGATCG', 'CGATCGATCGTA'), 40.45, places=2)
        self.assertAlmostEqual(annealing_temp('CCCC', 'GGGG'), 56.0, places=2)

    def test_index_with_lowest_at(self):
        self.assertEqual(index_with_lowest_at('ATCGGCGAT'), 5)
        self.assertEqual(index_with_lowest_at('GCGCTATAGCGC'), 6)
        self.assertEqual(index_with_lowest_at('TACGATCGATCGTAC'), 11)
        self.assertEqual(index_with_lowest_at('CCCCGGGG'), 4)

    def test_get_padlock_arms(self):
        arm_a, arm_b = get_padlock_arms('ATCGGCGAT')
        self.assertEqual(arm_a, 'ATCGG')
        self.assertEqual(arm_b, 'CGAT')

        arm_a, arm_b = get_padlock_arms('GCGCTATAGCGC')
        self.assertEqual(arm_a, 'GCGCTATA')
        self.assertEqual(arm_b, 'GCGC')

        arm_a, arm_b = get_padlock_arms('TACGATCGATCGTAC')
        self.assertEqual(arm_a, 'TACGATCGATCG')
        self.assertEqual(arm_b, 'CGTAC')

        arm_a, arm_b = get_padlock_arms('CCCCGGGG')
        self.assertEqual(arm_a, 'CCCC')
        self.assertEqual(arm_b, 'GGGG')

if __name__ == '__main__':
    unittest.main()
