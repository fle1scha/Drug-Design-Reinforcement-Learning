from Chemistry import Mol
import unittest


class ChemTest(unittest.TestCase):
    """Test class for the Chemistry.py class.
    """

    def test_atoms(self):
        """Asserts whether the atom space is correct or not."""
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        with self.subTest():
            self.assertEqual(test_molecule.get_atoms(), ["Cl", "C"])
        with self.subTest():
            self.assertEqual(test_molecule.get_atoms(), ["C", "Cl"])

    def test_bonds(self):
        """Asserts whether the bond space is correct or not."""
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()

    def test_add_front(self):
        """Tests the functionality of adding an atom to the front of a molecule."""
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()

    def test_add_back(self):
        """Tests the functionality of adding an atom to the back of a molecule."""
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()

    def test_remove(self):
        """Tests the functionality of removing an atom from the molecule."""
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()

    def test_history(self):
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()

    def test_get_mol(self):
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()

    def test_validity(self):
        test_molecule = Mol("C", "ClCC(Cl)(Cl)Cl")
        self.assertEqual()


if __name__ == '__main__':
    unittest.main()
