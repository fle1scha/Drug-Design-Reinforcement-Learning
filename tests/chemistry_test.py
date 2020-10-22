from Chemistry import Mol
import unittest


class test_chem(unittest.TestCase):
    """Test class for the Chemistry.py class.
    """

    def test_atoms(self):
        """Asserts whether the atom space is correct or not."""
        test_molecule = Mol()
        test_molecule.set_start_molecule("C")
        test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
        atom_set = set(test_molecule.get_atoms())
        self.assertEqual(atom_set, {"Cl", "C"})

    def test_bonds(self):
        """Asserts whether the bond space is correct or not."""
        test_molecule = Mol()
        test_molecule.set_start_molecule("C")
        test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
        self.assertEqual(test_molecule.get_bonds(), [1.0])

    def test_add_front(self):
        """Tests the functionality of adding an atom to the front of a molecule."""
        test_molecule = Mol()
        test_molecule.set_start_molecule("C")
        test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
        test_molecule.add_atom("C", "")
        self.assertEqual(test_molecule.start, "CC")

    def test_add_back(self):
        """Tests the functionality of adding an atom to the back of a molecule."""
        test_molecule = Mol()
        test_molecule.set_start_molecule("C")
        test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
        test_molecule.add_atom("", "C")
        self.assertEqual(test_molecule.start, "CC")

    def test_validity(self):
        test_molecule = Mol()
        test_molecule.set_start_molecule("CX")
        test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
        self.assertEqual(test_molecule.is_valid(test_molecule.start), False)


if __name__ == '__main__':
    unittest.main()
