import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
import random
import copy


# TODO: original atom included

class Mol:
    """
    A class used to create and modify a molecule object.

    ...
    Attributes
    ----------
    start : str
        The starting molecule of the environment.

    target : str
        The optimisation goal molecule of the environment.

    RAM : list
        A list of the 3 most recent molecule states.

    bondmap : dictionary
        A dictionary of the different bonds in the molecule.

    modifications : list
        A list of possible modifications that can be made to the molecule.

    storage : list
        A list of additions that have been made to the molecule.

    present : boolean
        A boolean modifier that determines if a molecule library is present.

    Methods
    ----------
    __init__(self, mol, goal)
        The constructor method for the Mol class.

    get_random_goal(self)
        Returns a random optimisation goal.

    get_atoms(self)
        Returns a list of the possible atoms to construct the mol with.

    get_bonds(self)
        Returns a list of the possible bonds to construct the mol with.

    add_atom(self, Atom, back)
        Adds an Atom to the molecule.

    remove_atom(self, index)
        Removes an Atom from the molecule.

    history(self)
        Returns the last two states of the molecule.

    get_mol(self)
        Returns the Mol object

    check_validity(self)
        Uses implicit sanitization to check chemical validity of a molecule.

    display_changes(self)
        Shows the current and previous state of the mol object.

    def save_modifications(self)
        Saves the modifications of the agent

    def set_start_molecule(self, start)
        The mutator method for the start molecule of the environment.

    def set_target_molecule(self, target)
        The mutator method for the target molecule of the environment.
    """

    def __init__(self):
        """
        Parameters
        ----------
        mol : str
            The starting molecule of the environment, as a String.

        goal : str
            The optimisation goal of the environment (and the agent).
        """

        self.start = ""
        self.target = ""
        self.bondmap = {1.0: "", 1.5: "", 2.0: "=", 3.0: "#"}
        self.modifications = [self.start]

        if os.path.isfile('./MoleculeLibrary.csv'):
            self.df = pd.read_csv('MoleculeLibrary.csv')
            self.present = True

        else:
            self.present = False

    def set_start_molecule(self, start):
        """The mutator method for the start molecule of the environment."""
        self.start = start

    def set_target_molecule(self, target):
        """The mutator method for the target molecule of the environment."""

        self.target = target

    def get_random_molecule(self, is_target=True):  # if target = True, then chooses random target molecule
        """
        Returns a random goal as the optimisation goal for the molecule.

        Returns
        ----------
        randomMol['Compound ID'] : ?


        False : Boolean
            The method was unable to return a random goal.
        """
        if is_target:
            random_mol = self.df.sample()
            return (str(random_mol['SMILES']).split("\n")[0].split()[1])
        else:
            return (random.choice(self.get_atoms()))

    def get_atoms(self):
        """Accessor method for the list of possible bonds to construct the molecule with.
        Returns
        ----------
        list(atoms) : list
            A list of the atoms in the optimisation goal.
        """

        temp_mol = Chem.MolFromSmiles(self.target)
        atoms = set()
        for a in temp_mol.GetAtoms():
            atoms.add(a.GetSymbol())
        return list(atoms)

    def get_bonds(self):
        """Accessor method for the list of possible bonds to construct the molecule with.

        Returns
        ----------
        list(bonds) : list
            A list of possible bonds to construct the molecule with.
        """
        temptarget = Chem.MolFromSmiles(self.target)
        bonds = set()
        for a in temptarget.GetBonds():
            bonds.add(a.GetBondTypeAsDouble())
        return list(bonds)

    # Adding an Atom to the molecule
    def add_atom(self, front, back):
        """
        Adds an atom to the front or back of the molecule.

        Parameters
        ----------
        front : str
            The string represetntaion of the atom being added to the front.

        back : str
            The string represetntaion of the atom being added to the back.
            
        Returns
        ----------
        false : boolean
            If the molecule is chemically invalid and a change has not been made.

        true : boolean
            If the molecule is chemically valid and a change has been made.
            
        """
        current_molecule = self.start
        new_molecule = front + current_molecule + back
        if self.is_valid(new_molecule) == True:
            self.modifications.append(back)
            self.modifications.insert(0, front)
            self.modifications = list(filter(None, self.modifications))

            self.start = new_molecule
            return True
        else:
            return False

    # Adding an Atom to the molecule
    def add_brackets(self, typeof, index):
        """
        Adds a bracket to a sepcified atom within the molecule.

        Parameters
        ----------
        typeof : int
            either 1 for round or 2 for square

        index : str
            The string represetntaion of the atom being added to the back.
            
        Returns
        ----------
        false : boolean
            If the molecule is chemically invalid and a change has not been made.

        true : boolean
            If the molecule is chemically valid and a change has been made.
            
        """
        
        current_modifications = copy.copy(self.modifications)
        if len(current_modifications[index])>2:
            return False
        if typeof == 1:
            current_modifications[index] = "(" + current_modifications[index] + ")"
        else:
            current_modifications[index] = "[" + current_modifications[index] + "]"
        new_molecule = "".join(current_modifications)
        if self.is_valid(new_molecule) == True:
            self.modifications = current_modifications
            self.start = new_molecule
            return True
        else:
            return False

    # functionality to remove an atom   
    def remove_atom(self, index):
        """
        Removes an atom from the molecule.

        Parameters
        ----------
        index : int
            The index of the atom being removed.
        """
        if len(self.modifications) < 2:
            return False
        current_modifications = copy.copy(self.modifications)
        del current_modifications[index]
        new_molecule = "".join(current_modifications)
        if self.is_valid(new_molecule) == True:
            self.modifications = current_modifications
            self.start = new_molecule
            return True
        else:
            return False
        pass

    # Returning the Mol   
    def get_mol(self):
        """Accessor method for the current molecule's Mol object representation.

        Returns
        ----------
        Chem.MolFromSmiles(self.start) : Chem.Mol
            A Chem.Mol object.
        """
        return Chem.MolFromSmiles(self.start)

    # Uses implicit sanitisation to check chemical validity

    def is_valid(self, mol):
        """
        Checks the chemical validity of the molecule state.

        Returns
        ----------
        false : boolean
            If the molecule is chemically invalid.

        true : boolean
            If the molecule is chemically valid.
        """
        try:
            molecule = Chem.MolFromSmiles(mol)
            smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
            mol2 = Chem.MolFromSmiles(smiles)
        except:
            return False
        else:
            return True

    # The Taninoto Similarity
    def get_similarity(self):
        """
        Checks the Taninoto Similarity between two molecules.

        Returns
        ----------
        value : float, int
            Either a positive float measuring similarity or -1 if invalid.
        """
        mol1 = Chem.MolFromSmiles(self.start)
        mol2 = Chem.MolFromSmiles(self.target)
        fingerprint1 = Chem.RDKFingerprint(mol1)
        fingerprint2 = Chem.RDKFingerprint(mol2)
        return round(DataStructs.TanimotoSimilarity(fingerprint1, fingerprint2) * 100, 4)

    def save_modifications(self):
        """Saves the modifications of the agent."""
        with open("data/molecule_modifications.txt", "w+") as f:
            f.write('\n'.join('%s' % element for element in self.modifications))


if __name__ == '__main__':
    Mol()
