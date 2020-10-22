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
    mol : str
        The starting molecule of the environment.

    goal : str
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
    """

    def __init__(self, mol, goal):
        """
        Parameters
        ----------
        mol : str
            The starting molecule of the environment, as a String.

        goal : str
            The optimisation goal of the environment (and the agent).
        """
        # Set instance variables
        self.mol = mol     # The starting canvas molecule
        self.goal = goal   # molecule to represent the optimised state
        self.bondmap = {1.0:"",1.5:"",2.0:"=",3.0:"#"} 
        self.modifications = [self.mol]    # building blocks of the molecule
        if self.mol == "1":
                self.mol = random.choice(self.get_atoms())

     
        #  Connect library if present
        if os.path.isfile('./MoleculeLibrary.csv'):
            self.df = pd.read_csv('MoleculeLibrary.csv')
            self.present = True

        else: 
            self.present = False
        
    # Returns Random Goal from library data frame   
    
    def get_random_molecule(self):
        """
        Returns a random goal as the optimisation goal for the molecule.
        Returns a random goal as the optimisation goal for the molecule.

        Returns
        ----------
        randomMol['Compound ID'] : ?


        False : Boolean
            The method was unable to return a random goal.
        """

        if self.present:
            random_mol = self.df.sample()

            self.goal = str(random_mol['SMILES']).split("\n")[0].split()[1]    # easy substring to remove object type   
            return random_mol['Compound ID']

        else:
            return False

    def get_atoms(self):
        """Accessor method for the list of possible bonds to construct the molecule with.
        Returns
        ----------
        list(atoms) : list
            A list of the atoms in the optimisation goal.
        """

        temp_mol = Chem.MolFromSmiles(self.goal)
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
        mol_bonds = Chem.MolFromSmiles(self.goal).GetBonds()
        bonds = set()
        for a in mol_bonds:
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
        current_molecule = self.mol
        new_molecule = front + current_molecule + back
        if self.is_valid(new_molecule) == True:
            self.modifications.append(back)
            self.modifications.insert(0,front)
            self.modifications = list(filter(None, self.modifications))

            self.mol = new_molecule
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
        if typeof == 1:
            current_modifications[index] = "(" + current_modifications[index] + ")"
        else:
            current_modifications[index] = "[" + current_modifications[index] + "]"
        new_molecule = "".join(current_modifications)
        if self.is_valid(new_molecule) == True:
            self.modifications = current_modifications
            self.mol = new_molecule
            return True
        else:
            return False
        
        
    # functionality to remove an atom   
    def remove_atom(self,index):
        """
        Removes an atom from the molecule.

        Parameters
        ----------
        index : int
            The index of the atom being removed.
        """
        current_modifications = copy.copy(self.modifications)
        del current_modifications[index]
        new_molecule = "".join(current_modifications)
        if self.is_valid(new_molecule) == True:
            self.modifications = current_modifications
            self.mol = new_molecule
            return True
        else:
            return False
        pass
        

    # Returning the Mol   
    def get_mol(self):
        """Accessor method for the current molecule's Mol object representation.

        Returns
        ----------
        Chem.MolFromSmiles(self.mol) : Chem.Mol
            A Chem.Mol object.
        """
        return Chem.MolFromSmiles(self.mol)

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
        mol1 = Chem.MolFromSmiles(self.mol)
        mol2 = Chem.MolFromSmiles(self.goal)
        fingerprint1 = Chem.RDKFingerprint(mol1)
        fingerprint2 = Chem.RDKFingerprint(mol2)
        return round(DataStructs.TanimotoSimilarity(fingerprint1,fingerprint2) * 100, 4) 
    
    def save_modifications(self):
        with open("data/molecule_modifications.txt", "w+") as f:
            f.write('\n'.join('%s' % element for element in self.modifications))
        


if __name__ == '__main__':
      Mol()
