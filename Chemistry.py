import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
import random


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

    df : ?

    present : boolean
        A boolean modifier that determines if a molecule library is present.

    Methods
    ----------
    __init__(self, mol, goal)
        The constructor method for the Mol class.

    GetRandomGoal(self)
        Returns a random optimisation goal.

    get_Atoms(self)
        Returns a list of the possible atoms to construct the mol with.

    get_Bonds(self)
        Returns a list of the possible bonds to construct the mol with.

    AddA(self, Atom, back)
        Adds an Atom to the molecule.

    RemoveA(self, index)
        Removes an Atom from the molecule.

    updateRAM(self, old)
     Updates the memory for recovering past states of the molecule.

    revertMol(self)
        Restores molecule to previous state.

    history(self)
        Returns the last two states of the molecule.

    GetMol(self)
        Returns the Mol object

    CheckValidity(self)
        Uses implicit sanitization to check chemical validity of a molecule.

    DisplayChanges(self)
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
        self.modifications = []    # building blocks of the molecule
        if self.mol == "1":
                self.mol = random.choice(self.get_Atoms())

     
        #  Connect library if present
        if os.path.isfile('./MoleculeLibrary.csv'):
            self.df = pd.read_csv('MoleculeLibrary.csv')
            self.present = True

        else: 
            self.present = False
        
    # Returns Random Goal from library data frame   
    
    def GetRandomMolecule(self):
        """Returns a random goal as the optimisation goal for the molecule.

        Returns
        ----------
        randomMol['Compound ID'] : ?


        False : Boolean
            The method was unable to return a random goal.
        """

        if self.present:
            randomMol = self.df.sample()

            self.goal = str(randomMol['SMILES']).split("\n")[0].split()[1]    # easy substring to remove object type   
            return randomMol['Compound ID']

        else:
            return False

    def get_Atoms(self):
        """Accessor method for the list of possible bonds to construct the molecule with.
        Returns
        ----------
        list(atoms) : list
            A list of the atoms in the optimisation goal.
        """

        tempmol = Chem.MolFromSmiles(self.goal)
        atoms = set()
        for a in tempmol.GetAtoms():
            atoms.add(a.GetSymbol())
        return list(atoms)

    def get_Bonds(self):
        """Accessor method for the list of possible bonds to construct the molecule with.

        Returns
        ----------
        list(bonds) : list
            A list of possible bonds to construct the molecule with.
        """
        molbonds = Chem.MolFromSmiles(self.goal).GetBonds()
        bonds = set()
        for a in molbonds:
            bonds.add(a.GetBondTypeAsDouble())
        return list(bonds)

    
     # Adding an Atom to the molecule
    def AddAtom(self, front, back):
        """
        Adds an atom to the front or back of the molecule.

        Parameters
        ----------
        Atom : str
            The string represetntaion of the atom being added.

        back : boolean
            Whether the atom should be added to the back of the molecule or not.
        """
        current_molecule = self.mol
        new_molecule = front + current_molecule + back
        if self.isValid(new_molecule) == True:
            self.modifications.append(front + back)
            self.mol = new_molecule
            return True
        else:
            return False
        
    # functionality to remove an atom   
    def RemoveA(self,index):
        """
        Removes an atom from the molecule.

        Parameters
        ----------
        index : int
            The index of the atom being removed.
        """
        pass
        

    # Returning the Mol   
    def GetMol(self):
        """Accessor method for the current molecule's Mol object representation.

        Returns
        ----------
        Chem.MolFromSmiles(self.mol) : Chem.Mol
            A Chem.Mol object.
        """

        return Chem.MolFromSmiles(self.mol)

    # Uses implicit sanitisation to check chemical validity

    def isValid(self, mol):
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

    def CheckGoal(self):
        """
        Unsure as to what this does.
        """
        try:
            molecule = Chem.MolFromSmiles(self.goal)
            smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
        except:
            self.GetRandomGoal()
        else:
            pass

    # The Taninoto Similarity
    def GetSimilarity(self):
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
        return round(DataStructs.TanimotoSimilarity(fingerprint1,fingerprint2) *100, 4) 
        


if __name__ == '__main__':
  Mol()
