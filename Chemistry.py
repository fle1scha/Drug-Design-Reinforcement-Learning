import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
import random
# TODO: original atom included

class Mol:    
    def __init__(self, mol, goal):
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
    def GetRandomGoal(self):
        if self.present:
            randomMol = self.df.sample()
            self.goal = str(randomMol['SMILES']).split("\n")[0].split()[1]    # easy substring to remove object type   
            return randomMol['Compound ID']
        else:
            return False
     
    # Returns a list of the possible atoms to construct the mol with 
    def get_Atoms(self):
        tempmol = Chem.MolFromSmiles(self.goal)
        atoms = set() 
        for a in tempmol.GetAtoms():
            atoms.add(a.GetSymbol())
        return list(atoms)
    
    # Returns a list of the possible bonds to construct the mol with 
    def get_Bonds(self):
        molbonds = Chem.MolFromSmiles(self.goal).GetBonds()
        bonds = set() 
        for a in molbonds:
            bonds.add(a.GetBondTypeAsDouble())
        return list(bonds)
    
     # Adding an Atom to the molecule
    def AddAtom(self, front, back):
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
        pass
        
       
    # Returning the Mol   
    def GetMol(self):
        return Chem.MolFromSmiles(self.mol)
      
    # Uses implicit sanitisation to check chemical validity
    def isValid(self, mol):
        try:
            molecule = Chem.MolFromSmiles(mol)
            smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
            mol2 = Chem.MolFromSmiles(smiles)
        except:
            return False
        else:
            return True
        
    def CheckGoal(self):
        try:
            molecule = Chem.MolFromSmiles(self.goal)
            smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
        except:
            self.GetRandomGoal()
        else:
            pass
    
     # The Taninoto Similarity   
    def GetSimilarity(self):
        mol1 = Chem.MolFromSmiles(self.mol)
        mol2 = Chem.MolFromSmiles(self.goal)
        fingerprint1 = Chem.RDKFingerprint(mol1)
        fingerprint2 = Chem.RDKFingerprint(mol2)
        return round(DataStructs.TanimotoSimilarity(fingerprint1,fingerprint2) *100, 4) 
        
        
if __name__ == '__main__':
       Mol()
   

