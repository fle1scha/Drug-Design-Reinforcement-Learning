import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
import random


class Mol:    
    def __init__(self, mol, goal):
        # Set instance variables
        self.mol = mol
        self.goal = goal
        self.RAM = [mol, mol, mol]
        self.bondmap = {1.0:"",1.5:"",2.0:"=",3.0:"#"}
        
        # store modifications added
        self.modifications = []
        
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
            self.mol = self.get_Atoms()[0]                      # Use different initial molecule
            self.RAM = [self.mol,self.mol,self.mol]
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
    
    def get_Bonds(self):
        molbonds = Chem.MolFromSmiles(self.goal).GetBonds()
        bonds = set() 
        for a in molbonds:
            bonds.add(a.GetBondTypeAsDouble())
        return list(bonds)
    
     # Adding an Atom to the molecule
    def AddA(self, Atom, back):
        self.modifications.append(Atom)
        self.updateRAM(self.mol)
        if back:
            newmol = self.mol + Atom
        else: 
            newmol =  Atom + self.mol 
        self.mol = newmol
        
    # The memory for recovering past states of the molecule
    def updateRAM(self, old):
        self.RAM[2] = self.RAM[1]
        self.RAM[1] = self.RAM[0]
        self.RAM[0] = old
        
    # Restore molecule to previous state
    def revertMol(self):
        self.mol = self.RAM[0]
        self.RAM[1] = self.RAM[2]
        self.RAM[0] = self.RAM[1]
                     
    # Return the last two states
    def history(self):
        return self.RAM
            
    # Returning the Mol   
    def GetMol(self):
        return Chem.MolFromSmiles(self.mol)
      
    # Uses implicit sanitisation to check chemical validity
    def CheckValidity(self):
        try:
            molecule = Chem.MolFromSmiles(self.mol)
            smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
            mol2 = Chem.MolFromSmiles(smiles)
        except:
            self.revertMol()
            return False
        else:
            return True
    
     # The Taninoto Similarity   
    def GetSimilarity(self):
        if self.CheckValidity():
            mol1 = Chem.MolFromSmiles(self.mol)
            mol2 = Chem.MolFromSmiles(self.goal)
            fingerprint1 = Chem.RDKFingerprint(mol1)
            fingerprint2 = Chem.RDKFingerprint(mol2)
            return round(DataStructs.TanimotoSimilarity(fingerprint1,fingerprint2) *100, 4) 
        else:
            return -1
    
    # Current and Previous state of the Mol
    def DisplayChanges(self):
        print("Molecule at state S-1: ")
        print(self.RAM[1])
        print("Molecule at state S: ")
        print(self.RAM[0])
        
        
if __name__ == '__main__':
       Mol()
   

