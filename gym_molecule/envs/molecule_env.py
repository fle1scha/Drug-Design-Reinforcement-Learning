import gym
from gym import error, spaces, utils
from gym.utils import seeding
from gym import spaces

from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import random as rd
from Chemistry import Mol
from PIL import Image
import numpy as np
import matplotlib.pyplot as plt


class MoleculeEnvironment(gym.Env):
    
    """
    Observation:
        Type: Box(4)
        Num     Observation                                  Min          Max
        0       Number of Atoms (starting molecule)          0            Number of Atoms (starting molecule)
        1       Number of Bonds (starting molecule)          0            Number of Bonds (starting molecule)
        2       Number of Atoms (target molecule)            0            Number of Atoms (target molecule)  
        3       Number of Bonds (target molecule)            0            Number of Bonds (target molecule) 
        
    Actions:
        Type: Discrete(6)
        Num   Action
        0     Add Atom to Back of Molecule
        1     Add Atom to Front of Molecule
        2     Add Atom with Bond to Back of Molecule
        3     Add Atom with Bond to Front of Molecule
        4     Add Bracketed Atom to Molecule
        5     Add Ring to Molecule
        6     Remove Atom from Molecule
        7     Remove Bond from Molecule
    """
        
    def __init__(self, start, target, goal):         
        if target == "1":
            self.molecule = Mol(start, "F")
            self.molecule.GetRandomGoal()
        else:
            self.molecule = Mol(start, target)

        self.mol = self.molecule.mol
        self.goal = self.molecule.goal
        
        self.similarity = goal             # Float comparison 
        self.currentReward = 0             # this iterations reward
        self.validstep = True 
         
        high = np.array([len(self.molecule.get_Atoms()), len(self.molecule.get_Bonds())], dtype=np.float32)
        
            
            
            
        # The action_space is defined. 
        self.action_space = spaces.Discrete(6)
        # Atom space determined from the goals atoms
        self.atom_space = self.molecule.get_Atoms()
        # Bond space based on goal, float values 1.0, 2.0, 1.5, 3.0
        self.bond_space = self.molecule.get_Bonds()
        # The observation_space is defined.
        self.observation_space = spaces.Box(0,high,dtype=np.float32)

        self.seed()
        self.state = 0
   
    def step(self, action): 
        self.currentReward = 0
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg
        # This switch determines how the environment changes given the agent's action
        if action == 0:                                 # 0: Add random atom to back
            self.validstep = self.molecule.AddAtom("", rd.choice(self.atom_space)) 
            
        elif action == 1:                               # Add random atom to the front
            self.validstep = self.molecule.AddAtom(rd.choice(self.atom_space), "")
            
        elif action == 2:                               # Add random atom with random bond to the back
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = bond + rd.choice(self.atom_space)
            self.validstep = self.molecule.AddAtom("", bondedatom )
            
        elif action == 3:                               # Add random atom with random bond to the front
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = rd.choice(self.atom_space) + bond
            self.validstep = self.molecule.AddAtom(bondedatom,"")
            
        elif action == 4:
            # print("add bracketed atom")
            self.molecule.modifications.append("bracketed atom")
        else:
            # print("add ring")
            self.molecule.modifications.append("ring")
            
        self.state = self.molecule.GetSimilarity()
        
        #In order to fix each Episode only iterating once. We must ensure this doesn't evaluate to true after one iteration.
        done = bool(self.molecule.GetSimilarity() >= (self.similarity))

        if done == True:
            reward = 100
        else:
            reward = self.CalculateReward()
        
        return np.array(self.state), reward, done, self.molecule.modifications


    def reset(self):
        self.molecule.goal = self.goal
        self.molecule.mol = self.mol
        self.molecule.modifications = []
        return self.molecule.GetSimilarity()

    def render(self):
        SMILES = self.molecule.GetMol()
        Image = Draw.MolToImage(SMILES, size=(300, 300))
        npFormat = np.asarray(Image)
        plt.imshow(npFormat)
        plt.draw()
        plt.pause(0.25) # pause how many seconds
        plt.close()


    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    def CalculateReward(self):
        if self.validstep == False:
            self.currentReward -= 10   # invalid penalty
        return self.currentReward
    
    def updatepolicy(self):
        return "policy received from agent"
    
       