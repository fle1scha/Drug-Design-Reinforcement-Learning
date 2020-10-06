import gym
from gym import error, spaces, utils
from gym.utils import seeding
from gym import spaces

from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import random as rd
from Chemistry import Mol


class MoleculeEnvironment(gym.Env):
    """
    Observation:
        Type: Box(4)
        Num     Observation               Min          Max
        0       Number of Atoms           0            Number of Atoms in Molecule
        1       Number of Bonds           0            Number of Bonds in Molecule
        2       Number of Conformers      0            Number of Conformers in Molecule
        
    Actions:
        Type: Discrete(6)
        Num   Action
        0     Add Atom
        1     Remove Atom
        2     Add Bond
        3     Remove Bond
        4     Add Conformer
        5     Remove Conformer
    """
    def __init__(self):
        pass
   
    def step(self, action):   
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg
        startstate = self.molecule.GetSimilarity()
        # This switch determines how the environment changes given the agent's action
        if action == 0:                                 # 0: Add random atom to back
            self.molecule.AddA(rd.choice(self.atom_space), True)
            self.molecule.CheckValidity() 
            
        elif action == 1:                               # Add random atom to the front
            self.molecule.AddA(rd.choice(self.atom_space), False)
            self.molecule.CheckValidity()
            
        elif action == 2:                               # Add random atom with random bond to the back
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = bond + rd.choice(self.atom_space)
            self.molecule.AddA(bondedatom, True)
            self.molecule.CheckValidity()
            
        elif action == 3:                               # Add random atom with random bond to the front
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = rd.choice(self.atom_space) + bond
            self.molecule.AddA(bondedatom, False)    
            self.molecule.CheckValidity()
            
        elif action == 4:
            print("add bracketed atom")
            self.molecule.modifications.append("bracketed atom")
        else:
            print("add ring")
            self.molecule.modifications.append("ring")
            
        print(self.molecule.modifications)

        # The observation state of the env is refreshed.
        self.state = self.molecule.GetSimilarity()
        if self.state < startstate:
            self.molecule.revertMol()
        
        #In order to fix each Episode only iterating once. We must ensure this doesn't evaluate to true after one iteration.
        done = bool(self.molecule.GetSimilarity() == 100.00)

        if done:
            reward = 0
        else:
            reward = self.CalculateReward()
        
        return np.array(self.state), reward, done, {}


    def reset(self):
        self.molecule.goal = self.goal
        self.molecule.mol = self.mol
        self.molecule.modifications = []

    def render(self):
        self.molecule.GetMol()

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    def CalculateReward(self):
        print("calculating Reward")
        return 1
    
    def SetState(self, mol, goal):
        if goal == "1": 
            self.molecule = Mol("F", "F")
            self.molecule.GetRandomGoal()
        else:
            self.molecule = Mol(mol,goal)
            self.molecule.CheckGoal()
            
        self.mol = self.molecule.mol
        self.goal = self.molecule.goal
        
        high = np.array([len(self.molecule.get_Atoms()),
                        len(self.molecule.get_Bonds())],
                        dtype=np.float32)
                
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

