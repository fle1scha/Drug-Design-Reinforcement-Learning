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
        self.goal = input("Enter The Optimisation goal: \n eg 'FC(F)(Cl)C(F)(Cl)Cl'")
        self.mol = input("Enter the starting molecule or atom: \n" )
        self.molecule = Mol(self.mol, self.goal)
        self.molecule.GetRandomGoal()
        
        high = np.array([len(self.molecule.get_Atoms()),
                        len(self.molecule.get_Bonds())],
                        dtype=np.float32)
        
        # The action_space is defined. 
        self.action_space = spaces.Discrete(6)
        # Atom space determined from the goals atoms
        self.atom_space = self.molecule.get_Atoms()
        # Bond space based on goal, float values 1.0, 2.0, 1.5
        self.bond_space = self.molecule.get_Bonds()
        # The observation_space is defined.
        self.observation_space = spaces.Box(0,high,dtype=np.float32)

        self.seed()
        self.state = 0
        
   
    def step(self, action):
        
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg
        
        #This switch determines how the environment changes given the agent's action. 
        if action == 0:                                     # 0: Add random atom 
            self.molecule.AddA(rd.choice(self.atom_space))
            self.molecule.CheckValidity()                   # if false: mol will revert                        
        elif action == 1:                                   # 1: Revert back to previous state
            self.molecule.revertMol()                      
        elif action == 2:
            self.molecule.AddA("C")
            self.molecule.CheckValidity()
        elif action == 3:
                print("3")                                  # Could add bonds, specific mols based on weighting of MDP
        elif action == 4:
                print("4")
        else:
                print("5")

        #The observation state of the env is refreshed.
        self.state = self.molecule.GetSimilarity()
        
        #In order to fix each Episode only iterating once. We must ensure this doesn't evaluate to true after one iteration.
        done = bool(self.molecule.GetSimilarity() == 100.00)

        if done:
            reward = 0
        else:
            reward = self.CalculateReward()
        
        return np.array(self.state), reward, done, {}


    def reset(self):
        self.state = self.molecule.GetSimilarity()
        return np.array(self.state)

    def render(self):
        self.molecule.GetMol()

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    def CalculateReward(self):
        print("calculating Reward")
        return 1

