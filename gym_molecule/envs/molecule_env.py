import gym
from gym import error, spaces, utils
from gym.utils import seeding
from gym import spaces

from rdkit import Chem
from rdkit.Chem import Draw

import numpy as np

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

        self.molecule = Chem.MolFromSmiles("C=C-C-C=C-C-O")
        
        #Set observation_space values based on self.molecule
        self.atoms = self.molecule.GetNumAtoms()
        self.bonds = self.molecule.GetNumBonds()
        self.conformers = self.molecule.GetNumConformers()

        #The high (max) values for the observation_space values
        high = np.array([self.atoms,
                        self.bonds,
                        self.conformers],
                        dtype=np.float32)

        #The action_space is defined. 
        self.action_space = spaces.Discrete(6)
        
        #The observation_space is defined.
        self.observation_space = spaces.Box(0, high, dtype=np.float32)

        self.seed()
        self.state = None
   
   
    def step(self, action):
        
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg

        atoms, bonds, conformers = self.state
        
        #This switch determines how the environment changes given the agent's action. 
        if action == 0:
            atoms += 1
        elif action == 1:
            atoms += 1
        elif action == 2:
            bonds += 1
        elif action == 3:
            bonds += 1
        elif action == 4:
            conformers += 1
        else:
            conformers += 1

        #The observation state of the env is refreshed.
        self.state = (atoms, bonds, conformers)
        
        #In order to fix each Episode only iterating once. We must ensure this doesn't evaluate to true after one iteration.
        done = bool(
            atoms < 0 or
            atoms >= self.atoms or
            bonds < 0 or
            bonds >= self.bonds or
            conformers < 0 or
            conformers >= self.conformers
        )

        if done:
            reward = 0
        else:
            reward = self.CalculateReward()
        
        return np.array(self.state), reward, done, {}


    def reset(self):
        self.state = self.np_random.uniform(low=0, high=min(self.atoms, self.bonds, self.conformers)/2, size=(3,))
        return np.array(self.state)

    def render(self):
        return Chem.MolFromSmiles('C1OC1')

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    
    # PROTOTYPE CODE methods 
    
    # init
    def Get_possible_Atoms():
        atoms = ['C', 'N', 'O', 'S', 'Cl'] 
        return atoms
    
    def Get_possible_Bonds():
        bonds = ["double", "Single", "Triple"]
        return bonds
    
    def Get_possible_actions():
        actions = [0,1,2,3,4,5]
        return 
    
    # step
    elements = Chem.GetPeriodicTable()
    def GetValency(element):
        return list(elements.GetValenceList(element))[0]
    
    def Check_Validity(state):
        print("checking Validity")
        valid = True
        return valid
    
    def CalculateReward(self):
        print("calculating Reward")
        return 1
        
    
    # 
    
    ## LUKE COMMENT