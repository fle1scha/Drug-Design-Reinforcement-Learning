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
        Type: Box(3)
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
        """
        The __init_ method for an env object.
        @parameters self:gym.Env
        """

        self.molecule = Chem.MolFromSmiles("C=C-C-C=C-C-O")
        
        
        # Set observation_space values based on self.molecule
        self.atoms = self.molecule.GetNumAtoms()
        self.bonds = self.molecule.GetNumBonds()
        self.conformers = self.molecule.GetNumConformers()
        
        print("Optimisation goal:")
        print("C=C-C-C=C-C-O")
        print("A:", self.atoms)
        print("B:", self.bonds)
        print("C:", self.conformers)
        print()
        
        

        # The high (max) values for the observation_space values
        high = np.array([self.atoms,
                        self.bonds,
                        self.conformers],
                        dtype=np.float32)
        
        # The action_space is defined. 
        self.action_space = spaces.Discrete(6)
        
        # The observation_space is defined.
        self.observation_space = spaces.Box(0, high, dtype=np.float32)

        # The following methods are just hollow method calls. They will be implemented as the prototype evolves.
        #----------------------------------------------------
        
        # Get the possible atoms to add to the molecule
        self.atom_space = self.Get_possible_Atoms()
        
        # Get the possible actions that the Agent can take
        self.action_options = self.Get_possible_Actions()
        
        # Get the possible atoms to add to the molecule
        self.bond_space = self.Get_possible_Bonds()
        
        #-----------------------------------------------------
        
        self.seed()
        self.state = None
   
   
    def step(self, action):
        """
        The step method considers the current state of the environment, and the chosen action of the agent, to adjust the env's 
        state accordingly.
        @parameters self:gym.Env, action: int 
        @returns np.array(self.state), reward, done
        """
        
        #Assert the validity of the agent's action. 
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg

        atoms, bonds, conformers = self.state
        
        # This switch determines how the environment changes given the agent's action. 
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

        # The observation state of the env is refreshed.
        self.state = (atoms, bonds, conformers)
        
        # Evaluate the terminality of the environment state. This aligns with the optimisation goal.
        done = bool(
            atoms < 0 or
            atoms >= self.atoms or
            bonds < 0 or
            bonds >= self.bonds or
            conformers < 0 or
            conformers >= self.conformers
        )
        
        #For prototyping purposes.
        done = False
        
        validity = self.Check_Validity()
        valency = self.GetValency()
        
        if done:
            reward = 0
        else:
            reward = self.CalculateReward(validity, valency)
        
        return np.array(self.state), reward, done, {}


    def reset(self):
        """
        The reset method refreshes the state of the environment.
        @parameters self: gym.env
        @returns np.array(self.state)
        """
        
        self.state = self.np_random.randint(low=0, high=(max(self.atoms, self.bonds, self.conformers)/2), size=(3,))
        print("Env original state: ")
        print(self.state)
        print()
        
        return np.array(self.state)

    def render(self):
        """
        The render method allows the user to intepret/visualise the env's change in state.
        @parameters self: gym.env
        @return [not yet defined]
        """
        return Chem.MolFromSmiles('C=C-C-C=C-C-O')

    def seed(self, seed=None):
        """
        The seed method maintains pseudo-random number generation to use throughout the training environment.
        @parameters self: gym.env, seed: float
        @returns seed
        """
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    
    """The below methods are still in development and are not needed to demonstrate the functionality of our prototype."""
    
    # __init__ method sub methods.
    def Get_possible_Atoms(self):
        atoms = ['C', 'N', 'O', 'S', 'Cl'] 
        return atoms
    
    def Get_possible_Bonds(self):
        bonds = ["double", "Single", "Triple"]
        return bonds
    
    def Get_possible_Actions(self):
        actions = [0,1,2,3,4,5]
        return 
    
    # step method sub methods.
    elements = Chem.GetPeriodicTable()
    def GetValency(self):
        print("... Calculating molecular valency.")
        return 1
    
    def Check_Validity(self):
        print("... Checking molecular validity.")
        valid = True
        return valid
    
    
    def CalculateReward(self, validity, valency):
        print("... Calculating reward based on chemical validity and molecule valency.")
        return 1
        
    