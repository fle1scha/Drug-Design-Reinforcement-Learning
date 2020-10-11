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
    A class used to create a molecule environment.
    ...

    Attributes
    ----------
    similarity : float

    currentReward : int

    validStep : Boolean

    molecule : Mol object

    mol : Chem.Mol object

    goal : Chem.Mol object

    self.action_space : spaces.Discrete(6)

    self.atom_space : list

    self.bond_space : list

    self.observation_space : spaces.Box()

    self.state : int

    Methods
    ----------
    __init__(self, mol, goal, similarity)
        The constructor method for this class. Called when gym.make() is called.

    step(self, action)
        Used to update the environment's state given an action from the agent.

    reset(self)
        Resets the environment's state.

    render(self)
        Graphically shows the environment's current state.

    seed(self, seed=None)
        Used to set a pseudo-random seed for the environment.

    CalculateReward(self)
        Determines the reward to be returned to the agent.

    updatepolicy(self)
        Changes the policy of the agent.
    """
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
 
    def __init__(self, mol, goal, similarity):
        """
        Parameters
        ----------
        mol : Chemistry.Mol object
            The starting molecule of the environment.

        goal: Chemistry.Mol object
            The optimisation goal of the environment (and the agent).

        similarity : float
            The percentage similarity of the mol to the goal, represented as a float.
        """

        self.similarity = similarity       # Float comparison 
        self.currentReward = 0             # this iterations reward
        self.validstep = True              # For reward calculation 
        if goal == "1": 
            self.molecule = Mol("F", "F")
            name = self.molecule.GetRandomGoal()
            print(name)
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
   
    def step(self, action):
        """
        The method to advance the environment by one step. A step takes in the agent's action
        and changes the environment's state accordingly.

        Parameters
        ----------
        action : int
            The action selected by the agent given the previous reward.


        """
        self.currentReward = 0
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg
        # This switch determines how the environment changes given the agent's action
        if action == 0:                                 # 0: Add random atom to back
            self.molecule.AddA(rd.choice(self.atom_space), True)
            self.validstep = self.molecule.CheckValidity() 
            
        elif action == 1:                               # Add random atom to the front
            self.molecule.AddA(rd.choice(self.atom_space), False)
            self.validstep = self.molecule.CheckValidity()
            
        elif action == 2:                               # Add random atom with random bond to the back
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = bond + rd.choice(self.atom_space)
            self.molecule.AddA(bondedatom, True)
            self.validstep = self.molecule.CheckValidity()
            
        elif action == 3:                               # Add random atom with random bond to the front
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = rd.choice(self.atom_space) + bond
            self.molecule.AddA(bondedatom, False)    
            self.validstep = self.molecule.CheckValidity()
            
            
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

    def render(self):
        if self.molecule.CheckValidity() == True:
            SMILES = self.molecule.GetMol()
            Image = Draw.MolToImage(SMILES, size=(300, 300))
            npFormat = np.asarray(Image)
            plt.imshow(npFormat)
            plt.draw()
            plt.pause(0.5) # pause how many seconds
            plt.close()
        else:
            print("Cannot render this molecule")

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    def CalculateReward(self):
        if self.validstep == False:
            self.currentReward -= 10   # invalid penalty
        return self.currentReward
    
    def updatepolicy(self):
        return "policy received from agent"
    
       