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
        The similarity between two Mol objects, as a float.

    current_reward : int
        The current reward to the agent.

    valid_step : Boolean
        True if the agent's action is valid.

    molecule : Mol object
        A user created class to hold the building methods for the environment mol object.

    mol : Chem.Mol object
        The current build state of the environment's molecule.

    goal : Chem.Mol object
        The optimisation state of the environment's molecule.

    self.action_space : spaces.Discrete(6)
        The action space of the environment.

    self.atom_space : list
        The atom space - a list of atoms contained in the goal.

    self.bond_space : list
        The bond space - a list of bonds that exist in the goal.

    self.observation_space : spaces.Box()
        The environment's observation space.

    self.state : float
        The current similarity of the current mol and the goal mol.

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

    calculate_reward(self)
        Determines the reward to be returned to the agent.

    update_policy(self)
        Changes the policy of the agent.
    """
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
        """
        Parameters
        ----------
        mol : str
            The starting molecule of the environment.

        goal: str
            The optimisation goal of the environment (and the agent).

        similarity : float
            The percentage similarity of the mol to the goal, represented as a float.
        """
        if target == "1":
            self.molecule = Mol(start, "F")
            self.molecule.get_random_molecule()
        else:
            self.molecule = Mol(start, target)

        self.mol = self.molecule.mol
        self.goal = self.molecule.goal
        self.similarity = goal
        
        self.similarity = goal       
        self.valid_step = True 
         
        high = np.array([len(self.molecule.get_atoms()), len(self.molecule.get_bonds())], dtype=np.float32)
            
        self.action_space = spaces.Discrete(6)
        self.observation_space = spaces.Box(0,high,dtype=np.float32)
        self.atom_space = self.molecule.get_atoms()
        self.bond_space = self.molecule.get_bonds()

        self.done = False
        self.seed()
        self.state = 0
   
    def step(self, action):
        """The method to advance the environment by one step. A step takes in the agent's action and changes the environment's state accordingly.

        Parameters
        ----------
        action : int
            The action selected by the agent given the previous reward.

        Returns
        ----------
        np.array(self.state)
            A numpy array of the state of the environment the done boolean.

        reward
            The reward to the agent.
        done
            The done boolean - it hows whether the environment has reached its goal.

        self.molecule.modifications
             The array of valid modifications for the agent.
        """
        
        """
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

        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg

        if action == 0:                              
            self.valid_step = self.molecule.add_atom("", rd.choice(self.atom_space)) 
            
        elif action == 1:                              
            self.valid_step = self.molecule.add_atom(rd.choice(self.atom_space), "")
            
        elif action == 2:                              
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = bond + rd.choice(self.atom_space)
            self.valid_step = self.molecule.add_atom("", bondedatom )
            
        elif action == 3:                            
            bond = self.molecule.bondmap[rd.choice(self.bond_space)]
            bondedatom = rd.choice(self.atom_space) + bond
            self.valid_step = self.molecule.add_atom(bondedatom,"")
            
        elif action == 4:
            self.validstep = self.molecule.add_brackets(rd.choice(range(1,3)), rd.choice(range(len(self.molecule.modifications))))
        
        elif action == 5:
            # TODO: remove atom from molecule
            pass   
        elif action == 6:
            pass    
        else:
            # TODO: remove bond from molecule
            pass
            
        self.state = self.molecule.get_similarity()
        self.done = bool(self.molecule.get_similarity() >= (self.similarity*100))
        reward = self.calculate_reward()
        
        return np.array(self.state), reward, self.done, self.molecule.modifications


    def reset(self):
        """Resets the state of the environment.
        """
     
        self.molecule.goal = self.goal
        self.molecule.mol = self.mol
        self.molecule.modifications = [self.mol]
        return self.molecule.get_similarity()

    def render(self):
        """
        Graphically renders the current state of the environment.
        """
        SMILES = self.molecule.get_mol()
        Image = Draw.MolToImage(SMILES, size=(300, 300))
        npFormat = np.asarray(Image)
        plt.imshow(npFormat)
        plt.draw()
        plt.pause(0.25) # pause how many seconds
        plt.close()


    def seed(self, seed=None):
        """Sets a pseudo-random seed for the environment.

        Parameters
        ----------
        seed : int, optional
            The seed value for the environment. None is the default.

        Returns
        ----------
        [seed] : int
            The pseudo-random number that is generated given the seed value.
        """

        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    
    def calculate_reward(self):
        """Returns the current reward, unless the step is invalid.

        Returns
        ----------
        self.current_reward : int
            The reward that will be given to the agent.
        """
        
        if self.done:
            return 100
        elif self.valid_step:
            return 10
        else:
            return -10
        
    
    def update_policy(self):
        """Shows that the policy of the agent has been updated.

        Returns
        ----------
        policy_str : str
            Confirmation message of the change.
        """

        policy_str = "policy received from agent"
        return policy_str
    
       