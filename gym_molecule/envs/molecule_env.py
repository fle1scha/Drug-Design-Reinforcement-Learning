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
        self.atoms = self.molecule.GetNumAtoms()
        self.bonds = self.molecule.GetNumBonds()
        self.conformers = self.molecule.GetNumConformers()

        # MAX values in observation space
        high = np.array([self.atoms,
                        self.bonds,
                        self.conformers],
                        dtype=np.float32)

        self.action_space = spaces.Discrete(6)
        self.observation_space = spaces.Box(0, high, dtype=np.float32)

        self.seed()
        self.state = None

    def step(self, action):
        err_msg = "%r (%s) invalid" % (action, type(action))
        assert self.action_space.contains(action), err_msg

        atoms, bonds, conformers = self.state
        if action == 0:
            atoms += 1
        elif action == 1:
            atoms -= 1
        elif action == 2:
            bonds += 1
        elif action == 3:
            bonds -= 1
        elif action == 4:
            conformers += 1
        else:
            conformers -= 1

        self.state = (atoms, bonds, conformers)
        print(self.state)

        done = bool(
            atoms < 0 or
            atoms >= self.atoms or
            bonds < 0 or
            bonds >= self.bonds or
            conformers < 0 or
            conformers >= self.conformers
        )

        if not done:
            reward = calculateReward()
        else:
            reward = 0

        return np.array(self.state), reward, done, {}


    def calculateReward():
        # OUT OF SCOPE
        pass

    def reset(self):
        self.state = self.np_random.uniform(low=0, high=min(self.atoms, self.bonds, self.conformers)/2, size=(3,))
        return np.array(self.state)

    def render(self):
        raise NotImplementedError

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
