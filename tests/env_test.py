import numpy as np
from .context import molecule_env
from Chemistry import Mol
from molecule_agent import MoleculeAgent

import gym
from gym import spaces

"""
Test class for integration testing for molecule_env.py and molecule_agent.py class.

Methods
__________

test_init_mol() 
    Tests the __init__ mol function of the molecule_env class. 
    
test_init_goal()
    Tests the __init__ goal function of the molecule_env class.
    
test_action_space()
    Tests the initiation of the action space of the environment. 
    
test_main()
    Tests that an iteration of the environment-agent interaction is correct.
"""


def test_init_mol():
    """Tests whether the env correctly initiates.
    """
    test_molecule = Mol()
    test_molecule.set_start_molecule("C")
    test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
    env = gym.make("gym_molecule:molecule-v0", mol=test_molecule, goal=100)
    assert env.start_molecule == "C"


def test_init_goal():
    """Tests whether the env correctly sets the goal.
    """
    test_molecule = Mol()
    test_molecule.set_start_molecule("C")
    test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
    env = gym.make("gym_molecule:molecule-v0", mol=test_molecule, goal=100)
    assert env.target_molecule == "ClCC(Cl)(Cl)Cl"


def test_action_space():
    """Tests the format of the action space."""
    test_molecule = Mol()
    test_molecule.set_start_molecule("C")
    test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
    env = gym.make("gym_molecule:molecule-v0", mol=test_molecule, goal=100)
    assert env.action_space == spaces.Discrete(6)


def test_main():
    """Tests that an iteration of the environment-agent interaction is correct."""
    test_molecule = Mol()
    test_molecule.set_start_molecule("C")
    test_molecule.set_target_molecule("ClCC(Cl)(Cl)Cl")
    env = gym.make("gym_molecule:molecule-v0", mol=test_molecule, goal=100)
    agent = MoleculeAgent(env.observation_space, env.action_space)
    state = env.reset()
    reward = 0
    done = False
    action = agent.act(state, reward, done)

    assert env.step(action) == (np.array(env.state), 10, False, env.molecule.modifications)
