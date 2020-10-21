import numpy as np

from molecule_agent import MoleculeAgent

import unittest
import gym
from gym import spaces


class AgentTest(unittest.TestCase):
    """Test class for molecule_agent.py class."""

    def test_action(self):
        """Tests whether the agent correctly returns an int action."""
        action_space = spaces.Discrete(6)
        observation_space = spaces.Box(0, 10, dtype=np.float32)
        agent = MoleculeAgent(action_space, observation_space)
        self.assertEqual(agent.act(1, 1, 1), int)
