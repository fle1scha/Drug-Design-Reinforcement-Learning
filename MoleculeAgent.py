import gym
from gym import wrappers, logger

class MoleculeAgent:
    """
    A class used to represent an Agent that implements a RL policy.

    ...
    Methods
    ----------
    __init__(action_space)
        the constructor method for a MoleculeAgent object.

    act(observation, reward, done)
        returns the chosen action based on the Agent's policy.
    """

    def __init__(self, action_space):
        """
        Parameters
        ----------
        action_space : spaces.Discrete(6)
            The actions available to the agent, represented numerically.
        """

        self.action_space = action_space
  
    def act(self, observation, reward, done):
        """Determines the MoleculeAgent's action to return to the environment.

        Parameters
        __________
        observation : 

        """
        #determine simply policy
        return self.action_space.sample()

