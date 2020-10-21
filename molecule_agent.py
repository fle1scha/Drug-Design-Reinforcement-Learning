import gym
from gym import wrappers, logger
from collections import deque

MEMORY_SIZE = 1000000


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

    def __init__(self, observation_space, action_space):
        """
        Parameters
        ----------
        action_space : spaces.Discrete(6)
            The actions available to the agent, represented numerically.
        """

        self.action_space = action_space
        self.memory = deque(maxlen=MEMORY_SIZE)

    def act(self, observation, reward, done):
        """Determines the MoleculeAgent's action to return to the environment.

        Parameters
        __________
        observation : list
            The current observation state of the environment.

        reward : int
            The environment's reward to the agent.

        done : Boolean
            Whether or not the environment has reached its optimisation goal.

        Returns
        __________
        action : int
            The chosen action of the agent, given its policy.
        """

        # determine simply policy
        action = self.action_space.sample()
        return action

    def remember(self, state, action, reward, next_state, done):
        """Used to add a previous iteration of the environment to memory.

        Parameters
        __________
        state : float
            The state of the environment.
        action : int
            The action chosen by the agent.
        reward : int
            The reward given to the agent.
        next_state : float
            The next state of the environment.
        done : boolean
            Whether the environment is terminal or not.
        """
        
        self.memory.append((state, action, reward, next_state, done))
        
    def save_memory(self):
        with open("memory.txt", "w+") as f:
            f.write('\n'.join('%s %s %s %s %s' % step for step in self.memory))
                
        


