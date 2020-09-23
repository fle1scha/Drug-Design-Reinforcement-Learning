import gym
from gym import wrappers, logger

class MoleculeAgent:
    def __init__(self, action_space):
        self.action_space = action_space
  
    def act(self, observation, reward, done):
        #determine simply policy
        return self.action_space.sample()

