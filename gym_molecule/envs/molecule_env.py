import gym
from gym import error, spaces, utils
from gym.utils import seeding


class MoleculeEnvironment(gym.Env):
    def __init__(self):
        print("works")
        super().__init__()
<<<<<<< HEAD
        
        self.action_space = spaces.Discrete( 3 )
        self.reset()
=======
        #We need to define the action space. It is a spaces.Discrete with actions 0..n-1. 
        self.action_space = spaces.Discrete("""Number of actions:""" 2) 
        
        self.observation_space = spaces.Box()
        print("Molecule Environment initiated.")
>>>>>>> c5e9137... Defined action and obs space.

    def step(self, action):
        print(action)
        
        """

        Parameters
        ----------
        action : {Option from policy}

        Returns
        -------
        Observation, reward, Done, info : tuple
             Observation (object) :
                 Observation of state of the environment.
                
            reward (float) :
                 Amount of reward associated with the previous action, goal is always to increase
                 your total reward.
                
            Done (bool) :
                 Whether it's time to reset the environment again.
                
            info (dict) :
                 Could contain the raw probabilities behind the environment's last state change.
                 
        """

    def reset(self):
        
        raise NotImplementedError

    def render(self):
        raise NotImplementedError

    def seed(self):
        raise NotImplementedError
