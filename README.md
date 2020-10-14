```
              _                 _      
  /\/\   ___ | | ___  ___ _   _| | ___ 
 /    \ / _ \| |/ _ \/ __| | | | |/ _ \
/ /\/\ \ (_) | |  __/ (__| |_| | |  __/
\/    \/\___/|_|\___|\___|\__,_|_|\___|
   __           _                                      _   
  /__\ ____   _(_)_ __ ___  _ __  _ __ ___   ___ _ __ | |_ 
 /_\| '_ \ \ / / | '__/ _ \| '_ \| '_ ` _ \ / _ \ '_ \| __|
//__| | | \ V /| | | | (_) | | | | | | | | |  __/ | | | |_ 
\__/|_| |_|\_/ |_|_|  \___/|_| |_|_| |_| |_|\___|_| |_|\__|

```                                                        
![testing](https://github.com/robmacc/capstone-molecule-environment/workflows/testing/badge.svg)

# Training Environment for Drug Design 

## Introduction 
This Python training environment (TEDD) allows an agent - essentially a learning piece of software - to modify drugs in order to generate novel structures with prescribed optimisation goals. The agent uses Reinforcement Learning - a machine learning technique - to make decisions, modify a molecule, and to reach the provided optimisation goal. The environment uses the AIGym framework, and the RDKit Chemistry package. 


Although the AIGym environment is open source, the extensions of it that are made by this project are not. If you would like access to the repo, please contact the [team lead.](mailto:MLSLUK002@myuct.ac.za)

## Installation
Use the environment management system [Conda,](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
and the package manager [pip](https://pip.pypa.io/en/stable/) to install the environment and its dependencies. 

In the TEDD project directory on your local computer:
```bash
conda env create -f environment.yml
```
Activate environment:
```bash
conda activate mol-env
```

Once the environment is activated, install all dependencies: 
```bash
pip install -e .
```

## Usage
To run the TEDD environment:
```bash
python MainMolcule.py
```

The user will be presented with the option to set the initial state of the molecule.  For first time users, it is suggested that you input '1' to allow the program to select its initial state. This will allow you to become familiar with the environment. 
```bash
Welcome to TEDD, the Training Environment for Drug Discovery.
In order to run the environment, you will need to input a starting molecule, a target molecule and an optimisation goal.
Step 1) Choose starting molecule: 
Input ‘1’ to choose a CARBON molecule.
Input ‘2’ to choose a RANDOM molecule.

```


Next, the user will be presented with the option to set the optimisation target molecule. This is essentially what the learning agent will try and reach. Again, for first time users, it is suggested to input '2' and allow the program to demonstrate its functionality. 
```bash
Step 2) Choose target molecule: 
Input ‘1’ to SPECIFY a molecule.
Input ‘2’ to choose a RANDOM molecule.

```

Lastly, the user will be prompted to input a desired similarity value. This value is the optimisation goal of the environment, and the environment will stop iterating once the initial molecule reaches aforementioned similarity. For example, an inputted a value of '0.6', will cause the environment to terminate when the starting molecule is 60% similar to its target molecule. 
```bash
Step 3) Choose optimisation goal: 
Input a floating-point value between 0 and 1.

```

The environment will then run through its iterations, returning information on its state, the agent's action and the Reinforcement Algorithm. 
