# Reopt Heuristic

This repository contains a modular implementation of a hybrid genetic algorithm for energy system optimization, originally developed for solving sizing and dispatch problems involving electric and thermal energy systems.

## Features
- Hybrid GA with tournament, roulette wheel, and hybrid selection
- Latin Hypercube initialization
- Warm-start and cold-start AMPL integrations
- Battery dispatch via linear programming (including thermal model)
- Modular and extensible architecture for future improvements

## Repository Structure
```
reopt_heuristic/
├── core/               # Core models and system logic
├── data/               # Twelves case study examples in pickled files
├── heuristics/         # Genetic algorithm logic
├── optimization/       # AMPL interface and output modules
├── packages/           # Required due to unpickling of older files
├── utils/              # Helper utilities, seeding, timing
├── main.py             # Source for program
├── requirements.txt    # Python dependencies
├── .gitignore          # Ignored files and folders
```

## Setup
1. **Clone the repo**
```bash
git clone https://github.com/YOUR_USERNAME/reopt_heuristic.git
cd reopt_heuristic
```

2. **Install dependencies**
```bash
pip install -r requirements.txt
```

3. **Install AMPL**
You will need to create an account and obtain a license you can then run the following commands to activate.

```bash
# Install Python API for AMPL:
$ python -m pip install amplpy --upgrade

# Install solver modules:
$ python -m amplpy.modules install highs gurobi

# Activate your license (e.g., free ampl.com/ce or ampl.com/courses licenses):
$ python -m amplpy.modules activate <your-license-uuid>
```

4. **Run the main script**
```bash
python main.py
```

## Example Usage
To run the genetic algorithm on a given problem instance:
```python
from heuristics import run_heur
from core import Top

model = Top(case_file="example_case.dat")
best_solution = run_heur(model, num_iter=30, pop_size=25, trial=1, selection="tournament")
```

## Requirements
- Python 3.8+
- AMPL (with Gurobi or GLPK solver)

## Contributors
Developed by Jamie Grymes and collaborators.

## Acknowledgments
This work was supported by research into hybrid energy systems with a special thanks to Dr. Alexander Zolan, Dr Alexxandra Newman, and Dr. Dinesh Mehta along with support from the National Renewable Energy Labratory
