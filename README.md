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
├── heuristics/         # Genetic algorithm logic
├── optimization/       # AMPL interface and output modules
├── utils/              # Helper utilities, seeding, timing
├── model/              # AMPL export/import handlers
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

3. **Run the main script**
```bash
python core/run_heuristic.py
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

## License
MIT License (add actual license if needed)

## Contributors
Developed by Jamie Grymes and collaborators.

## Acknowledgments
This work was supported by research into hybrid energy systems and military applications of optimization. Special thanks to contributors from West Point and OptTek Systems.
