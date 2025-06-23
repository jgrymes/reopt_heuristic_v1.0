# File: heuristics/ga_runner.py
# Purpose: Wrapper functions for running the GA heuristic with or without time limit
import sys
from core.system_model import Top
from heuristics.ga_utils import run_ga


def run_heur(case_name, num_iter, batt, pop_size, mut, trial, selection="org", file='log_chatter.txt'):
    _redirect_stdout(file)
    print(case_name, selection)

    case = f'{case_name}_1.pkl'
    d = Top(case, percent=0.3, batt=batt)
    d.mut = mut

    sol, d = run_ga(d, num_iter, pop_size, trial, selection)
    _restore_stdout()
    return sol, d





# --- Internal Helper Functions ---

_orig_stdout = sys.stdout

def _redirect_stdout(filename):
    global _orig_stdout
    sys.stdout = open(filename, 'a')

def _restore_stdout():
    global _orig_stdout
    sys.stdout.close()
    sys.stdout = _orig_stdout