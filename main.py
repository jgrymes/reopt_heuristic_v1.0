# # Directory: reopt_heuristic
# # File: main.py
# # Purpose: Entry point for running heuristic optimization across multiple cases

# from heuristics.ga_runner import run_heur
# from optimization.ampl_interface import check_feas
# import multiprocessing
# import numpy as np

# cases = ['cdee0b']#,'e06d54','a60a2d','08f599','73506c','5407d0','9c42ca','559dac','bf6213','c7fbe7','960bfe','c2384e']
# cpu = multiprocessing.cpu_count()
# population_size = int(np.floor(cpu * 0.9))
# population_size=3
# selection_method = "ts"  # Tournament selection

# for case in cases:
#     for trial_num in range(1):
#         # if __name__ == '__main__':
#     # run_heur("case", num_iter=50, batt="lp", pop_size=30, mut=0.2, trial=1, selection="hybrid")
#         sol, d = run_heur(
#             case_name=case,
#             num_iter=2,
#             batt='alp',
#             pop_size=population_size,
#             mut=0.2,
#             trial=trial_num + 1,
#             selection=selection_method
#         )
#         check_feas(d, case)
        
        
# Directory: reopt_heuristic
# File: main.py
# Purpose: Entry point for running heuristic optimization across multiple cases

from heuristics.ga_runner import run_heur
from optimization.ampl_interface import check_feas
import multiprocessing
import numpy as np

def main():
    cases = ['9c42ca']#,'559dac','cdee0b']  # Add more cases as needed
    cpu = multiprocessing.cpu_count()
    population_size = int(np.floor(cpu * 0.9))
    print(population_size)
    # population_size = 9
    num_iter=2
    selection_method = "ts"  # Tournament selection

    for case in cases:
        for trial_num in range(1):
            sol, d = run_heur(
                case_name=case,
                num_iter=num_iter,
                batt='alp',
                pop_size=population_size,
                mut=0.2,
                trial=trial_num + 1,
                selection=selection_method
            )
            # check_feas(d, case)

if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")  # Ensures compatibility on Windows
    main()