# File: optimization/ampl_interface.py
# Purpose: Interface functions for calling AMPL models and processing feasibility/warm start runs

from amplpy import AMPL
from optimization.ampl_output import outputALL
import os


def check_feas(d, case):
    outputALL(d, "warm_now.dat")
    curr = os.getcwd()
    os.chdir("/home/jgrymes/reopt_heuristic/reopt")

    ampl = AMPL()
    ampl.eval(f"param case symbolic default '{case[0:6]}';")
    ampl.read("REoptLT.run")
    ampl.read("check_feas.run")
    os.chdir(curr)
    ampl.eval("option presolve_eps 0.000001;")
    ampl.eval("option presolve 0;")
    ampl.eval(f'option solver gurobi; option gurobi_options "writesol {case}_solution.sol writeprob {case}_model.mps timelim 300 mipgap 1e-4 outlev 0 bestbound 1 logfile check_fesibility_output.txt";')
    ampl.solve()
    return ampl


def output_model(d, case):
    outputALL(d, "warm_now.dat")
    curr = os.getcwd()
    os.chdir("/home/jgrymes/reopt_heuristic/reopt")

    ampl = AMPL()
    ampl.eval(f"param case symbolic default '{case[0:6]}';")
    ampl.read("REoptLT.run")
    ampl.read("solve_warm_start.run")
    os.chdir(curr)
    ampl.eval("option presolve_eps 0.000001;")
    ampl.eval("option presolve 0;")
    ampl.eval(f'option solver gurobi; option gurobi_options "writesol {case}_solution.sol writeprob {case}_model.mps timelim 1 mipgap 1e-4 outlev 0 bestbound 1 logfile check_fesibility_output.txt";')
    ampl.solve()
    return ampl


def solve_warm_start(d, case):
    outputALL(d, "warm_now.dat")
    curr = os.getcwd()
    os.chdir("/home/jgrymes/reopt_heuristic/reopt")

    ampl = AMPL()
    ampl.eval(f"param case symbolic default '{case[0:6]}';")
    ampl.read("REoptLT.run")
    ampl.read("solve_warm_start.run")
    os.chdir(curr)
    ampl.eval("option presolve_eps 0.000001;")
    ampl.eval("option presolve 1;")
    ampl.eval(f'option solver gurobi; option gurobi_options "timelim 20 mipgap 1e-4 outlev 1 bestbound 1 logfile ../revisions_rd_2/code2/logfiles/gurobi/{case[0:6]}_warm.log";')
    ampl.solve()
    return ampl


def solve_basic(case):
    curr = os.getcwd()
    os.chdir("/home/jgrymes/reopt_heuristic/reopt")

    ampl = AMPL()
    ampl.eval(f"param case symbolic default '{case[0:6]}';")
    ampl.read("REoptLT.run")
    os.chdir(curr)
    ampl.eval("option presolve_eps 0.000001;")
    ampl.eval("option presolve 1;")
    ampl.eval(f'option solver gurobi; option gurobi_options "timelim 300 mipgap 1e-4 outlev 1 bestbound 1 logfile ../revisions_rd_2/code2/logfiles/gurobi/{case[0:6]}_full.log";')
    ampl.solve()
    return ampl
