# File: heuristics/ga_utils.py
# Purpose: Selection, fitness evaluation, crossover, mutation functions for the GA

import numpy as np
from random import random, uniform, randrange, sample
from copy import copy
from timeit import default_timer as timer
import pickle
import cloudpickle
import csv
from multiprocessing import Pool, get_context
from itertools import repeat
from heuristics.initialization import latin_hyper

def system_input(d,person):
    prod_size=dict()
    batt_size=dict()
    techs=list(set(d.T_CHP).union(set(d.T_td)))
    for t in d.T:
        if t in techs:
            prod_size[t]=person[t[0:2]] 
            person[t[0:2]]=0   
        else:
            prod_size[t]=0
    for b in d.B:
        batt_size[b]=person[b]

    for b in d.B_e:
        batt_size[b*10]=person[b*10]
    return (prod_size,batt_size)

def fitness1(d, person, pickled_mod):
    disp_time = timer()
    prod_size, batt_size = system_input(d, person)
    if d.batt == 'lp':
        model = pickle.loads(pickled_mod)
    else:
        model = ""
    d.init_system(prod_size, batt_size)
    d.run_dispatch(model)
    return d.lcc(), d.batt_timer


def parent_selection(fit_val, population):
    selected_parents = int(len(population) / 2)
    idx = sorted(range(len(fit_val)), key=lambda i: fit_val[i])[:selected_parents]
    parents = [population[i] for i in idx]
    curr_obj = fit_val[idx[0]]
    curr_sol = population[idx[0]]
    return parents, curr_obj, curr_sol


def parent_selection_roulette(fit_val, population):
    fit_val1 = [1 / i for i in fit_val]
    list1 = np.cumsum(fit_val1) / sum(fit_val1)
    parents = []
    new_fit = []
    for _ in range(int(len(population) / 2)):
        num = random()
        i = 0
        while num > list1[i]:
            i += 1
        parents.append(population[i])
        new_fit.append(fit_val[i])
    best = np.argmin(new_fit)
    curr_obj = fit_val[best]
    curr_sol = population[best]
    return parents, curr_obj, curr_sol


def parent_selection_tournament(fit_val, population):
    foo = [(i, j, k) for i, j, k in zip(range(len(fit_val)), fit_val, population)]
    parent_size = int(len(foo) / 2)
    k = 4
    parents = []
    for _ in range(parent_size):
        tournament = sample(foo, k)
        best_fit = float('inf')
        for idx, fit, gene in tournament:
            if fit < best_fit:
                best_fit, best_idx, best_gene = fit, idx, gene
        parents.append((best_idx, best_fit, best_gene))
    best = min(parents, key=lambda p: p[1])
    parents = [i[2] for i in parents]
    curr_obj = best[1]
    curr_sol = best[2]
    return parents, curr_obj, curr_sol


def hybrid_selection(fit_val, population, best_obj):
    if np.min(fit_val) <= best_obj:
        return parent_selection_tournament(fit_val, population)
    if uniform(0, 1) <= 0.7:
        return parent_selection_tournament(fit_val, population)
    return parent_selection_roulette(fit_val, population)


def mating(parents, children, bounds, mut):
    for i in range(len(parents) - 1):
        for j in range(i + 1, len(parents)):
            child1, child2 = copy(parents[i]), copy(parents[j])
            for t in child1.keys():
                if random() > 0.5:
                    child1[t] = parents[j][t]
                    child2[t] = parents[i][t]
            child1 = mutation1(child1, bounds, mut)
            child2 = mutation1(child2, bounds, mut)
            if child1 not in parents and child1 not in children:
                children.append(child1)
            if child2 not in parents and child2 not in children:
                children.append(child2)
    return children


def crossover(parents, pop_size, bounds, mut):
    children = []
    count = 0
    while len(children) <= pop_size:
        children = mating(parents, children, bounds, mut)
        count += 1
    return sample(children, pop_size)


def mutation1(child, bounds, mut):
    for t in child.keys():
        if (isinstance(t,str)):
            child_range=(bounds[t][1]-bounds[t][0])/bounds['iteration']
            if random()<mut:
                if child[t]==bounds[t][0]:
                    delta_x=uniform(0,child_range)
                elif child[t]==bounds[t][1]:
                    delta_x=uniform(-child_range,0)
                else:
                    delta_x=uniform(-child_range,child_range)
                child[t]+=delta_x
        else:
            if t > 9:
                child_range=.4/bounds['iteration']
                if random()<mut:
                    if child[t]==.2:
                        delta_x=uniform(0,child_range)
                    elif child[t]==.6:
                        delta_x=uniform(-child_range,0)
                    else:
                        delta_x=uniform(-child_range,child_range)
                    child[t]+=delta_x
            else:
                child_range=(bounds[t][1]-bounds[t][0])/bounds['iteration']
                if random()<mut:
                    if child[t]==bounds[t][0]:
                        delta_x=uniform(0,child_range)
                    elif child[t]==bounds[t][1]:
                        delta_x=uniform(-child_range,0)
                    else:
                        delta_x=uniform(-child_range,child_range)
                    child[t]+=delta_x
    return child



def series_process(d,population,models):
    fit_val=[]
    for i in range(len(population)):
        person=copy(population[i])
        fit_val+=[fitness1(d,person,models[i])]
    fit_vals=[i[0] for i in fit_val]
    d.batt_timers=fit_val[-1][1]
    # print(fit_vals,d.batt_timers)
    return fit_vals


def multi_process(d, population, models):
    pop_size = len(population)
    ctx = get_context("spawn")  # Use "spawn" to be safe on Windows
    with ctx.Pool(pop_size) as p:
        values = list(zip(repeat(d), population, models))
        fit_val = p.starmap(fitness1, values)
    
    fit_vals = [i[0] for i in fit_val]
    for i in fit_val:
        d.batt_timers += i[1]
    return fit_vals




def run_ga(d,num_iter,pop_size,trial,selection,run_type='parallel',opt_val=0):
    counter=0
    obj_vals=dict()
    d.it_timer=[]
    d.batt_obj=[]
    d.run_type=run_type
    best_obj=100000000000000 #Initializes best objective found
    population=latin_hyper(d,pop_size)

    bounds=dict()

    for c in d.C:
        for t in d.Tc[c]:
            bounds[t[0:2]]=(d.bubar_sigma[c],d.bbar_sigma[t])
    # print(bounds)
    for b in d.B:
        bounds[b]=(d.wubar_bkWh[b],d.wbar_bkWh[b])
            
    model=""
    start=timer()

    x = np.linspace(1, num_iter, num_iter)
    z =pop_size/(1- .98/(1 + np.exp(-x+int(num_iter*0.5))))
    total_time=0
    for i in range(num_iter): #loops over the number of generations provided by the user
        start_time=timer()
        bounds['iteration']=z[i]
        tot_iter=timer()
        start=timer()
        if d.batt=='lp':
            if __name__ == "__main__":
                fit_val = multi_process(d, population, models)

        else:
            models=[j for j in range(pop_size)]
            if run_type=='series':
                fit_val=series_process(d,population,models)
            else:
                fit_val = multi_process(d, population, models)
            
        pop_creation=timer()
        
        match selection:
            case "org":
                parent_info=parent_selection(fit_val,population) #Selects the best performing members of the population, returns a tuple of the selected parents, and best solution/obj_fn of all parents
            
            case "ts":    
                parent_info=parent_selection_tournament(fit_val,population) #Selects the best performing members of the population, returns a tuple of the selected parents, and best solution/obj_fn of all parents
        
            case "rw":
                parent_info=parent_selection_roulette(fit_val,population) #Selects the best performing members of the population, returns a tuple of the selected parents, and best solution/obj_fn of all parents        
                
            case "hybrid":
                parent_info=hybrid_selection(fit_val,population,best_obj)

            case _:
                parent_info=parent_selection(fit_val,population) #Selects the best performing members of the population, returns a tuple of the selected parents, and best solution/obj_fn of all parents

        parents=parent_info[0] 
        curr_obj=parent_info[1]
        curr_sol=parent_info[2]
        # print(parent_info)
        
        change=0
        if curr_obj<=best_obj: #Checks if there was an improvement on the best solution
            change=abs((curr_obj-best_obj)/curr_obj)    
            best_obj=curr_obj
            best_sol=curr_sol
        if change<=0.0005:
            counter+=1
            # print(f"no change {counter}")
        else:
            counter=0
        
        #####Keeps elitist solution in the population.
        #####Removing

        obj_vals[f"iteration_{i}"]=[population,fit_val,best_obj,best_sol] #Stores the generations solutions
        print(f"iteration_{i+1}: Best Obj: {best_obj}, Best Solution: {best_sol} and Time: {timer()-tot_iter}")
        d.it_timer+=[timer()-tot_iter]
        filename = "iteration_results.csv"
        with open(filename, 'a', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter=',')
            # csvwriter.writerows([[d.name[0:6],pop_size,d.mut,trial,f"iteration_{i+1}",best_obj,timer()-tot_iter]])
            csvwriter.writerows([[d.name[0:6],pop_size,d.mut,trial,f"iteration_{i+1}",best_obj,timer()-tot_iter,selection]])
        children=crossover(parents,pop_size,bounds,d.mut) #Creates the next set of children from the selected parents
        # children+=[best_sol]
        
                
        for child in children:
            for t in list(set(d.T_CHP).union(set(d.T_td))):
                t=t[0:2]
                if child[t]<bounds[t][0]:
                    child[t]=bounds[t][0]
                elif child[t]>bounds[t][1]:
                    child[t]=bounds[t][1]
                
            for b in d.B:
                if child[b]<bounds[b][0]:
                    child[b]=bounds[b][0]
                elif child[b]>bounds[b][1]:
                    child[b]=bounds[b][1]

            for b in d.B_e:
                if child[b*10]<0.2:
                    child[b*10]=0.2
                elif child[b*10>.6]:
                    child[b*10]=.6

        population=children

        if counter>=15: #terminates if the best does not improve after 5 iterations
            break
        total_time+=(timer()-start_time)
        print(total_time)
        if (total_time)>130:
            break

    prod_size,batt_size=system_input(d,best_sol)#population[i]
    d.init_system(prod_size,batt_size)
    d.run_dispatch(model)
    return (obj_vals,d)

