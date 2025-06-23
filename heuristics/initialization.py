# File: heuristics/initialization.py
# Purpose: Functions for creating initial population of GA individuals

import numpy as np
from random import randrange
from scipy.stats import qmc


def create_size(d):
    bubar = {}
    for c in d.C:
        for t in d.Tc[c]:
            bubar[t] = d.bubar_sigma[c]

    person = {}
    for t in list(set(d.T_CHP).union(set(d.T_td))):
        person[t[0:2]] = max(
            int(bubar[t]),
            randrange(
                int(bubar[t]),
                int(max(d.delta_d)),
                max(1, int((int(max(d.delta_d)) - int(bubar[t])) / 20))
            )
        )
        person[t[0:2]] = min(int(d.bbar_sigma[t]), person[t[0:2]])

    for b in d.B:
        person[b] = max(
            int(d.wubar_bkWh[b]),
            randrange(
                int(d.wubar_bkWh[b]),
                int(max(d.delta_d)),
                max(1, int((int(max(d.delta_d)) - int(d.wubar_bkWh[b])) / 20))
            )
        )
        person[b] = min(int(d.wbar_bkWh[b]), person[b])

    return person


def latin_hyper(d, pop_size):

    bubar=dict()
    for c in d.C:
        for t in d.Tc[c]:
            bubar[t]=d.bubar_sigma[c]
    
    techs=[]
    sizes=dict()
    for t in list(set(d.T_CHP).union(set(d.T_td))):
        techs+=[t[0:2]]
        if t[0:2]=='PV':
            sizes[t[0:2]]=[bubar[t],min(max(d.delta_d/.5),d.bbar_sigma[t])]
        else:    
            sizes[t[0:2]]=[bubar[t],d.bbar_sigma[t]] #pv upper bound might be too high
    for b in d.B:#.B_e
        techs+=[b]
        sizes[b]=[d.wubar_bkWh[b],min(max(d.delta_d),d.wbar_bkWh[b])]

    techs=list(sizes.keys())
    sampler = qmc.LatinHypercube(d=len(sizes)+1)
    sample = sampler.random(n=pop_size-1)
    pop=[]
    for i in range(pop_size-1):
        pop+=[dict([(techs[j],int(sizes[techs[j]][0]+(sizes[techs[j]][1]-sizes[techs[j]][0])*sample[i][j])) for j in range(len(techs))])]
    pop+=[dict([(techs[j],sizes[techs[j]][0]) for j in range(len(techs))])]
    # pop1=[]
    for i in range(pop_size-1):
        for b in d.B_e:
            pop[i][b*10]=.2+.4*sample[i][-1]
    for b in d.B_e:
        pop[-1][b*10]=.2+.4*sample[-1][-1]
    return pop
