# File: utils/helper_function.py
# Purpose: Contains miscellaneous helper functions for GA and model evaluation workflows

import numpy as np
import pickle
from copy import copy
from timeit import default_timer as timer
import numpy as np
import pandas as pd
from packages.exportAMPL import loadSolution1 as ls
from copy import copy
import gzip
import pickle

def initialize(vals,name,percent,batt):#prod_sizes,bat_energy,bat_power,
    ####################################
    ###Do once during initialization####
    ####################################
    start=timer()
    sol=ls(f"data/{name}")#a.loadSolution(name)
    [setattr(vals, attr, getattr(sol,attr)) for attr in dir(sol) if attr.startswith('__') is False]
    del sol
    vals.name=name
    print(vals.T.to_list())  
    # b=vals.B_e.to_list()[0]
    # vals.Tb=vals.Tb[b].to_list()
    Tb=dict()
    for b in vals.B:
        Tb[b]=vals.Tb[b].to_list()
    vals.Tb=Tb
    ####General Params#####
    vals.batt_timer=[]
    vals.batt_timers=[]
    vals.it_timer=[]
    vals.batt_obj=[]
    vals.batt=batt
    vals.H=[h for h in range(0,8760)]
    Hm=dict()
    for m in vals.M:
        Hm[m]=vals.Hm[m].values.astype(dtype=np.int32)-1
    vals.Hm=Hm
    del Hm
    Hr=dict()
    for r in vals.R:
        Hr[r]=vals.Hr[r].values.astype(dtype=np.int32)-1
    vals.Hr=Hr
    del Hr
    vals.percent=percent
    vals.periods=8760
    fp=dict()
    for t in  vals.f_p.index.unique(0).to_list():
        fp[t]=vals.f_p[t].values
    vals.f_p=fp
    zto=dict()
    for t in  vals.Z_to.index.unique(0).to_list():
        zto[t]=np.zeros(8760)
    vals.Z_to=zto
    xptw=dict()
    for t in  vals.X_ptw.index.unique(0).to_list():
        xptw[t]=np.zeros(8760)
    vals.X_ptw=xptw
    xtp=dict()
    for t in vals.X_tp.index.unique(0).to_list():    
        xtp[t]=np.zeros(8760)
    vals.X_tp=xtp
    
    xf=dict()
    for t in vals.X_f.index.unique(0).to_list():    
        xf[t]=np.zeros(8760)
    vals.X_f=xf
    
    # xdfs=dict()
    # for t in vals.X_dfs.index.unique(0).to_list():    
    #     xdfs[t]=np.zeros(8760)
    # vals.X_dfs=xdfs
    fha=dict()
    fht=dict()
    ffa=dict()
    for t in vals.f_ha.index.unique(0).to_list():    
        fha[t]=vals.f_ha[t].values
        fht[t]=vals.f_ht[t].values
        ffa[t]=vals.f_fa[t].values
    vals.f_ha=fha
    vals.f_ht=fht
    vals.f_fa=ffa
    
    vals.f_tow=vals.f_tow[1]
    vals.f_tot=vals.f_tot[1]
    vals.f_om=vals.f_om[1]
    comp=dict()
    for t in vals.c_omp.index.unique(0).to_list():
        comp[t]=vals.c_omp[t]
    vals.c_omp=comp
    
    vals.Mh=np.full(vals.periods,1)
    for m in vals.M:
        for h in vals.Hm[m]:
            vals.Mh[h]=m
    vals.Rh=np.full(vals.periods,1)
    for r in vals.R:
        for h in vals.Hr[r]:
            vals.Rh[h]=r
            
            
    if len(vals.T_cl)>0:  
        for t in vals.T_cl:
            if t == 'ELECCHL':
                vals.eta_ec=vals.eta_ec[1]  
            if t == 'ABSORPCHL':
                vals.eta_ac=vals.eta_ac[1]        
        vals.cool_load_unmet=vals.delta_c*vals.eta_ec
        max_cl=max(vals.cool_load_unmet)
        if max_cl>0:
            for t in ['ELECCHL','ABSORPCHL']:
                if t in vals.T_cl.to_list():
                    if max_cl>0:
                        if max_cl<vals.bbar_sigma[t]:
                            vals.X_sigma[t]=np.ceil(max_cl)
                            max_cl=0
                        else:
                            vals.X_sigma[t]=np.floor(vals.bbar_sigma[t])
                            max_cl-=np.floor(vals.bbar_sigma[t])
    
        
    vals.delta_h=vals.delta_h.values
    vals.delta_d=vals.delta_d.values
    vals.delta_c=vals.delta_c.values
    vals.elec_load_unmet=copy(vals.delta_d)
    vals.heat_load_unmet=vals.delta_h*vals.eta_b[1]
    vals.curtailed_heat=np.full(8760,0,dtype=np.float64)
    
    #####Variables#####
    
    for c in vals.C:
        for t in vals.Tc[c]:
            vals.X_sigma[t]=vals.bubar_sigma[c]
        vals.xrp=dict()
    vals.dfs=dict() 
    w_d=dict()
    for b in vals.B:   
        vals.dfs[b]=np.zeros(8760)
        try:
            w_d[b]=vals.w_d[b]
        except:
            w_d[b]=0
    vals.w_d=w_d
    vals.gts=np.zeros(8760)
    vals.ptg=dict()
    for t in vals.T_e:
        vals.ptg[t]=np.zeros(8760)
    vals.xse=dict()
    for b in vals.B:
        vals.xse[b]=np.full(8761,vals.w_zero[b]*vals.X_bkWh[b])
    vals.remaining_prod=dict()
    for t in vals.X_rp.index.unique(0).to_list():    
        vals.xrp[t]=np.zeros(8760)
        vals.remaining_prod[t]=np.full(8760,vals.X_sigma[t],dtype='float')    
    vals.pts=dict()
    for b in vals.B:
        for t in vals.Tb[b]:
            vals.pts[(b,t)]=np.zeros(8760)
    for t in vals.Z_to.keys():        
        vals.Z_to[t]*=0
    for t in vals.X_tp.keys():  
        vals.X_tp[t]*=0
    for t in vals.X_f.keys():  
        vals.X_f[t]*=0
    for t in vals.Z_sigmas.keys():  
        vals.Z_sigmas[t]*=0
        
    ####Utility params
    c_g=copy(vals.c_g)
    vals.c_g=dict()
    for u in c_g.index.unique(0).to_list():
        vals.c_g[u]=np.array(c_g[u])
    del c_g
    c_e=copy(vals.c_e)
    vals.c_e=dict()
    for u in c_e.index.unique(0).to_list():
        vals.c_e[u]=np.array(c_e[u])
    del c_e
    vals.energy_cost_param=dict()
    for u in vals.U_P:
        vals.energy_cost_param[u]=(1-vals.f_tot)*vals.f_e[1]*vals.Delta[1]*vals.c_g[u]
    #create indicator array for on and off peak energy price. Below the mean is off-peak (0) and above the mean is on-peak (1)        
    average=np.mean(vals.c_g[1])#.mean()
    vals.peak_ind=np.zeros(8760)
    for h in vals.H:
        if vals.c_g[1][h]>average: vals.peak_ind[h]=1
    # vals.utility_sum=np.full(13,0)
    if vals.Tf is not None:
        vals.F_t=dict()
        for key,val in vals.Tf.to_dict().items():
            vals.F_t[val.to_list()[0]]=key[0]

    #####CHP Params######
    for t in vals.T_ht:
        if t =='CHP':
            kte=dict() 
            ktp=dict()
            kte[t]=vals.k_te[t]
            ktp[t]=vals.k_tp[t]
            vals.k_te=kte 
            vals.k_tp=ktp
            mfbm=dict()
            mfbm[t]=vals.m_fbm[t]
            vals.m_fbm=mfbm
            vals.extra_chp=np.zeros(8760)
            f=vals.F_t[t]
            vals.chp_om=(1-vals.f_tow)*vals.f_om*(vals.c_omp[t])
            vals.chp_fuel_slope=vals.f_fa[t]*vals.f_p[t]*vals.m_fm[t]
            vals.chp_fuel_cost=(1-vals.f_tot)*vals.Delta[1]*vals.c_u[f]*vals.f_pf[t]
            vals.td=vals.fubar_td[t]*vals.X_sigma[t]
            vals.chp_fuel_intercept=vals.f_fa[t]*vals.m_fbm[t]*vals.X_sigma[t]
            vals.chp_heat_slope=vals.f_ha[t]*vals.f_ht[t]*(vals.k_te[t])*vals.f_p[t]
            print(t)
            vals.chp_heat_intercept=vals.f_ha[t]*vals.f_ht[t]*vals.k_tp[t]*vals.X_sigma[t]
            # print(vals.X_sigma[t])
        if t =='BOILER':
            vals.boiler_cost_param=(1-vals.f_tot)*vals.c_u[f]*vals.f_pf[t]*vals.m_fm[t]
    
    vals.extra_elec=np.zeros(8760)
    vals.timer=dict()
    
    ####Renewables####
    if len(vals.T_td)>0:
        vals.extra_renew=dict()
        for t in vals.T_td:
            vals.extra_renew[t]=np.zeros(8760)
    
    for b in vals.B:
        vals.X_bkWh[b]=vals.wubar_bkWh[b]
        vals.X_bkW[b]=vals.X_bkWh[b]/2
    
    if vals.batt == 'alp':
        print("writing output")
        f = open("optimization/batt_data.dat", "w")
        for m in vals.M:
            f.write(f"set Hm[{int(m)}]:= ")
            for h in vals.Hm[m]:
                f.write(f"{h+1} ")
            f.write(f";\n")
        for r in vals.R:
            f.write(f"set Hr[{int(r)}]:= ")
            for h in vals.Hr[r]:
                f.write(f"{h+1} ")
            f.write(f";\n")
        f.close()
    
    
def print_time(self,start):
    stop=timer()
    print(round(stop-start,4))
    return stop

        
def fill_xg(self):
    delta_d=copy(self.elec_load_unmet)
    if len(self.U_P)>1:
        dbar_tu=np.array(self.deltabar_tu)
        dbar_cum=np.insert(np.array(self.deltabar_tu.cumsum()),0,0)
        for u in range(len(dbar_tu)):
            xg=np.zeros(self.periods)
            for h in range(self.periods):
                if delta_d[h]>dbar_cum[u+1]:
                    xg[h]=dbar_tu[u]
                else:
                    xg[h]=max(0,delta_d[h]-dbar_cum[u])
            self.X_g.loc[(u+1,slice(None))]=xg
    else:
        arrays=[np.full(8760,1),self.H]
        idx=pd.MultiIndex.from_arrays(arrays,names=('U_P','H'))
        self.X_g=pd.Series(data=delta_d,index=idx)
        
def order_demand(self,vals):
    return vals.sort_values(ascending=False).index            

def update_z_sigma(self):
    for t in self.T:
        if self.X_sigma[t]>0:
            self.Z_sigmas[t,1,1]=1
        else:
            self.Z_sigmas[t,1,1]=0

def update_fuel(self,t,h):
    if t == 'CHP':
        self.X_f[t,h]=self.f_fa[t,h]*(self.m_fbm[t]*self.X_sigma[t] + self.f_p[t,h]*self.m_fm[t]*self.X_rp[t,h])
    else:
        self.X_f[t,h]=self.m_fm[t]*self.X_tp[t,h]
    #TODO Consider making this a single routine done at the end of the algorithm


def loadSolution(filename):
    """
    Load a pickled solution from file
    
    :param filename: Name of a file containing a pickled solution 
    return 
    """    
    # open the file for reading
    with gzip.open(filename, 'rb') as f:
        solution = pickle.load(f)
        
    return solution


def saveSolution(solution,filename):
    with gzip.open(filename, 'wb') as f:
        pickle.dump(solution, f)
        

def c_s(d1,cap):
    xrp=np.empty(0)
    xg=np.empty(0) 
    target_m=dict()
    target_r=dict()    
    for m in d1.M:
        load=d1.elec_load_unmet[d1.Hm[m]]
        target_m[m]=max(0,max(load)-cap)
        xrp=np.concatenate((xrp,np.maximum(load-target_m[m],0)))
        xg=np.concatenate((xg,np.minimum(load,target_m[m])))
    for r in d1.R:
        if len(d1.Hr[r])>0:
            target_r[r]=max(xg[d1.Hr[r]])
        else:
            target_r[r]=max(xg)
    return target_m,target_r, xrp, xg
    
def monthly_tiers(d1,xg):
        x_g=dict()
        for u in d1.U_P:
            x_g[u]=np.zeros(len(d1.H))
        for m in d1.M:
            u=1
            limit=copy(d1.deltabar_tu[u])
            bucket=0
            for h in d1.Hm[m]:
                remain=limit-bucket
                load=xg[h]
                while load>0:
                    if load>remain:
                        x_g[u][h]+=remain
                        load-=remain
                        remain=0
                        u+=1
                        limit=d1.deltabar_tu[u]
                        bucket=0
                        remain=limit-bucket
                    else:
                        x_g[u][h]+= load
                        bucket+=load
                        load=0
        return x_g
def month_demand(d1,xg,xrp):
    xdn=dict()
    for m in d1.M:
        grid=xg[d1.Hm[m]]
        if d1.T_rdc.empty:
            max_val=max(grid)
        else:
            prod=sum(xrp[d1.Hm[m]] for t in d1.T_rdc)
            max_val=max(grid+prod)
        xdn[m]=np.zeros(len(d1.N))
        for n in d1.N:
            if max_val>d1.deltabar_mt[n]:
                xdn[m][int(n-1)]=d1.deltabar_mt[n]
                max_val-=d1.deltabar_mt[n]
            else:
                xdn[m][int(n-1)]=max_val
                break
    return xdn


def get_targets(d1,cap):
    target_m,target_r,xrp,xg = c_s(d1,cap)
    x_g=monthly_tiers(d1,xg)
    eng_cost=(1-d1.f_tot)*d1.f_e[1]*d1.Delta[1]*(sum(x_g[u][h]*d1.c_g[u][h] for u in d1.U_P for h in d1.H))
    xdn=month_demand(d1,xg,xrp)
    demand_cost=(1-d1.f_tot)*d1.f_e[1]*d1.Delta[1]*sum(d1.c_rm[m][n]*xdn[m][int(n-1)] for m in d1.M for n in d1.N)
    zto=[1 if i > 0 else 0 for i in xrp]
    fuel_cost=d1.chp_fuel_cost*(sum(d1.chp_fuel_intercept*zto)+sum(xrp*d1.chp_fuel_slope))
    total_cost=fuel_cost+demand_cost+eng_cost
    return total_cost,target_m,target_r