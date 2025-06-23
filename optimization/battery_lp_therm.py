# File: core/battery_lp_therm.py
# Purpose: Extended LP model for battery dispatch including thermal/efficiency constraints (if any)

import pandas as pd
import numpy as np
from amplpy import AMPL

def ret_np(obj):
    return np.array(list(obj.get_values().values()))

def set_param(ampl,name,val):
    if np.isscalar(val):
        param=ampl.get_parameter(name)
        param.set(val)
        return ampl
    else:
        param=ampl.get_parameter(name)
        param.set_values(val)
        return ampl
    
def set_set(ampl,name,val):
    new_set=ampl.get_set(name).set_values(val)
    return ampl



def build_ampl_model_therm(d):
    ampl=AMPL()
    ampl.read("optimization/battery_lp_heat.mod")
    ampl.get_set('H').set_values(np.array(d.H)+1)#np.array(range(1,8761))
    ampl.get_set('M').set_values(d.M)#np.array(range(1,13))
    ampl.get_set('R').set_values(d.R)#np.array(range(1,30))
    for m in d.M:
        ampl.set['Hm'][m]=d.Hm[m]+1
    for r in d.R:
        ampl.set['Hr'][r]=d.Hr[r]+1
    ampl.get_set('B').set_values(d.B.values)#np.array(range(1,30))
    ampl.get_set('B_e').set_values(d.B_e.values)
    ampl.get_set('B_c').set_values(d.B_c.values)
    ampl.get_set('B_h').set_values(d.B_h.values)
    ampl.get_set('T').set_values(d.T.values)
    ampl.get_set('T_e').set_values(d.T_e.values)
    ampl.get_set('T_h').set_values(d.T_ht.values)
    ampl.get_set('T_c').set_values(d.T_cl.values)
    ampl.get_set('T_CHP').set_values(d.T_CHP.values)
    ampl.get_set('T_th').set_values(np.concatenate([d.T_ht.values,d.T_cl.values]))
    ampl.get_set('T_ec').set_values(d.T_ec.values)
    ampl.get_set('T_ac').set_values(d.T_ac.values)
    for b in d.B:
        ampl.set['Tb'][b]=d.Tb[b]
    ampl.get_parameter('w_d').set_values(d.w_d)
    ampl.get_parameter('xbkwh').set_values(d.X_bkWh)
    ampl.get_parameter('xbkw').set_values(d.X_bkW)
    ampl.get_parameter('w_zero').set_values(d.w_zero)
    ampl.get_parameter('wubar_mcp').set_values(d.wubar_mcp)
    ampl.get_parameter('eta_minus').set_values(d.eta_minus)
    ampl.get_parameter('load_e').set_values(d.delta_d)#np.maximum(0,d.delta_d-sum(d.f_p[t]*d.f_l[t]*d.X_sigma[t] for t in d.T_td)))
    if len(d.T_cl)>0:
        load_c=d.delta_c*d.eta_ec
        ampl.get_parameter('load_c').set_values(load_c)
    else:
        load_c=np.zeros(8760)
        ampl.get_parameter('load_c').set_values(np.zeros(8760))
    ampl.get_parameter('load_h').set_values(d.delta_h*d.eta_b[1])
    
    if len(d.T_CHP)>0:
        ampl.get_parameter('chp_heat_slope').set_values(d.chp_heat_slope)
        ampl.get_parameter('chp_heat_intercept').set_values(d.Z_to['CHP']*d.chp_heat_intercept)
    else:
        ampl.get_parameter('chp_heat_slope').set_values(np.zeros(len(d.H)))
        ampl.get_parameter('chp_heat_intercept').set_values(np.zeros(len(d.H)))
        min_td=np.zeros(8760)
        
    df=pd.DataFrame()
    for t in d.T_e:
        if t == 'CHP':
            df[t]=d.Z_to[t]*(d.X_sigma[t])#-d.xrp[t]
            min_td=d.Z_to[t]*d.X_sigma[t]*d.fubar_td[t]
        else:
            df[t]=d.f_p[t]*d.f_l[t]*d.X_sigma[t]
        
    df.index=np.array(d.H)+1
    ampl.get_parameter('cap_e').set_values(df.unstack())
    ampl.get_parameter('min_td').set_values(min_td)
    df=pd.DataFrame()
    for t in np.concatenate([d.T_ht.values,d.T_cl.values]):
        if t == 'CHP':
            # df[t]=d.Z_to[t]*(d.X_sigma[t]-d.xrp[t])
            df[t]=d.Z_to[t]*(d.X_sigma[t]*d.chp_heat_slope+d.chp_heat_intercept)
        else:
            if t=='BOILER':
                max_load=min(d.bbar_sigma[t],max(d.delta_h*d.eta_b[1] + max(load_c)))
            else:
                max_load=min(d.bbar_sigma[t],max(d.delta_c*d.eta_ec))
            df[t]=np.full(len(d.H),max_load)
    df.index=np.array(d.H)+1
    free=-np.minimum(0,d.delta_h*d.eta_b[1]-sum(d.X_tp[t] for t in d.T_CHP))
    ampl.get_parameter('cap_free').set_values(free)
    if len(np.concatenate([d.T_ht.values,d.T_cl.values]))>0:
        ampl.get_parameter('cap_th').set_values(df.unstack())
        
    kappa_d=dict()
    for m in d.M:
        kappa_d[m]=(1-d.f_tot)*d.f_e[1]*d.c_rm[m][1]
    kappa_r=dict()
    for r in d.R:
        if d.c_r is None:
            kappa_r[r]=0    
        else:
            kappa_r[r]=(1-d.f_tot)*d.f_e[1]*d.c_r[r][1]
    df=pd.DataFrame()
    for t in d.T_ht:#list(set(d.F_t)):
        if t in set(d.T_CHP):
            df[t]=d.chp_fuel_slope*d.chp_fuel_cost
        else:
            df[t]=d.boiler_cost_param
    df.index=np.array(d.H)+1
    ampl.get_parameter('kappa_e').set_values((1-d.f_tot)*d.f_e[1]*np.array(d.c_g[1]))
    ampl.get_parameter('kappa_d').set_values(kappa_d)
    ampl.get_parameter('kappa_r').set_values(kappa_r)
    if len(d.T_ht)>0:
        ampl.get_parameter('kappa_c').set_values(df.unstack())
    df=pd.DataFrame()
    for b in d.B:
        for t in d.Tb[b]:
                ampl.eval(f"let eta_plus[{b},'{t}']:={d.eta_plus[b,t]};")
    ampl.param['eta_gplus']=d.eta_gplus[1]
    if len(d.T_cl)==0:
        ampl.param['eta_ac']=0
        ampl.param['eta_ec']=0
    else:
        ampl.param['eta_ac']=d.eta_ac
        ampl.param['eta_ec']=d.eta_ec
    T_td=[t for t in d.T_td]
    ampl.get_set('T_td').set_values(T_td)
    if len(T_td)>0:
        renew=pd.DataFrame()
        for t in T_td:
            renew[t]=d.f_l[t]*d.f_p[t]*d.X_sigma[t]
        renew.index=np.array(d.H)+1
        ampl.get_parameter('renew').set_values(renew.unstack())


    #T_s
    T_s=[t for t in d.Tv[1] if d.X_sigma[t]>0]
    ampl.get_set('T_s').set_values(T_s)
    #sales_cap
    if len(T_s)>0:
        sales_cap=pd.DataFrame()
        for t in T_s:
            sales_cap[t]=d.f_p[t]*d.f_l[t]*d.X_sigma[t]
        sales_cap.index=np.array(d.H)+1
        ampl.get_parameter('sales_cap').set_values(sales_cap.unstack())
        #deltabar_tu
        ampl.param['deltabar_gs']=d.deltabar_gs[d.U_nm[1]]
    #kappa_s
    kappa_s=(1-d.f_tot)*d.f_e[1]*d.Delta[1]*sum(d.c_e[u] for u in d.U_nm)
    ampl.get_parameter('kappa_s').set_values(kappa_s)
    
    ampl.eval(f"option gentimes 0; option solver gurobi; option gurobi_options 'outlev 0 method 1 threads 1  timing 0';")
   
    return ampl

def rt_np_ampl(ampl,name):
    x=ampl.get_variable(name).get_values().to_list()
    x=np.array([i[1] for i in x])
    return x

def rt_pts_ampl(ampl):
    T=ampl.get_set("T").get_values().to_list()
    try:
        x=ampl.get_variable("pts").get_values().to_dict()
        pts=dict()
        for t in T:
            pts[t]=[]        
        for key,val in x.items():
            pts[key[0]]+=[val]
        for t in T:
            pts[t]=np.array(pts[t])
    except:
        pts=dict()
    return pts


    
def return_data(ampl,var_name):
    x=ampl.get_variable(var_name)    
    if x.indexarity()==2:
        return x.get_values().to_dict()
    else:
        # for key,val in x1.items():
        return x.get_values().to_pandas()[var_name+".val"].values
    

def update_ampl_solution_heat(d,ampl):
    xsc=return_data(ampl,'xsc')
    dfs=return_data(ampl,'dfs')
    for b in d.B:
        if d.X_bkWh[b]>0:
            d.xse[b][0]=xsc[b,0]
            for h in d.H:
                d.xse[b][h+1]=xsc[b,h+1]
                d.dfs[b][h]=dfs[b,h+1]
    
    pts_e=return_data(ampl,'pts_e')
    if len(d.T_CHP)>0:
        xrp=return_data(ampl,'xrp')
    
    
    T_s=[t for t in d.Tv[1] if d.X_sigma[t]>0]
    if len(T_s)>0:
        ptg=return_data(ampl,'ptg')
        for t in T_s:
            for h in d.H:
                d.ptg[t][h]=ptg[t,h+1]

    for b in d.B_e:
        # if d.X_bkWh[b]>0:
        for t in d.Tb[b]:
            if d. X_sigma[t]>0:
                for h in d.H:
                    if t == 'CHP':
                        d.xrp[t][h]=xrp[t,h+1]
                        d.X_f[t][h]=d.chp_fuel_slope[h]*xrp[t,h+1]+d.Z_to[t][h]*d.chp_fuel_intercept[h]
                    d.pts[(b,t)][h]=pts_e[t,h+1]
    if len(np.concatenate([d.T_ht.values,d.T_cl.values]))>0:
        xtp=return_data(ampl,'xtp')
        pts_th=return_data(ampl,'pts_th')
        # pts_extra=return_data(ampl,'xtp_f')
    for b in d.B_th:
        # if d.X_bkWh[b]>0:
        for t in d.Tb[b]:
            # if d. X_sigma[t]>0:
            for h in d.H:
                d.pts[(b,t)][h]=pts_th[t,h+1]
                if t =='CHP':
                    d.X_tp[t][h]=xtp[t,h+1]
                    # d.X_tp[t][h]+=pts_extra[h]
                else:
                    d.X_tp[t][h]=xtp[t,h+1]


    d.stg=return_data(ampl,'stg')
    d.gts=return_data(ampl,'gts')
    d.elec_load_unmet=return_data(ampl,'xg')
    return d