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


def build_ampl_model(d):
    import uuid

    unique_filename = str(uuid.uuid4())
    ampl=AMPL()
    ampl.read("/home/jgrymes/reopt_heuristic/code/packages/battery_lp.mod")
    ampl.read_data("/home/jgrymes/reopt_heuristic/code/packages/batt_data.dat")
    T=[]
    T_td=[]
    T_CHP=[]
    for t in d.T_e:
        if t in d.T_CHP.to_list():
            if d.X_sigma[t]>0:
                T+=[t]
                T_CHP+=[t]
        if t in d.T_td.to_list():
            if d.X_sigma[t]>0:
                T+=[t]
                T_td+=[t]
    ampl=set_set(ampl,"T",T)
    ampl=set_set(ampl,"T_td",T_td)
    ampl=set_set(ampl,"T_CHP",T_CHP)
    ampl=set_set(ampl,"U_P",d.U_P)
    ampl=set_set(ampl,"R",d.R)
    ampl=set_param(ampl,'eta_gplus',d.eta_gplus[1])
    ampl=set_param(ampl,'elec_load_unmet',d.elec_load_unmet)
    
    for b in d.B_e:
        ampl=set_param(ampl,'xbkwh',d.X_bkWh[b])
        ampl=set_param(ampl,'w_zero',d.w_zero[b])
        ampl=set_param(ampl,'wubar_mcp',d.wubar_mcp[b])
        ampl=set_param(ampl,'eta_minus',d.eta_minus[b])
    eta_plus=dict()
    prod_remain=dict()
    for t in d.T_CHP:
        for h in d.H:
            prod_remain[(t,h+1)]=max(0,d.Z_to[t][h]*(d.X_sigma[t]-d.xrp[t][h]))
    for t in T_td:
        for h in d.H:
            prod_remain[(t,h+1)]=d.f_l[t]*d.f_p[t][h]*d.X_sigma[t]#max(0,d.extra_renew[t][h])
    for b in d.B_e:
        for t in T:
            eta_plus[t]=d.eta_plus[b][t] 


    ampl=set_param(ampl,'eta_plus',eta_plus)
    ampl=set_param(ampl,'prod_remain',prod_remain)

    kappa_d=dict()
    for m in d.M:
        kappa_d[m]=(1-d.f_tot)*d.f_e[1]*d.c_rm[m][1]
    # for u in d.U_P:
    #     print(d.c_g[u])
    kappa_r=dict()
    for r in d.R:
        if d.c_r is None:
            kappa_r[r]=0    
        else:
            kappa_r[r]=(1-d.f_tot)*d.f_e[1]*d.c_r[r][1]
    
    kappa_c=dict()
    for t in d.T_CHP:
        for h in d.H:
            kappa_c[(t,h+1)]=d.chp_fuel_slope[h]*d.chp_fuel_cost
    # print(kappa_r)
    ampl=set_param(ampl,'kappa_e',(1-d.f_tot)*d.f_e[1]*np.array(d.c_g[1]))
    ampl=set_param(ampl,'kappa_d',kappa_d)
    ampl=set_param(ampl,'kappa_c',kappa_c)
    ampl=set_param(ampl,'kappa_r',kappa_r)
    
    ampl.eval(f"option solver gurobi; option gurobi_options 'outlev 0 method 1 threads 1  timing 1 logfile dual3.txt';")#option solver_msg 0; option log_file ''; 
    return ampl


def build_ampl_model_therm(d):
    import uuid

    unique_filename = str(uuid.uuid4())
    ampl=AMPL()
    ampl.read("/home/jgrymes/reopt_heuristic/code/packages/battery_lp.mod")
    ampl.read_data("/home/jgrymes/reopt_heuristic/code/packages/batt_data.dat")
    T=[]
    T_td=[]
    T_CHP=[]
    T_f=[]
    T_th=[]
    T_e=[]
    Tb=dict()
    B=[]
    B_e=[]
    T_ht=d.T_ht.values
    load=dict()
    eta_plus=dict()
    prod_remain_e=dict()
    prod_remain_th=dict()
    for b in d.B_e:
        Tb[b]=[]
        # if d.X_bkWh[b]>0:
        load[b]=d.elec_load_unmet
        B+=[b]
        B_e+=[b]    
        for t in d.Tb[b]:
            # if d.X_sigma[t]>0:
            Tb[b]+=[t]
            eta_plus[(b,t)]=d.eta_plus[b][t] 
            T+=[t]
            T_e+=[t]
            if t in d.T_CHP.to_list():
                T_CHP+=[t]
                T_f+=[t]
                chp_heat_slope=d.chp_heat_slope
                for h in d.H:
                    prod_remain_e[(t,h+1)]=max(0,d.Z_to[t][h]*(d.X_sigma[t]-d.xrp[t][h]))
                
            if t in d.T_td.to_list():
                T_td+=[t]
                for h in d.H:
                    prod_remain_e[(t,h+1)]=d.f_l[t]*d.f_p[t][h]*d.X_sigma[t]
    for b in d.B_h:
        B_h=[]
        Tb[b]=[]
        if d.X_bkWh[b]>0:
            B+=[b]
            B_h+=[b]
            load[b]=d.heat_load_unmet# d.delta_h*d.eta_b[1]
            max_load=max(load[b])
            for t in d.Tb[b]:
                # if d.X_sigma[t]>0:
                T+=[t]
                Tb[b]+=[t]
                T_th+=[t]
                if t in d.T_CHP.to_list():
                    for h in d.H:
                        prod_remain_th[(t,h+1)]=prod_remain_e[(t,h+1)]*chp_heat_slope[h]
                # if t in d.T_ac.to_list():
                else:
                    T_f+=[t]
                    for h in d.H:
                            prod_remain_th[(t,h+1)]=np.minimum(max_load,np.floor(d.bbar_sigma[t]))#max_load#d.X_sigma[t]
                    # if t in d.T_ec.to_list():
                    #     for h in d.H:
                    #             prod_remain_th[(t,h+1)]=d.X_sigma[t]
    for b in d.B_c:
        Tb[b]=[]
        T_ec=[]
        T_ac=[]
        B_c=[]
        if d.X_bkWh[b]>0:
            B+=[b]
            B_c+=[b]
            load[b]=d.delta_c*d.eta_ec
            max_load=max(load[b])
            for t in d.Tb[b]:
                if t in d.T_ac.to_list():
                    T_ac+=[t]
                if t in d.T_ec.to_list():
                    T_ec+=[t]
                # if d.X_sigma[t]>0:
                Tb[b]+=[t]
                T+=[t]
                T_th+=[t]
                for h in d.H:
                    prod_remain_th[(t,h+1)]=np.minimum(max_load,np.floor(d.bbar_sigma[t]))
    # for f in d.F:
    #     T_f+=[d.Tf[f][1]]      
    ampl=set_set(ampl,"T",set(T))
    ampl=set_set(ampl,"T_td",T_td)
    ampl=set_set(ampl,"T_CHP",T_CHP)
    ampl=set_set(ampl,"T_f",T_f)
    ampl=set_set(ampl,"T_e",T_e)
    ampl=set_set(ampl,"T_th",T_th)
    ampl=set_set(ampl,"T_ac",T_ac)
    ampl=set_set(ampl,"T_ec",T_ec)
    ampl=set_set(ampl,"T_ht",T_ht)
    ampl.display("T_ht")
    ampl=set_set(ampl,"B",B)
    ampl=set_set(ampl,"B_e",B_e)
    ampl=set_set(ampl,"B_c",B_c)
    ampl=set_set(ampl,"B_h",B_h)
    for b in B:
        ampl.set['Tb'][b]=Tb[b]
    ampl=set_set(ampl,"U_P",d.U_P)
    ampl=set_set(ampl,"R",d.R)
    ampl=set_param(ampl,'eta_gplus',d.eta_gplus[1])
    ampl=set_param(ampl,'eta_ec',d.eta_ec)
    ampl=set_param(ampl,'eta_ac',d.eta_ac)
    ampl=set_param(ampl,'chp_heat_slope',d.chp_heat_slope)
    ampl=set_param(ampl,'chp_xtra',d.curtailed_heat)
    
    loads=dict()
    for key,val in load.items():
        for i in range(len(val)):
            loads[(key,i+1)]=val[i]
    ampl=set_param(ampl,'load',loads)
    
    eta_plus=dict()
    for b in B:
        ampl.eval(f'let xbkwh[{b}]:={d.X_bkWh[b]};')
        ampl.eval(f'let w_zero[{b}]:={d.w_zero[b]};')
        ampl.eval(f'let wubar_mcp[{b}]:={d.wubar_mcp[b]};')
        ampl.eval(f'let eta_minus[{b}]:={d.eta_minus[b]};')
        ampl.eval(f'let w_d[{b}]:={d.w_d[b]};')
    
    for b in B:
        for t in Tb[b]:
            eta_plus[(b,t)]=d.eta_plus[b][t] 
    #     for h in range(len(d.remaining_prod[t])):
    #         prod_remain[(t,h+1)]=d.remaining_prod[t][h]
    ampl=set_param(ampl,'eta_plus',eta_plus)
    ampl=set_param(ampl,'prod_remain_e',prod_remain_e)

    ampl=set_param(ampl,'prod_remain_th',prod_remain_th)

    kappa_d=dict()
    for m in d.M:
        kappa_d[m]=(1-d.f_tot)*d.f_e[1]*d.c_rm[m][1]
    # for u in d.U_P:
    #     print(d.c_g[u])
    kappa_r=dict()
    for r in d.R:
        if d.c_r is None:
            kappa_r[r]=0    
        else:
            kappa_r[r]=(1-d.f_tot)*d.f_e[1]*d.c_r[r][1]
    
    kappa_c=dict()
    for t in T_f:#list(set(d.F_t)):
        if t in set(d.T_CHP):
            for h in d.H:
                kappa_c[(t,h+1)]=d.chp_fuel_slope[h]*d.chp_fuel_cost
        else:
            for h in d.H:
                kappa_c[(t,h+1)]=d.m_fm[t]
                
    # print(kappa_r)
    ampl=set_param(ampl,'kappa_e',(1-d.f_tot)*d.f_e[1]*np.array(d.c_g[1]))
    ampl=set_param(ampl,'kappa_d',kappa_d)
    ampl=set_param(ampl,'kappa_c',kappa_c)
    ampl=set_param(ampl,'kappa_r',kappa_r)
    
    ampl.eval(f"option solver gurobi; option gurobi_options 'outlev 1 method 1 threads 1 timing 1 logfile dual1.txt';")#option solver_msg 0; option log_file ''; 
    print(T_f,T,T_th,T_e)
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


def update_ampl_solution(dispatch,ampl):
    b =dispatch.B_e.values[0]
    dispatch.elec_load_unmet=rt_np_ampl(ampl,'xg')
    dispatch.xse[b]=rt_np_ampl(ampl,'xsc')
    
    # temp=[i[1] for i in ampl.get_variable('xsc').get_values().to_list()]
    # for h in range(8761):
    #     print(h,dispatch.xse[h],temp[h])
    dispatch.dfs[b]=rt_np_ampl(ampl,'dfs') 
    dispatch.gts=rt_np_ampl(ampl,'gts')
    # ampl.display("xsc,dfs,gts,xg")
    for u in dispatch.U_sb:
        if u == 3:
            stg=rt_np_ampl(ampl,'stg')
            for h in dispatch.H:
                dispatch.X_stg[u,h]=stg[h]
    pts=rt_pts_ampl(ampl)

    if len(dispatch.T_CHP)>0:
        for h in dispatch.H:
            x_rp=pts['CHP'][h]
            if x_rp>0:
                t='CHP'
                
                #Energy
                dispatch.xrp[t][h]+=x_rp
                val=min(dispatch.elec_load_unmet[h],x_rp)
                # print(h,dispatch.elec_load_unmet[h],x_rp)
                dispatch.extra_chp[h]+=x_rp-val
                dispatch.extra_elec[h]+=x_rp-val
                dispatch.remaining_prod[t][h]-=x_rp

                #Heat
                chp_x_tp=dispatch.chp_xtp(x_rp,h)
                heat_unmet=dispatch.heat_load_unmet[h]
                # if h in range(8739,8746):
                #     print(h,dispatch.X_tp[t][h],chp_x_tp,heat_unmet)
                            
                if chp_x_tp > heat_unmet:
                    # print(h,chp_x_tp, heat_unmet)
                    dispatch.X_tp[t][h]+= chp_x_tp
                    dispatch.heat_load_unmet[h] = 0
                    dispatch.Z_to[t][h]=1
                    dispatch.X_ptw[t][h] += chp_x_tp - heat_unmet
                    dispatch.curtailed_heat[h] += chp_x_tp - heat_unmet
                else:
                    dispatch.heat_load_unmet[h] -= chp_x_tp
                    dispatch.X_tp[t][h] += chp_x_tp
                    dispatch.Z_to[t][h]=1
                # a=dispatch.X_f[t][h]
                dispatch.update_fuel(t,h)
                # b=dispatch.X_f[t][h]
    T=ampl.get_set("T").get_values().to_list()

    
    for t in T:
        for h in dispatch.H:
            dispatch.pts[(b,t)][h]=pts[t][h]
       
    return dispatch
    
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
                d.xse[b][h]=xsc[b,h]
                d.dfs[b][h]=dfs[b,h+1]
    
    pts_e=return_data(ampl,'pts_e')
    for b in d.B_e:
        if d.X_bkWh[b]>0:
            for t in d.Tb[b]:
                if d. X_sigma[t]>0:
                    for h in d.H:
                        if t == 'CHP':
                            d.xrp[t][h]+=pts_e[t,h+1]
                            d.X_f[t][h]+=d.chp_fuel_slope[h]*pts_e[t,h+1]
                        d.pts[(b,t)][h]=pts_e[t,h+1]
    
    xtp=return_data(ampl,'xtp')
    pts_th=return_data(ampl,'pts_th')
    pts_extra=return_data(ampl,'pts_extra')
    for b in d.B_th:
        if d.X_bkWh[b]>0:
            for t in d.Tb[b]:
                # if d. X_sigma[t]>0:
                for h in d.H:
                    d.pts[(b,t)][h]=pts_th[t,h+1]
                    d.X_tp[t][h]=xtp[t,h+1]
                    if t =='CHP':
                        d.pts[(b,t)][h]+=pts_extra[t,h+1]
    d.gts=return_data(ampl,'gts')
    d.elec_load_unmet=return_data(ampl,'xg')
    
    
       
    return d