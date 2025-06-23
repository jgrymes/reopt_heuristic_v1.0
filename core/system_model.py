# File: core/system_model.py
# Purpose: Defines the main system model class and supporting functions for initializing and simulating system operation

import numpy as np
from copy import copy
from timeit import default_timer as timer
from optimization.battery_lp_therm import build_ampl_model_therm
from optimization.battery_lp_therm import update_ampl_solution_heat as update_ampl_solution_therm
from utils.helper_function import get_targets
from utils.helper_function import initialize

class Top:
    def __init__(self, name, percent, batt):
        s = timer()
        self = initialize(self, name, percent, batt)

    def update_system(self,prod_sizes,bat_energy):#,bat_power
        for b in self.B:
            val=bat_energy[b]
            if val > self.wbar_bkWh[b]:
                print(f"Battery {b} too large. Battery resized to {self.wbar_bkWh[b]}")
                val=self.wbar_bkWh[b]
            if val%2>0:
                val=max(0,val-1)
            self.X_bkWh[b]=val
            if b in self.B_e.values:
                self.X_bkW[b]=val*bat_energy[b*10]
            else:
                self.X_bkW[b]=val
                
        for t in self.T:
            self.X_sigma[t]=prod_sizes[t]
    ##############################
    ####Do at every iteration#####
    ##############################
    def init_system(self,prod_sizes,bat_energy):#,bat_power
        # self.it_timer=[]
        # self.batt_obj=[]
        s=timer()
        self.update_system(prod_sizes,bat_energy)
        self.check_NM()
        self.elec_load_unmet=copy(self.delta_d)
        self.extra_renew=dict()
        for t in self.T_td:
            self.extra_renew[t]=np.zeros(8760)
        self.heat_load_unmet=copy(self.delta_h*self.eta_b[1])
        if len(self.T_cl)>0:  
            self.cool_load_unmet=copy(self.delta_c*self.eta_ec)
        
        self.curtailed_heat=np.full(8760,0,dtype=np.float64)
        
        self.gts=np.zeros(8760)
        self.stg=np.zeros(8760)
        for b in self.B:
            self.xse[b]=np.full(8761,self.w_zero[b]*self.X_bkWh[b])
            self.dfs[b]*=0
        for t in self.xrp.keys():
            self.xrp[t]*=0
        for t in self.remaining_prod.keys():
            self.remaining_prod[t]=np.full(8760,self.X_sigma[t],dtype='float')
        for t in self.pts.keys():
            self.pts[t]*=0
        for t in self.Z_to.keys():        
            self.Z_to[t]*=0
        for t in self.X_tp.keys():  
            self.X_tp[t]*=0
        for t in self.X_f.keys():  
            self.X_f[t]*=0
        for t in self.Z_sigmas.keys():  
            self.Z_sigmas[t]*=0
        for t in self.X_ptw.keys():  
            self.X_ptw[t]*=0
        for m,n in self.X_dn.keys():  
            self.X_dn[m,n]*=0
        
        for t in self.T_CHP:
            self.td=self.fubar_td[t]*self.X_sigma[t]
            self.chp_fuel_intercept=self.f_fa[t]*self.m_fbm[t]*self.X_sigma[t]
            self.extra_elec=np.zeros(8760)
            self.chp_heat_intercept=self.f_ha[t]*self.f_ht[t]*self.k_tp[t]*self.X_sigma[t]
        # print("Initialization Complete!!")
        # print(f"Initialize full system: {timer()-s}")
    ##########
    
    
    def run_dispatch(self,batt_mod):
        self.timer=dict()
        timer_segs=['reduce_demand','storage','energy','boiler','run_chp','energy_boiler','energy_seg1','full_chp']
        for t in timer_segs:
            self.timer[t]=0
        #TODO Add an initialization function that checks for feasibility
            
        start=timer()
        #TODO confirm technologies otherwise don't need to do this
        if len(self.T_td)>0:
            self.dispatch_renew()
        stop=timer()
        # print(f"renewables: {stop-start}")
        
        start=timer()
        #TODO confirm cooling load otherwise don't need to do this
        if sum(self.delta_c)>0:
            self.initialize_cool()
        
        max_month=0
        for m in self.M:
            curr_max=sum(self.elec_load_unmet[self.Hm[m]])
            if curr_max > max_month:
                max_month=curr_max
        self.deltabar_tu[len(self.U_P)]=np.ceil(max_month)+sum(self.bbar_sigma[t]/self.eta_ec for t in self.T_ec)
        stop=timer()
        # print(f"cooling: {stop-start}")
        
        
        start=timer()
        #TODO Maybe consider checking if there are demand charges??
        if len(self.T_CHP)>0:
            self.set_demand_targets()
            stop=timer()
            
            start=timer()
            self.fuel_lookup()
            stop=timer()
            
            start=timer()
            self.dispatch(batt_mod)
            stop=timer()
        stor_ops=0
        if sum(self.X_bkWh[b] for b in self.B) > 0:
            self.storage_ops(batt_mod)
            stor_ops=1
        boiler=(self.delta_h*self.eta_b[1]+sum(self.X_tp[t1]/self.eta_ac for t1 in self.T_ac)+
                sum(self.pts[b,t] for t in self.T_ht for b in self.B_h))-(sum(self.X_tp[t] for t in self.T_ht) + 
                                                                          sum(self.dfs[b] for b in self.B_h))
        
        if stor_ops == 0:
            boil=np.where(boiler<0,0,boiler)
            self.X_tp['BOILER']=boil
            if len(self.T_CHP)>0:
                self.use_extra_heat()
        
        ptw=np.where(boiler>0,0,-boiler)
        if 'BOILER' in self.T.to_list():
            self.X_f['BOILER']=self.X_tp['BOILER']*self.m_fm['BOILER']
        for t,val in self.X_tp.items():
            if t in ['BOILER','ELECCHL','ABSORPCHL']:
                self.X_sigma[t]=max(val)
        if len(self.T_CHP)>0:
            self.X_ptw['CHP']=ptw
        
            
        self.fill_peaks()
        s=timer()
        max_month=0
        for m in self.M:
            curr_max=sum(self.elec_load_unmet[self.Hm[m]])
            if curr_max > max_month:
                max_month=curr_max
        self.deltabar_tu[len(self.U_P)]=np.ceil(max_month)
        self.monthly_demand()
                
    def fuel_lookup(self):
        t='CHP'
        self.fuel_conv=dict()
        intercept=self.f_fa[t]*self.m_fbm[t]*self.X_sigma[t]
        slope= self.f_fa[t]*self.f_p[t]*self.m_fm[t]
        self.fuel_conv[t] = tuple(zip(intercept,slope))
    
      
    ###########
    ###########
    #Renewables
    ###########
    ########### 

    def dispatch_renew(self):
        if self.T_td.empty: return
        else:
            for t in self.T_td:
                size=self.X_sigma[t]
                if size>0:
                    prod=np.full(8760,size)
                    self.xrp[t]=prod
                    self.remaining_prod[t]-=prod
                    val=self.f_p[t]*self.f_l[t]*prod
                    dispatch=np.minimum(val,self.elec_load_unmet)
                    # self.extra_elec=val-dispatch
                    self.extra_renew[t]+=val-dispatch#self.extra_elec
                    self.elec_load_unmet-=dispatch
                    self.Z_to[t]=np.full(8760,1)
       
    ###########
    ###########
    #Cooling Techs
    ###########
    ########### 
      
    def initialize_cool(self):
        cool_load_max=max(self.cool_load_unmet)
        # for t in self.T:
        #     print(t,self.bbar_sigma[t])
        if cool_load_max>0:
            eta_ec=self.eta_ec
            eta_ac=self.eta_ac
            tcl=self.T_cl.values

            t='ELECCHL'

            if t in tcl:
                fp=self.f_p[t]
                load=self.cool_load_unmet/fp
                dispatch=np.minimum(load,np.floor(self.bbar_sigma[t]))
                self.cool_load_unmet-=dispatch
                self.X_tp[t]+=dispatch
                self.X_sigma[t]=max(dispatch)
                self.elec_load_unmet+=dispatch/eta_ec
                _max=0

            if max(self.cool_load_unmet)>0:
                t='ABSORPCHL'
                if t in tcl:
                    fp=self.f_p[t]
                    load=self.cool_load_unmet/fp
                    if load.max()>self.bbar_sigma[t]:
                        print("can't meet cooling load")
                    else:
                        self.cool_load_unmet-=load
                        self.X_tp[t]+=load
                        self.X_sigma[t]=max(load)
                        self.heat_load_unmet+=load/eta_ac
            if max(self.cool_load_unmet)>0:
                print("Did not meet cooling load. Better Check. Could be bug or could be the upper bound on system size is too small.")
    
    ###########
    ###########
    #CHP Functions
    ###########
    ###########

    

    def dispatch(self,batt_mod):
        for m in self.M:
            u=1
            util_bin_max=self.deltabar_tu[u]
            util_sum=0
            energ_cost=self.energy_cost_param[u]
            for h in self.Hm[m]:
                start=timer()
                chp_remain=self.reduce_demand(h)     #First Green Box
                stop=timer()
                self.timer['reduce_demand']+=(stop-start)

                start=timer()
                if  (chp_remain>0):                         #Third Diamond
                    self.turn_on_up_CHP(h,chp_remain,energ_cost)       #Second and Third Green Box and Fourth Diamond
                unmet_heat=self.heat_load_unmet[h]
                stop=timer()
                self.timer['energy']+=(stop-start)

                util_sum+=self.elec_load_unmet[h]
                if util_sum > util_bin_max:
                    u+=1
                    util_bin_max=copy(self.deltabar_tu[u])
                    energ_cost=copy(self.energy_cost_param[u])
            
    
        
    def use_extra_heat(self):
        val=sum(self.X_tp[t] for t in self.T_CHP) - ((self.delta_h*self.eta_b[1]) + sum(self.X_tp[t1]/self.eta_ac for t1 in self.T_ac) +sum(self.X_tp[tb] for tb in self.T_ht if tb not in self.T_CHP.to_list()))
        val= np.where(val<0, 0, val)
        # val=copy(self.X_ptw['CHP'])
        if max(val)>0.0000001:
            if 'ABSORPCHL' in self.T_cl.to_list():
                
                xac=copy(self.X_tp['ABSORPCHL'])
                ab_chl_remain=copy(self.X_sigma['ABSORPCHL']-self.X_tp['ABSORPCHL'])
                xec=copy(self.X_tp['ELECCHL'])
                
                cool_load=self.delta_c*self.eta_ec
                a_chill_add=np.amin([self.eta_ac*val,cool_load,ab_chl_remain],axis=0)
                a_chill_reduce=np.maximum(0,xac+a_chill_add-cool_load)
                e_chill=np.maximum(0,cool_load-(xac+a_chill_add))
                e_chill_reduce=xec-e_chill
                self.X_tp['ABSORPCHL']+=a_chill_add-a_chill_reduce
                self.X_tp['ELECCHL']-=e_chill_reduce
                
                self.X_sigma['ABSORPCHL']=max(self.X_tp['ABSORPCHL'])
                self.heat_load_unmet-=a_chill_reduce/self.eta_ac
                self.X_sigma['ELECCHL']=max(self.X_tp['ELECCHL'])
                self.elec_load_unmet-=e_chill_reduce/self.eta_ec
                self.elec_load_unmet=np.maximum(self.elec_load_unmet,self.gts)
                self.X_ptw['CHP']=sum(self.X_tp[t] for t in self.T_CHP) - ((self.delta_h*self.eta_b[1]) + sum(self.X_tp[t1]/self.eta_ac for t1 in self.T_ac) +sum(self.X_tp[tb] for tb in self.T_ht if tb not in self.T_CHP.to_list()))
            else:
                self.X_ptw['CHP']=val
        self.heat_load_unmet=self.delta_h*self.eta_b[1]-sum(self.X_tp[t] for t in self.T_CHP) + self.X_ptw['CHP'] + sum(self.X_tp[t1]/self.eta_ac for t1 in self.T_ac)
        
        

            
            
#Dispatch Functions
    def monthly_demand(self):
        self.x_g=dict()
        self.z_ut=dict()
        for u in self.U_P:
            self.x_g[u]=np.zeros(len(self.H))
        elu=copy(self.elec_load_unmet)
        for m in self.M:
            u=1
            self.z_ut[m]=np.zeros(len(self.U_P))
            self.z_ut[m][u-1]=1
            limit=copy(self.deltabar_tu[u])
            bucket=0
            for h in self.Hm[m]:
                remain=limit-bucket
                load=elu[h]
                while load>0:
                    if load>remain:
                        self.x_g[u][h]+=remain
                        load-=remain
                        remain=0
                        u+=1
                        try:
                            self.z_ut[m][u-1]=1
                        except:
                            print(load,remain,u,self.z_ut)
                        limit=copy(self.deltabar_tu[u])
                        bucket=0
                        remain=limit-bucket
                    else:
                        self.x_g[u][h]+= load
                        bucket+=load
                        load=0

        

    def update_model(self,model):
        model.xbkwh=sum(self.X_bkWh[b] for b in self.B_e)
        model.elec_load_unmet =self.elec_load_unmet
        prod_remaining=dict()
        prod_remaining['CHP']=self.Z_to['CHP']*(self.X_sigma['CHP']-self.xrp['CHP'])
        for t in self.T_td:
            prod_remaining[t]=self.extra_renew[t]
        model.prod_remain = prod_remaining
        
        
    def partition_hours(self, H):
        vals=np.zeros((2,len(H)))
        charge_hours=[]
        for h in range(len(H)):
            if 0< self.remaining_prod['CHP'][h]<= self.X_sigma['CHP']:
                charge_hours+=h
            vals[1,h]=self.c_g[1][h+1]
        vals[1]=np.sort(vals[1])[::-1]
        return vals 
            
    def reduce_demand(self,h):
        m=self.Mh[h]
        r=self.Rh[h]
        target=min(self.targ_m[m],self.targ_r[r])
        load=self.elec_load_unmet[h]
        if load > target:
            dispatch=min(self.X_sigma['CHP'],max(load-target,self.td))
            self.run_chp(dispatch,h)
        return self.remaining_prod['CHP'][h]

    
    def turn_on_up_CHP(self,h,chp_remain,energ_cost):
        a=timer()
        t='CHP'
        elec_load=self.elec_load_unmet[h]
        heat_load=self.heat_load_unmet[h]
        chp_from_heat=self.chp_from_heat(heat_load,h)
        vals=[0,0]
        if self.Z_to[t][h]==0:
            vals[0]=max(self.td,min(chp_remain,chp_from_heat,elec_load))
        else:
            vals[0]=max(min(chp_remain,chp_from_heat,elec_load),0)
        vals[1]=max(min(chp_remain-vals[0],elec_load-vals[0]),0)
        b=timer()
        self.timer['energy_seg1']+=b-a
        for dispatch in vals:
            if self.Z_to[t][h]==0:
                dispatch=max(dispatch,self.td)
            a=timer()
            xtp=min(self.chp_xtp(dispatch,h),self.heat_load_unmet[h])
            b=timer()
            self.timer['energy_seg1']+=b-a
            
            a=timer()
            cost_chp=self.full_chp_cost(dispatch,h) #Second Green Box
            b=timer()
            self.timer['full_chp']+=b-a
            
            a=timer()
            cost_utility=energ_cost[h]*min(dispatch,self.elec_load_unmet[h]) + self.boiler_cost(xtp,h)#self.energy_cost(dispatch,h) + self.boiler_cost(xtp,h) #Second Green Box
            b=timer()
            self.timer['energy_boiler']+=b-a

            a=timer()
            
            if cost_chp <= cost_utility:                  #Fourth Diamond .... Cost without considering heat
                self.run_chp(dispatch,h)

            b=timer()
            self.timer['run_chp']+=b-a
                
    def dispatch_boiler(self,h,unmet_heat):
        t='BOILER'
        x_tp = unmet_heat/self.f_p[t][h]
        self.heat_load_unmet[h]=0
        self.X_tp[t][h]=x_tp #Boiler size could be determined during dispatch by computing optimal CHP dispatch and remaining is boiler size
        self.update_fuel(t,h)
        if x_tp>self.X_sigma['BOILER']: 
            self.X_sigma['BOILER']=x_tp

    
    def run_chp(self,kW,h):
        #TODO Can likely vectorize this and post process at the end
            t='CHP'
            x_rp=kW
           
            #Energy
            self.xrp[t][h]+=x_rp
            dispatch=min(self.elec_load_unmet[h],x_rp)
                        
            self.elec_load_unmet[h]-=dispatch
            self.extra_chp[h]+=x_rp-dispatch
            self.extra_elec[h]+=x_rp-dispatch
            self.remaining_prod[t][h]-=x_rp

            #Heat
            chp_x_tp=self.chp_xtp(x_rp,h)
            heat_unmet=self.heat_load_unmet[h]
                        
            if chp_x_tp > heat_unmet:
                self.X_tp[t][h]+= chp_x_tp
                self.heat_load_unmet[h] = 0
                self.Z_to[t][h]=1
                # print(h,chp_x_tp, heat_unmet)
                self.X_ptw[t][h] = chp_x_tp - heat_unmet
                self.curtailed_heat[h] = chp_x_tp - heat_unmet
            else:
                self.heat_load_unmet[h] -= chp_x_tp
                self.X_tp[t][h] += chp_x_tp
                self.Z_to[t][h]=1

            self.update_fuel(t,h)
    
    def storage_ops(self,batt_mod):
        if self.batt=='policy' and sum(self.X_bkWh[b] for b in self.B_e)>0:
            if not hasattr(self, 'battery_peak'):
                self.battery_peak=dict()
                for m in self.M:
                    self.battery_peak[m]=max(self.elec_load_unmet[self.Hm[m]])-self.percent*sum(self.X_bkW[b] for b in self.B_e)            
            self.policy_storage_ops()
            
        if self.batt=='alp':
            ampl=build_ampl_model_therm(self)
            start=timer()
            # print(self.X_sigma,self.X_bkWh)
            ampl.solve() 
            self.batt_timer+=[timer()-start]
            self=update_ampl_solution_therm(self,ampl)
   
    def policy_storage_ops(self):
        for b in self.B_e:
            for m in self.M:
                for h in self.Hm[m]:
            # m=self.Mh[h]
                    discharge=0
                    charge_cap=0
                #discharge
                    if (self.elec_load_unmet[h]>self.battery_peak[m]) & (self.xse[b][h]>0):
                        dispatch=min(sum(self.eta_minus[b]*self.X_bkW[b] for b in self.B_e),self.eta_minus[b]*(self.xse[b][h]-self.wubar_mcp[b]*self.X_bkWh[b]),self.elec_load_unmet[h])
                        self.elec_load_unmet[h]-=dispatch
                        self.xse[b][h:]-=dispatch/self.eta_minus[b]
                        self.dfs[b][h]+=dispatch
                        discharge=1
                        
                    if (self.peak_ind[h]==1) & (self.xse[b][h]>0)& (self.dfs[b][h]<self.X_bkW[b]):
                        dispatch=min(sum(self.eta_minus[b] for b in self.B_e)*(self.X_bkW[b]-self.dfs[b][h]),self.eta_minus[b]*(self.xse[b][h]-self.wubar_mcp[b]*self.X_bkWh[b]),self.elec_load_unmet[h])
                        self.elec_load_unmet[h]-=dispatch
                        self.xse[b][h:]-=dispatch/self.eta_minus[b]
                        self.dfs[b][h]+=dispatch
                        discharge=1
                        
                #charge
                    if (self.peak_ind[h]==0) & (self.elec_load_unmet[h]<=self.battery_peak[m]) & (discharge==0):
                        for t in self.T_e:
                            if self.xse[b][h]>=self.X_bkWh[b]: break
                            if charge_cap>=self.X_bkW[b]: break
                            if self.remaining_prod[t][h]>0:
                                charge=min(self.X_bkW[b]-charge_cap,self.X_bkWh[b]-self.xse[b][h],self.remaining_prod[t][h])
                                if (self.Z_to['CHP'][h] == 0) & (charge<self.td): break
                                charge_cap+=charge
                                if t == 'CHP':
                                    if self.Z_to['CHP'][h] == 0:
                                        dispatch=max(charge,self.td)
                                    else:
                                        dispatch=charge
                                    self.xrp[t][h]+=dispatch
                                    self.extra_elec+=dispatch-charge
                                    heat=self.chp_xtp(dispatch,h)
                                    self.X_tp[t][h]+=heat
                                    self.heat_load_unmet[h]=max(self.heat_load_unmet[h]-heat,0) 
                                else:
                                    self.xrp[t][h]+=charge
                                                        
                                self.pts[(b,t)][h]+=charge
                                self.remaining_prod[t][h]-=charge
                                self.xse[b][h:]+=charge*self.eta_plus[b,t]
                                self.Z_to[t][h]=1
                                self.update_fuel(t,h) #might need to exclude PV
                    
                        if (self.xse[b][h]<self.X_bkWh[b]) & (charge_cap<self.X_bkW[b]):
                            charge=min(self.X_bkWh[b]-self.xse[b][h],self.X_bkW[b]-charge_cap)
                            self.gts[h]+=charge
                            self.elec_load_unmet[h]+=charge
                            self.xse[b][h:]+=charge*self.eta_gplus[1]
   
    #Helper Functions
    
    def full_chp_cost(self,kW,h):
        t='CHP'
        # f=self.F_t[t]
        om_cost=self.chp_om*kW
        fuel=self.chp_fuel_intercept[h]*(1-self.Z_to[t][h]) + self.chp_fuel_slope[h]*kW
        fuel_cost=self.chp_fuel_cost*fuel
        chp_cost=om_cost+fuel_cost
        
        return chp_cost

    def chp_from_heat(self,heat,h):
        t='CHP'
        if self.Z_to[t][h]==1:
            x_tpb = 0
        else:
            x_tpb = self.k_tp[t]*self.X_sigma[t]
        kwh=((heat/(self.f_ha[t][h]*self.f_ht[t][h]) - x_tpb)/self.k_te[t])/self.f_l[t]
        return kwh
    

    def chp_xtp(self,x_rp,h):
        #TODO Need to pre-compute the slope and intercept########
        t='CHP'
        x_tp=self.f_ha[t][h]*self.f_ht[t][h]*((self.k_te[t])*self.f_p[t][h]*x_rp + self.k_tp[t]*self.X_sigma[t]*(1-self.Z_to[t][h]))
        return x_tp

    ###########
    ###########
    #Boiler
    ###########
    ###########


    def boiler_cost(self,kw,h):
        return self.boiler_cost_param*kw
    
    ###########
    ###########
    #Lifecycle Costs
    ###########
    ###########


    def cap_cost(self):
        cost=0
        for t in self.T:
            cost+=(self.X_sigma[t]*self.c_cm[t,1])
            cost+=(self.Z_sigmas[t,1,1]*self.c_cb[t,1])
        for b in self.B:
            cost+=(self.X_bkW[b]*self.c_kW[b])
            cost+=(self.X_bkWh[b]*(self.c_kWh[b]+self.c_omb[b]))
        return cost
    #TODO Update capital cost with the segments, should be able to move away from Z_sigmas
            
    def fix_om_cost(self):
        cost =0
        for t in self.T:
            cost += (1-self.f_tow)*self.f_om*(self.X_sigma[t]*self.c_omsigma[t])
        return cost

    def var_om_cost(self):
        cost =0
        for t in self.T:
            cost += (1-self.f_tow)*self.f_om*(self.c_omp[t] * self.xrp[t].sum())
        return cost    

    def fuel_cost(self):    
        cost=0
        for f in self.F:
            for t in self.Tf[f]:
                cost+=(1-self.f_tot)*(self.Delta[1]* self.c_u[f]*self.f_pf[t]* self.X_f[t].sum())
        return cost

    def energy_fixed_cost(self):
        cost = (1-self.f_tot)*self.f_e[1]*self.Delta[1]*(+self.X_mc[1] + self.c_afc[1])
        return cost

    def total_energy_cost(self):
        cost=(1-self.f_tot)*self.f_e[1]*self.Delta[1]*(sum(self.x_g[u][h]*self.c_g[u][h] for u in self.U_P for h in self.H))
        return cost
    #TODO confirm this is accurate

    def month_demand_cost(self):
        if self.c_rm is None or self.c_rm.sum()==0:
            return 0
        else:
            cost=(1-self.f_tot)*self.f_e[1]*(self.c_rm*self.X_dn).sum()
            return cost
        #TODO extend to multiple tiers in each month deltabar_mt

    def tou_demand_cost(self):
        if self.c_r is None or self.c_r.sum()==0:
            return 0
        else: 
            cost=(1-self.f_tot)*self.f_e[1]*(self.c_r*self.X_de).sum()
            return cost
        #TODO extend to multiple tiers in each month deltabar_t

    def exports_cost(self):
        cost=(self.c_e[1] * self.stg).sum()
        for t in self.T_e:
            cost+=(self.c_e[1]*self.ptg[t]).sum()
        return (1-self.f_tot)*self.f_e[1]*self.Delta[1]*cost    
    #
    def pi_cost(self):
        cost=(1-self.f_tow)* self.X_pi.sum()
        return cost

    def lcc(self):
        self.update_z_sigma()
        cost=self.cap_cost()+self.fix_om_cost() + self.var_om_cost()+self.fuel_cost()+self.total_energy_cost()+self.month_demand_cost() +self.tou_demand_cost() +self.energy_fixed_cost() -self.exports_cost() - self.pi_cost()     
        return cost
    
#####################
#####################
##Helper Functions###
#####################
#####################                        

    def display_solution(self):
        xrp=0
        for key,val in self.xrp.items():
            xrp+=np.sum(val)
        return ('case1',#"heuristic_"+str(self.X_sigma['PV1'])+"_"+str(self.X_sigma['CHP']),
                round(self.lcc(),0),
                round(self.delta_d.sum(),0) ,
                round(self.elec_load_unmet.sum(),0)  ,
                round(self.total_energy_cost(),0) ,
                round(self.tou_demand_cost(),0) ,
                round(self.month_demand_cost(),0),
                round(xrp,0),
                # self.Z_to.sum(),
                round(self.fuel_cost(),0),
                round(self.cap_cost(),0),
                round(self.fix_om_cost() + self.var_om_cost(),0))

    def print_time(self,start):
        stop=timer()
        print(round(stop-start,4))
        return stop


    def update_z_sigma(self):
        for t in self.T:
            if self.X_sigma[t]>0:
                self.Z_sigmas[t,1,1]=1
            else:
                self.Z_sigmas[t,1,1]=0
    
    def fuel_conversion(self,amt,t,h):
        if t=='CHP':
            fuel = self.f_fa[t][h]*(self.m_fbm[t]*self.X_sigma[t] + self.f_p[t][h]*self.m_fm[t]*amt)
        else:
            fuel = self.m_fm[t]*amt
        return fuel

    def update_fuel(self,t,h):
        if t == 'CHP':
            self.X_f[t][h]=self.fuel_conversion(self.xrp[t][h],t,h)
        else:
            self.X_f[t][h]=self.fuel_conversion(self.X_tp[t][h],t,h)

                    
    def demand_cost(self,amt,tiers):
        m,r=tiers
        for n in self.N:
            if amt <= self.deltabar_mt[n]:
                month = (1-self.f_tot)*self.f_e[1]*self.Delta[1]*self.c_rm[m,n]*amt
                break
        for e in self.E:
            if amt <= self.deltabar_t[e]:
                tou = (1-self.f_tot)*self.f_e[1]*self.Delta[1]*self.c_r[r,e]*amt        
                break
        return month+tou
    ##############################################################################
    ######Helper function to populate X_dn and X_de###############################
    
    def fill_peaks(self):
        for m in self.M:
            grid=self.elec_load_unmet[self.Hm[m]]
            if self.T_rdc.empty:
                max_val=max(grid)
            else:
                prod=sum(self.xrp[t][self.Hm[m]] for t in self.T_rdc)
                max_val=max(grid+prod)
            for n in self.N:
                # print(m,n,max_val,self.deltabar_mt[n])
                if max_val>self.deltabar_mt[n]:
                    self.X_dn[m,n]=self.deltabar_mt[n]
                    max_val-=self.deltabar_mt[n]
                else:
                    self.X_dn[m,n]=max_val
                    break
        for r in self.R:
            if len(self.Hr[r])>=1:
                grid=self.elec_load_unmet[self.Hr[r]]
                if self.T_rdc.empty:
                    max_val=max(grid)
                else:
                    prod=sum(self.xrp[t][self.Hr[r]] for t in self.T_rdc)
                    max_val=max(grid+prod)
                self.X_de[r,1]=max_val
    
    
    ##################################################################################
    #####Sets demand targets as a function of CHP capacity and storage capacity#######
    ##################################################################################    
    def set_demand_targets(self):
        s=timer()
        capacity=self.X_sigma['CHP']
        self.targ_m = dict()
        self.targ_r = dict()
        for m in self.M: 
            if len(self.Hm[m])>0:
                self.targ_m[m]=max(self.elec_load_unmet[self.Hm[m]])
            else:
                self.targ_m[m]=max(self.elec_load_unmet)
        for r in self.R: 
            if len(self.Hr[r])>0:
                self.targ_r[r]=max(self.elec_load_unmet[self.Hr[r]])
            else:
                self.targ_r[r]=max(self.elec_load_unmet)
        best=1000000000000
        # print(f"Part 1:{timer()-s}")
        s=timer()
        #############################################################################################################
        ####Next loop is intended to choose the right demand targets based on how variable the load profile is#######
        ####Low variance likely leads to less peak shaving and high variance leads to high peak shaving##############
        #############################################################################################################
        #TODO Fix this part, the cost structure in create_segments is off. Use the new parameters to compute costs
        for perc in np.arange(self.fubar_td['CHP'],1.1,0.1):        
            # vals=self.create_segments(perc*capacity)
            cost,target_m,target_r=get_targets(self,perc*capacity)
            if cost<best:
                self.targ_m=target_m
                self.targ_r=target_r
                best=cost
        
        self.battery_peak=dict()
        for m in self.M:
            self.battery_peak[m]=self.targ_m[m]-self.percent*sum(self.X_bkW[b] for b in self.B_e)


    def check_NM(self):
        cap = sum(self.X_sigma[t] for v in self.V for t in self.Tv[v])
        self.Z_nmil*=0
        self.Z_nmil[self.V[1]]=1
        for v in self.V:
            if cap<=self.i_n[v]:
                for t in self.T_e:
                # for t in self.Tv[v]:
                    if self.X_sigma[t]>0:
                        val=self.X_sigma[t]
                        for tv in self.Tv[v]:
                            if tv[0:2] == t[0:2]:
                                self.X_sigma[t]=0
                                self.X_sigma[tv]=val
                self.Z_nmil*=0
                self.Z_nmil[v]=1
                return