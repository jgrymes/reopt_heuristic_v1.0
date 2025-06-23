param horizon_length default 8760;
set H := 1..horizon_length;
set M := 1..12;
set R;
set Hm{M};
set Hr{R};
set T;
set T_td;
set T_CHP;
set U_P;

param xbkwh; 
param xbkw :=xbkwh/2;
param elec_load_unmet{H};
param prod_remain{T,H};
param kappa_e{H};
param kappa_d{M};
param kappa_r{R};
param kappa_c{T_CHP,H};
param eta_plus{T};
param eta_gplus;
param eta_minus;
param w_zero;
param wubar_mcp;

var xsc{0..card(H)} >=wubar_mcp*xbkwh <= xbkwh;
var dfs{H} >=0, <=xbkw;
var gts{H} >= 0, <=xbkw;
var xg{H} >= 0;
var xd{M} >= 0;
var xr{R} >= 0;
var pts{t in T, h in H} >= 0 , <= prod_remain[t,h]; 


minimize lcc: sum{h in H}kappa_e[h]*xg[h] + sum{m in M}kappa_d[m]*xd[m]+ sum{r in R}kappa_r[r]*xr[r] + sum{t in T_CHP,h in H}kappa_c[t,h]*pts[t,h];

    # energy_cost=(1-dispatch.f_tot)*dispatch.f_e[1]*sum(dispatch.c_g[1][h]*model.xg[h] for h in model.H) 
    # demand_cost=(1-dispatch.f_tot)*dispatch.f_e[1]*sum(dispatch.c_rm[m][1]*model.xd[m] for m in dispatch.M) 
    # battery_cost=(dispatch.c_kW[3]*model.xbkw) 
    # chp_cost=sum(dispatch.chp_fuel_slope[h]*dispatch.chp_fuel_cost*model.pts['CHP',h] for h in model.H)

s.t. state_of_charge{h in H}: xsc[h] = xsc[h-1] + sum{t in T} eta_plus[t]*pts[t,h] + eta_gplus*gts[h] - dfs[h]/eta_minus;
s.t. state_of_charge_0: xsc[0] = w_zero*xbkwh;
s.t. peak_load{m in M, h in Hm[m]}: xd[m]>= xg[h];
s.t. peak_ratchet{r in R, h in Hr[r]}: xr[r]>= xg[h];
s.t. loadBalance{h in H}: xg[h] - gts[h] + dfs[h] >= elec_load_unmet[h]+ sum{t in T_td}pts[t,h];
s.t. max_gts{h in H}: xg[h]>=gts[h];
# s.t. max_gts{h in H}: xg[h]>=gts[h];
s.t. max_power_in{h in H}: sum{t in T} pts[t,h] +gts[h] +dfs[h] <= xbkw;