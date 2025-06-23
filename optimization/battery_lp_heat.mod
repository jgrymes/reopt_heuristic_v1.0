param horizon_length default 8760;
set H;# := 1..horizon_length;
set M;# := 1..12;
set R;
set Hm{M};
set Hr{R};
set B;
set B_e;
set B_c; #default {2};
set B_h;# default {1};
set T; #Renews,CHP,AC,EC,BOIL
set T_e;
set T_h;
set T_c;
set T_CHP;
set T_th;
set T_ec;
set T_ac;
set Tb{B};
#Neeed
set T_s;
set T_td;


param xbkwh{B}; 
param w_zero{B};
param wubar_mcp{B};
param w_d{B};
param eta_minus{B};
param xbkw{B};# :=xbkwh[b]/2;
param load_e{H} default 0;
param load_h{H} default 0;
param load_c{H} default 0;
param cap_e{T_e,H} default 0;
param cap_th{T_th,H} default 0;
param cap_free{H} default 0;
param kappa_e{H} default 0;
param kappa_d{M} default 0;
param kappa_r{R} default 0;
param kappa_c{T_h,H} default 0;
param eta_plus{b in B,Tb[b]};
param eta_gplus;
param eta_ac;
param eta_ec;
param chp_heat_slope{H} default 0;
param min_td{H} default 0;
param chp_heat_intercept{H} default 0;

#need to add
param deltabar_gs default 0;
param kappa_s{H} default 0;
param sales_cap{T_s,H} default 0;
param renew{T_td,H} default 0;

#Storage
var xsc{b in B,0..card(H)} >=wubar_mcp[b]*xbkwh[b] <= xbkwh[b];
var dfs{b in B,H} >=0, <=xbkw[b];
var pts_e{t in T_e, h in H} >= 0 , <= cap_e[t,h]; 
var pts_th{t in T_th, h in H} >= 0 , <= cap_th[t,h]; 
var gts{H} >= 0, <=sum{b in B_e} xbkw[b];
var stg{H} >=0, <= sum{b in B_e} xbkw[b];

#Utility
var xg{H} >= 0;
var xd{M} >= 0;
var xr{R} >= 0;

#Production
var xrp{t in T_CHP, h in H} >= min_td[h], <= cap_e[t,h]; 
var xtp{t in T_th, h in H} >= 0 , <= cap_th[t,h]; 
var ptg{t in T_s, h in H} >= 0 , <= sales_cap[t,h];

var ab_chl >=0, <= max{h in H} load_c[h];

minimize lcc:  sum{h in H}kappa_e[h]*xg[h] + sum{m in M}kappa_d[m]*xd[m]+ sum{r in R}kappa_r[r]*xr[r] + sum{t in T_CHP,h in H}kappa_c[t,h]*xrp[t,h]+ sum{t in T_h,h in H:t not in T_CHP}kappa_c[t,h]*xtp[t,h]-sum{h in H}(kappa_s[h]*(sum{t in T_s}ptg[t,h] + stg[h])) + 569.142173*ab_chl;

#Load Balance Constraints
s.t. loadBalanceElec{h in H}: sum{t in T_CHP} xrp[t,h] + sum{t in T_td} renew[t,h] + xg[h] + sum{b in B_e} dfs[b,h] >= load_e[h]+ sum{t in T_e}pts_e[t,h] + sum{t in T_s}ptg[t,h] + stg[h] + gts[h] + sum{t in T_ec}xtp[t,h]/eta_ec;
s.t. loadBalanceHeat{h in H}: sum{t in T_h} xtp[t,h] + sum{b in B_h} dfs[b,h] >= load_h[h]+ sum{t in T_h}pts_th[t,h]+ sum{t in T_ac}xtp[t,h]/eta_ac; #+ xtp_f[h]
s.t. loadBalanceCool{h in H}: sum{t in T_c} xtp[t,h] + sum{b in B_c} dfs[b,h] >= load_c[h] + sum{t in T_c}pts_th[t,h];

#Battery State of Charge
s.t. state_of_chargeElec{b in B_e,h in H}: xsc[b,h] = xsc[b,h-1] + sum{t in T_e} eta_plus[b,t]*pts_e[t,h] + eta_gplus*gts[h] - dfs[b,h]/eta_minus[b];
s.t. state_of_chargeThem{b in B,h in H: b not in B_e}: xsc[b,h] = xsc[b,h-1]+ sum{t in Tb[b]} eta_plus[b,t]*pts_th[t,h] - dfs[b,h]/eta_minus[b]-w_d[b]*xsc[b,h];
s.t. state_of_charge_0 {b in B}: xsc[b,0] = w_zero[b]*xbkwh[b];

#Utility
s.t. peak_load{m in M, h in Hm[m]}: xd[m]>= xg[h];
s.t. peak_ratchet{r in R, h in Hr[r]}: xr[r]>= xg[h];
s.t. grid_to_store{h in H}: xg[h]>=gts[h];
s.t. store_to_grid{h in H}: sum{b in B_e} dfs[b,h] >= stg[h];
s.t. max_sales1: sum{h in H} (sum{t in T_s} ptg[t,h] + stg[h]) <= deltabar_gs;
s.t. grid_arbitrage: sum{h in H} (sum{t in T_s} ptg[t,h] + stg[h]) <= sum{h in H} xg[h];

#Microgrid Prodution
s.t. max_ptsElec{t in T_CHP, h in H}: xrp[t,h]>=pts_e[t,h];# + ptg[t,h]
s.t. max_ptsElec_a{t in T_CHP, h in H: t in T_s}: xrp[t,h]>=pts_e[t,h] + ptg[t,h];
s.t. max_ptsTherm{t in T_th, h in H:t not in T_CHP}: xtp[t,h]>=pts_th[t,h];


#Storage Boundry 
s.t. max_power_inElec{b in B_e,h in H}: sum{t in Tb[b]} pts_e[t,h] +gts[h] +dfs[b,h] <= xbkw[b];
s.t. max_power_inTherm{b in B,h in H: b not in B_e}: sum{t in Tb[b]} pts_th[t,h] +dfs[b,h] <= xbkw[b];

#CHP
s.t. chp_thermal{t in T_CHP,h in H}: xtp[t,h] = chp_heat_slope[h]*xrp[t,h]+chp_heat_intercept[h];
s.t. chp_storage{t in T_CHP,h in H}: xtp[t,h] >= pts_th[t,h];#+ xtp_f[h] 

#ab_chl
s.t. ab_chl_siz{t in T_c, h in H: t <> 'ELECCHL'}: ab_chl >= xtp[t,h];








# s.t. max_bat {b in B}: xbkw[b] <=xbkwh[b];
# s.t. max_SoC {b in B,h in 0..card(H)}: xsc[b,h]<=xbkwh[b];
# s.t. min_SoC {b in B,h in 0..card(H)}: xsc[b,h]>=wubar_mcp[b]*xbkwh[b];

# s.t. discharge{h in H}: sum{t in T_e} pts_e[t,h] + gts[h] <=xbkw[3]* z[h];

# s.t. charge{h in H}: sum{b in B_e} dfs[b,h] <=xbkw[3] * (1-z[h]);

