# File: optimization/ampl_output.py
# Purpose: Format and write decision variables from system model to AMPL .dat file format

import numpy as np

def outputALL(d, file_name):
    sigma_out(d.X_sigma, file_name, "x_sigma")
    x_sigmas(d.X_sigma, file_name, "x_sigmas")
    z_sigmas(d.X_sigma, file_name, "z_sigmas")
    num_h_out(d.dfs, file_name, 'x_dfs')
    num_h_out(d.x_g, file_name, 'x_g')
    zut_out(d.z_ut, file_name, 'z_ut')
    num_h_out(d.xse, file_name, 'x_se')
    b_out(d.X_bkW, file_name, "x_bkw")
    b_out(d.X_bkWh, file_name, "x_bkwh")
    th_out(d.xrp, file_name, "x_rp")
    th_out(d.X_tp, file_name, "x_tp")
    th_out(d.Z_to, file_name, "z_to")
    th_out(d.X_f, file_name, "x_f")
    pts_out(d, file_name, "x_pts")
    th_out(d.ptg, file_name, "x_ptg")
    month_out(d.X_dn, file_name, "x_dn")
    tou_out(d.X_de, file_name, "x_de")
    intercept(d, file_name, "x_fb")
    gts(d.gts, file_name, "x_gts")
    gts(d.stg, file_name, "x_stg")
    ptw(d.X_ptw, file_name, "x_ptw")
    z_nmil(d.Z_nmil, file_name, "z_nmil")

def sigma_out(var, name, var_name):
    with open(name, "w") as f:
        f.write(f'param {var_name} := \n')
        for t, val in var.items():
            if val == 0: continue
            f.write(f"{str(t[0])} {val}\n")
        f.write(";\n\n")

def x_sigmas(var, name, var_name):
    pass  # placeholder if needed

def z_sigmas(var, name, var_name):
    pass  # placeholder if needed

def num_h_out(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for key, val in var.items():
            for i in range(len(val)):
                if val[i] == 0: continue
                f.write(f"{key} {i+1} {val[i]}\n")
        f.write(";\n\n")

def zut_out(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for key, val in var.items():
            for u in range(len(val)):
                if val[u] == 0: continue
                f.write(f"{int(key)} {u+1} {val[u]}\n")
        f.write(";\n\n")

def b_out(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for b, val in var.items():
            if val == 0: continue
            f.write(f"{b[0]} {val}\n")
        f.write(";\n\n")

def th_out(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for t, arr in var.items():
            for h in range(8760):
                if arr[h] < 1e-9: continue
                if var_name == 'x_pts':
                    f.write(f"3 {t} {h+1} {arr[h]}\n")
                elif var_name == 'x_ptg':
                    f.write(f"{t} 1 {h+1} {arr[h]}\n")
                else:
                    f.write(f"{t} {h+1} {arr[h]}\n")
        f.write(";\n\n")

def pts_out(d, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for (b, t), val in d.pts.items():
            for h in range(8760):
                if val[h] < 1e-9: continue
                f.write(f"{b} {t} {h+1} {val[h]} \n")
        f.write(";\n\n")

def month_out(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for (m, n), val in var.items():
            if val == 0: continue
            f.write(f"{m} {n} {val}\n")
        f.write(";\n\n")
        f.write("param z_dmt := \n")
        for (m, n), val in var.items():
            if val == 0: continue
            f.write(f"{m} {n} 1\n")
        f.write(";\n\n")

def tou_out(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for tier in var.index.unique(0).to_list():
            if var[tier, 1] == 0: continue
            f.write(f"{tier} 1 {var[tier,1]}\n")
        f.write(";\n\n")
        f.write("param z_dt := \n")
        for tier in var.index.unique(0).to_list():
            if var[tier, 1] > 0:
                f.write(f"{tier} 1 1\n")
        f.write(";\n\n")

def intercept(d, name, var_name):
    pass  # implement as needed

def gts(arr, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for h in range(len(arr)):
            if arr[h] > 0:
                f.write(f"{h+1} {arr[h]}\n")
        f.write(";\n\n")

def ptw(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for h in range(len(var['CHP'])):
            if var['CHP'][h] > 0:
                f.write(f"CHP {h+1} {var['CHP'][h]}\n")
        f.write(";\n\n")

def z_nmil(var, name, var_name):
    with open(name, "a") as f:
        f.write(f'param {var_name} := \n')
        for k, v in var.items():
            f.write(f"{k} {v}\n")
        f.write(";\n\n")
