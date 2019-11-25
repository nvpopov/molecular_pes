#!/usr/bin/python3

from math import *
import os
import sys
import aux_data
from operator import itemgetter
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

def pes_extract_a(batch_dir_inp, output_filename):
    
    task_list = os.listdir(batch_dir_inp)
    print("Total steps = {}".format(task_list.__len__()))

    out_data = []
    r_domain = []
    theta_domain = []

    for task_id, task in enumerate(task_list):
        #print ("{} {}".format(task_id, task))
        task_data = {}

        batch_dir = "{}/{}".format(batch_dir_inp, task)
        #open DATA file
        data_file = open("{}/DATA".format(batch_dir), "r")
        data_file_lines = data_file.readlines()[0].split()
        task_r = float(data_file_lines[0])
        task_theta = float(data_file_lines[1])

        outfile = open("{}/orca.out".format(batch_dir), "r")

        toten_hf = []
        toten_ccsdt = []
        cbs_scf = []
        cbs_mdci = []
        cbs_final = []
        dipople_moment = []

        for line in outfile.readlines():
        
            if "Total Energy       :" in line:
                toten_hf.append(float(line.split()[3]))    
                pass
        
            if "E(CCSD(T))" in line:
                toten_ccsdt.append(float(line.split()[2]))    
                pass

            if "Extrapolated CBS SCF energy (4/5)" in line:
                cbs_scf.append(float(line.split()[6]))    
                pass

            if "Extrapolated CBS correlation energy (4/5) :" in line:
                cbs_mdci.append(float(line.split()[6]))    
                pass

            if "Estimated CBS total energy (4/5) : " in line:
                cbs_final.append(float(line.split()[6]))    
                pass

            if "Total Dipole Moment" in line:
                lsp = line.split()
                dipople_moment.append([float(lsp[4]), float(lsp[5]), float(lsp[6])])

        if task_r not in r_domain:
            r_domain.append(task_r)        

        if task_theta not in theta_domain:
            theta_domain.append(task_theta)

        task_data["r"] = task_r
        task_data["theta"] = task_theta
        task_data["E_HF_AB"] = toten_hf[0]
        task_data["E_HF_A"] = toten_hf[1]
        task_data["E_HF_B"] = toten_hf[2]   
        task_data["E_HF_BSSE"] = toten_hf[0] - toten_hf[1] - toten_hf[2]
        task_data["E_CCSDT_AB"] = toten_ccsdt[0]         
        task_data["E_CCSDT_A"] = toten_ccsdt[1]
        task_data["E_CCSDT_B"] = toten_ccsdt[2]
        task_data["E_CCSDT_B"] = toten_ccsdt[0] - toten_ccsdt[1] - toten_ccsdt[2]
        task_data["E_MDCI_CBS_AB"] = cbs_final[0]
        task_data["E_MDCI_CBS_A"] = cbs_final[1]
        task_data["E_MDCI_CBS_B"] = cbs_final[2]
        task_data["E_FINAL"] = aux_data.ha2cminv*(cbs_final[0] - cbs_final[1] - cbs_final[2])
        task_data["DIPOLE_MOMENT_MDCI"] = dipople_moment[1]

        out_data.append(task_data)

    out_data_sorted = sorted(out_data, key = lambda i: (i['theta'], i['r'])) 
    r_domain_s = sorted(r_domain)
    theta_domain_s = sorted(theta_domain)

    # print(r_domain_s, theta_domain_s)
    # z = []
    # for elem in out_data_sorted:
    #     z.append(elem["E_FINAL"])
    # print(z)

    # print("Min Z = {}".format(min(z)))

    # x = np.array(r_domain_s)
    # y = np.array(theta_domain_s)
    # X,Y = np.meshgrid(x,y)
    # Z = np.array(z).reshape(len(theta_domain_s), len(r_domain_s))
    # h = plt.contourf(X,Y,Z, 1000, levels = [-100 + a*2 for a in range(0, 1000)])
    # plt.show()

    output_data = open(output_filename, "w")
    output_data.write("r_angs,theta_rad,e_final_cminv,dm_x_au,dm_y_au,dm_z_au\n")
    for elem in out_data_sorted:
        output_data.write("{},{},{},{},{},{}\n".format(elem["r"], elem["theta"], elem["E_FINAL"],
        elem["DIPOLE_MOMENT_MDCI"][0], elem["DIPOLE_MOMENT_MDCI"][1], elem["DIPOLE_MOMENT_MDCI"][2]))

    #generate data for plots
    for t_i in range(0, len(theta_domain_s)):

        cur_r = open("r{}.data".format(t_i), "w")
        
        data_per_theta = [elem for elem in out_data_sorted if elem["theta"] == theta_domain_s[t_i]]

        for elem2 in data_per_theta:
            cur_r.write("{},{}\n".format(elem2["r"], elem2["E_FINAL"]))
   

    pass