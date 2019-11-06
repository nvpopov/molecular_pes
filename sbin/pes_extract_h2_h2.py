#!/usr/bin/python3
import os
import sys
from operator import itemgetter, attrgetter, methodcaller
from math import degrees

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

arg_batchdir = ""
arg_outfile = ""
arg_task_per_file = 21

if len(sys.argv) > 2:
    arg_batchdir = sys.argv[1]
    arg_outfile = sys.argv[2]
    arg_task_per_file = int(sys.argv[3])

    print("batchdir = {}, outfile = {}, task_per_file = {}".format(arg_batchdir, arg_outfile, arg_task_per_file))
else:
    print("pes_extract.py <batchdir> <outfile>")
    exit()

task_list = os.listdir(arg_batchdir)

task_data_accum = []
task_data_energy = []

total_bad_tasks = 0
for task_id, task in enumerate(task_list):
    print("Processing {} of {} ...".format(task_id, len(task_list)))
    task = "{}/{}".format(arg_batchdir, task)

    for orca_inp in open("{}/orca.inp".format(task)):
        if '#' in orca_inp:
            task_data = list()
            splt = orca_inp.split()
            #print(splt)
            #R, TT, T, P, do_efield, Q, bsse
            r = float(splt[1])
            tt = float(splt[2])
            t = float(splt[3])
            p = float(splt[4])
            field = str2bool(splt[5])
            q = float(splt[6])
            bsse = splt[7]
          #
          #  task_data["r"] = r
          #  task_data["tt"] = tt
          #  task_data["t"] = t
          #  task_data["p"] = p
          # task_data["field"] = field
          #  task_data["q"] = q
          #  task_data["bsse"] = bsse
          #  task_data["t_normally"] = False
            # 0_r 1_rr 2_t 3_p 4_field 5_q 6_bsse 7_tnorm 8_en 9dhx 10dhy 11dhz 12dcx 13dcy 14dcz
            task_data=[r, tt, t, p, field, q, bsse, False, -100.0, 0, 0, 0, 0, 0, 0]
            task_data_accum.append(task_data)

    oi_cnt = 0
    dipole_counter = 0
    for orca_out in open("{}/orca.out".format(task)):


        if "FINAL SINGLE POINT ENERGY" in orca_out:
            splt = orca_out.split()
            task_data_energy.append(float(splt[4]))
            #print(task_data_accum[task_id*36 + oi_cnt][8],task_id + oi_cnt )
            task_data_accum[task_id*arg_task_per_file + oi_cnt][8] = float(splt[4])
            task_data_accum[task_id*arg_task_per_file + oi_cnt][7] = True
            oi_cnt += 1

        if "Total Dipole Moment" in orca_out:
            #print("#######")
            #print(orca_out.split())
            #print(task)
           # print(dipole_counter)
            #print(task_data_accum[task_id * 36 + oi_cnt-1], oi_cnt)
           # print("#######")

            dipole_moment_vector = [float(i) for i in orca_out.split()[4:7]]
            if dipole_counter == 0:
                task_data_accum[task_id * arg_task_per_file + oi_cnt-1][9]  = dipole_moment_vector[0]
                task_data_accum[task_id * arg_task_per_file + oi_cnt-1][10] = dipole_moment_vector[1]
                task_data_accum[task_id * arg_task_per_file + oi_cnt-1][11] = dipole_moment_vector[2]
                dipole_counter+=1
            elif dipole_counter == 1:
                task_data_accum[task_id * arg_task_per_file + oi_cnt-1][12] = dipole_moment_vector[0]
                task_data_accum[task_id * arg_task_per_file + oi_cnt-1][13] = dipole_moment_vector[1]
                task_data_accum[task_id * arg_task_per_file + oi_cnt-1][14] = dipole_moment_vector[2]
                dipole_counter = 0

print(len(task_data_accum))
print(len(task_data_energy))
task_data_accum_sol = []

for i in range(0, int(len(task_data_accum)/3)):
    t1 = task_data_accum[i * 3]
    t2 = task_data_accum[i * 3 + 1]
    t3 = task_data_accum[i * 3 + 2]
    en = (t1[8] - t2[8] - t3[8])*220000
    t1[8] = en
    task_data_accum_sol.append(t1)

tdas_sorted = sorted(task_data_accum_sol, key=itemgetter(0,1,2))

out_pes = open(arg_outfile, "w")

out_pes.write("r(A),tt(D),t(D),p(D),ec(cminv),dx(au),dy(au),dz(au)\n")
for task_ in tdas_sorted:
    out_pes.write("{0:.5f},{1:.5f},{2:.5f},{3:.5f},{4:.5f},{5:.5f},{6:.5f},{7:.5f}\n".format(
        task_[0], task_[1], task_[2], task_[3], task_[8], task_[12], task_[13], task_[14]))




