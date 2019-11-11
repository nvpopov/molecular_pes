#!/usr/bin/python3

from math import *
import os
import sys
from operator import itemgetter

batchdir = "batch"

task_list = os.listdir(batchdir)

output_batch = []

for task_id, task in enumerate(task_list):
    #print ("{} {}".format(task_id, task))

    batch_dir = "batch/{}".format(task)
    #open DATA file
    data_file = open("{}/DATA".format(batch_dir), "r")

    total_batch_in_file = 0
    batch_data = []
    for data_file_line in data_file.readlines():
        total_batch_in_file += 1
        data_file_ls = data_file_line.split()
        batch_data_line = []
        for df_ls in data_file_ls:
            batch_data_line.append(float(df_ls))
        batch_data.append(batch_data_line)
    
    #print(batch_data)
    
    outfile = open("{}/orca.out".format(batch_dir), "r")
    
    reads_hit = 0
    batch_id = 0

    energies_bsse = []
    dipole_moments = []

    for out_line in outfile.readlines():
        
        if "Total Energy       :" in out_line:    
            out_line_split = out_line.split()
            en = float(out_line_split[3])
            reads_hit+=1
            energies_bsse.append(en)

        #if "Total Dipole Moment    :" in out_line:
        if reads_hit > 2:
            data_for_current_batch = [elem for elem in batch_data[batch_id]]

            for i in range(0,3):
                data_for_current_batch.append(energies_bsse[i])

            output_batch.append(data_for_current_batch)
            reads_hit = 0
            batch_id+=1
            energies_bsse = []
            dipole_moments = []

sorted_output_batch = sorted(output_batch, key=itemgetter(0,1))

print("r,phi,e_bsse(cminv),e(AB,Hartree),e(ABg,Hartree),e(AgB,Hartree)")

for srt_b in sorted_output_batch:
    bsse_corrected_en = (srt_b[2] - srt_b[3] - srt_b[4]) * 220000.00000000003
    #220000.00000000003
    print("{},{},{},{},{},{}".format(srt_b[0], srt_b[1], bsse_corrected_en, srt_b[2], srt_b[3], srt_b[4]))