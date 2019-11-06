from math import *
import os
import h2_aux
import sys

do_efield = False
Q = 0.47

if len(sys.argv) > 1:
    if 'doef' in sys.argv[1]:
        do_efield = True
        Q = float(sys.argv[2])

######
##
##  parameters
##

r = 0.7369617067 /2e0  ## H2 half bond length

R0 = 1.7
Rend=10.80
dR = .25
NR = int((Rend - R0)/dR)


T0 = 0
Tend=90
dT = 15
NT = int( (Tend-T0) / dT )

TT0 = 0
TTend=90
dTT = 15
NTT = int( (TTend-TT0) / dTT )

P0 = 0
Pend=180
dP = 30
NP = int( (Pend - P0)/dP)


Hmol1 = [[ r, 0, 0],
      [  -r, 0, 0]]


# print(NT)
# print(NP)
# print(NTT)
# print(NR)
# print(NT * NR * NTT * NP)
# exit(0)

nchgx=8
nchgy=8
dx=dy=5

offsetx = (nchgx - 1)*dx /3
offsety = (nchgy - 1)*dy /2
chg=[]

charge_file_data = "{}\n".format(8*8*2)

for i in range(nchgx):
    for j in range(nchgy):
        chg.append([ -Q, i*dx - offsetx,-10, j*dy - offsety ])
        chg.append([  Q, i*dx - offsetx,  25, j*dy - offsety ])
        charge_file_data += "{} {} {} {}\n".format(-Q, i*dx - offsetx,-10, j*dy - offsety )
        charge_file_data += "{} {} {} {}\n".format( Q, i*dx - offsetx,  25, j*dy - offsety )

counter = 0

batch_dir = "./batch"
if not os.path.exists(batch_dir):
    os.makedirs(batch_dir)

for IR in range(NR):
    R =  dR * IR + R0
    for ITT in range(NTT+1):

        TT = pi/180e0*( dTT*ITT + TT0 )
        O = [ R * cos(TT), R * sin(TT), 0]

        
        """
        print(4 + len(chg))
        print(R,T0,TT0,P0)
        print("H  ", *Hmol1[0]      )
        print("H  ", *Hmol1[1]      )
        print("H  ", *[  r + O[0] , O[1], O[2]])
        print("H  ", *[- r + O[0] , O[1], O[2]])
        [print("X ", *c) for c in chg]
        """
        for IT in range(NT+1):
            P_enter_again = False
            orca_gen_name = "r_{0}_tt_{1}_t_{2}".format(IR, ITT, IT)
            orca_gen_filename = "{0}/{1}/orca.inp".format(batch_dir, orca_gen_name)
            batch_exact_dir_name = "{0}/{1}".format(batch_dir, orca_gen_name)


            if not os.path.exists(batch_exact_dir_name):
                os.makedirs(batch_exact_dir_name)

            if do_efield:
                point_charge_filename = "{0}/{1}/ntc_charges.chg".format(batch_dir, orca_gen_name)
                file_charge = open(point_charge_filename, 'w')
                file_charge.write(charge_file_data)

            file_orca = open(orca_gen_filename, "w")

            for IP in range(NP+1):
                P = pi / 180e0 * (dP * IP + P0)
                T = pi / 180e0 * (dT * IT + T0)

                counter += 1

                if (P_enter_again):
                    file_orca.write("$new_job\n")

                # write dimer
                file_orca.write("# {0} {1} {2} {3} {4} {5} d\n".format(
                    R, degrees(TT), degrees(T), degrees(P), do_efield, Q
                ))
                file_orca.write("! RHF aug-cc-pVQZ RIJK AutoAux  VeryTightSCF CCSD(T)-F12/RI  DIRECT PMODEL")
                file_orca.write(h2_aux.h2cabs)
                if do_efield:
                    file_orca.write("% pointcharges \"ntc_charges.chg\"\n")

                file_orca.write(h2_aux.gen_orca_aux(3, 2800))

                Hmol2 = [[  r * cos(T) + O[0],
                            r *sin(T) * cos(P) + O[1],
                            r * sin(T) * sin(P)+ O[2]],
                         [ -r * cos(T) + O[0],
                           -r *sin(T) * cos(P) + O[1],
                           -r * sin(T) * sin(P)+ O[2]],
                         ]
                for mol_c in range(0,2):
                    file_orca.write("H {0:10.4f} {1:10.4f} {2:10.4f}\n".format(
                        Hmol1[mol_c][0], Hmol1[mol_c][1], Hmol1[mol_c][2]))

                for mol_c in range(0,2):
                    file_orca.write("H {0:10.4f} {1:10.4f} {2:10.4f}\n".format(
                        Hmol2[mol_c][0], Hmol2[mol_c][1], Hmol2[mol_c][2]))



                file_orca.write("*\n\n")

                #write monomer 1
                file_orca.write("# {0} {1} {2} {3} {4} {5} m1\n".format(
                    R, degrees(TT), degrees(T), degrees(P), do_efield, Q
                ))
                file_orca.write("$new_job\n")
                file_orca.write("! RHF aug-cc-pVQZ RIJK AutoAux  VeryTightSCF CCSD(T)-F12/RI  DIRECT PMODEL")
                file_orca.write(h2_aux.h2cabs)
                if do_efield:
                    file_orca.write("% pointcharges \"ntc_charges.chg\"\n")

                file_orca.write(h2_aux.gen_orca_aux(1, 2800))
                for mol_c in range(0, 2):
                    file_orca.write("H {0:10.4f} {1:10.4f} {2:10.4f}\n".format(
                        Hmol1[mol_c][0], Hmol1[mol_c][1], Hmol1[mol_c][2]))

                for mol_c in range(0, 2):
                    file_orca.write("H : {0:10.4f} {1:10.4f} {2:10.4f}\n".format(
                        Hmol2[mol_c][0], Hmol2[mol_c][1], Hmol2[mol_c][2]))


                file_orca.write("*\n\n")

                #write monomer 2
                file_orca.write("# {0} {1} {2} {3} {4} {5} m2\n".format(
                    R, degrees(TT), degrees(T), degrees(P), do_efield, Q
                ))
                file_orca.write("$new_job\n")
                file_orca.write("! RHF aug-cc-pVQZ RIJK AutoAux  VeryTightSCF CCSD(T)-F12/RI  DIRECT PMODEL")
                file_orca.write(h2_aux.h2cabs)
                if do_efield:
                    file_orca.write("% pointcharges \"ntc_charges.chg\"\n")

                file_orca.write(h2_aux.gen_orca_aux(1, 2800))
                for mol_c in range(0, 2):
                    file_orca.write("H : {0:10.4f} {1:10.4f} {2:10.4f}\n".format(
                        Hmol1[mol_c][0], Hmol1[mol_c][1], Hmol1[mol_c][2]))

                for mol_c in range(0, 2):
                    file_orca.write("H  {0:10.4f} {1:10.4f} {2:10.4f}\n".format(
                        Hmol2[mol_c][0], Hmol2[mol_c][1], Hmol2[mol_c][2]))

                file_orca.write("*\n\n")

                """
                print(4+len(chg))
                print(R,T,TT,P)
                print("H  ", *Hmol1[0])
                print("H  ", *Hmol1[1])
                print("H  ", *Hmol2[0])
                print("H  ", *Hmol2[1])
                [print("X ", *c) for c in chg]
                """
                P_enter_again = True

            file_orca.close()

print("Counter = {0}".format(counter))
