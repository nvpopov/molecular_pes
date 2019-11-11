#!/usr/bin/python3

from math import pi,sin,cos,sqrt,pow
from enum import Enum
import sys, os
import itertools

d_N2 = 1.099998  
d_N2_SolidState = 1.06  

DEFAULT_ORCA_HEADER = """! RHF  DLPNO-CCSD(T) aug-cc-pVQZ aug-cc-pVQZ/C RIJCOSX def2/J VERYTIGHTSCF  DIRECT PMODEL
%maxcore 3000
%pal nprocs 3 
end"""

class ModelType(Enum):
    UNKNOWN = 1
    LINEAR_MOLECULE_VERSUS_SINGLE_ATOM = 2
    LINEAR_MOLECULE_VERSUS_LINEAR_MOLECULE = 3

class Atom:
    def __init__(self, atom_name, atom_coord):
        self.name = atom_name
        self.coord = atom_coord

class Molecule:
    def __init__(self):
        self.atoms = []

    def add_atom(self, atom_name, atom_coord):
        self.atoms.append(Atom(atom_name, atom_coord))

"""
  cb - create batch
"""
def cb_unknown(linear_molecule, single_atom, outdir, batch_per_input):
    raise Exception("Unknown model format")

def cb_linear_mol_vs_single_atom(linear_molecule, single_atom, outdir, batch_per_input,
                                 R0 = 1.7, Rmax = 10, NR = 50, phi0 = 0e0, phiM = pi/2e0, Nphi = 15):

    # phi0 = 0e0      # angle
    # phiM = pi/2e0   # angle
    # Nphi = 15
    Dphi = (phiM - phi0) / Nphi

    # R0 = 1.7        # Angstrom
    # Rmax = 10       # Angstrom
    # NR = 50        
    DR = (Rmax - R0) / NR

    #compose tasks root directory

    task_root_dir = "{}/batch".format(outdir)
    os.mkdir(task_root_dir)

    range_R = [i for i in range(NR)]
    range_phi = [i for i in range(Nphi)]

    batch_id = 0
    batch_file = None
    batch_file_domain = None

    xyz_file = open("anim.xyz", "w")

    for q, i in enumerate(itertools.product(range_R, range_phi)):

        if (q % batch_per_input == 0):
            batch_dir = "{}/batch/t_{}".format(outdir, batch_id)
            os.mkdir(batch_dir)

            batch_file_path = "{}/orca.inp".format(batch_dir)
            batch_file_domain_path = "{}/DATA".format(batch_dir)
            batch_file_domain = open(batch_file_domain_path, "w")
            batch_file = open(batch_file_path, "w")

            batch_id += 1

        ir = i[0]
        iphi = i[1]

        xyz_file.write("{}\n\n".format(len(linear_molecule.atoms)+1))
        #compose task directory
       
        #N1 = [d_N2/2, 0, 0]
        #N2 = [-d_N2/2, 0, 0]
        
        #Emit dimer AB
        if (batch_per_input > 1):
            pass
        batch_file.write("{}\n".format(DEFAULT_ORCA_HEADER))
        batch_file.write("* xyz 0 1\n")
        for atom in linear_molecule.atoms:
            batch_file.write("{} {} {} {}\n".format(atom.name, atom.coord[0], atom.coord[1], atom.coord[2]))
            xyz_file.write("{} {} {} {}\n".format(atom.name, atom.coord[0], atom.coord[1], atom.coord[2]))
        
        SA_R = R0 + ir*DR
        SA_phi = phi0 + iphi*Dphi
        SA = [SA_R * cos(SA_phi), SA_R * sin(SA_phi),0]
        
        dists_to_lm = []
        for atom_lm in linear_molecule.atoms:
            dist = pow(atom_lm.coord[0] - SA_R * cos(SA_phi), 2) + \
            pow(atom_lm.coord[1] - SA_R * sin(SA_phi), 2) + \
            pow(atom_lm.coord[2] - 0, 2) 
            dist_sqrt = sqrt(dist)
            dists_to_lm.append(dist_sqrt)

        dists_to_lm_str = " ".join(map(str, dists_to_lm))

        batch_file_domain.write("{} {} {}\n".format(SA_R, SA_phi, dists_to_lm_str))

        SA = [SA_R * cos(SA_phi), SA_R * sin(SA_phi),0]
        batch_file.write("{} {} {} {}\n".format(single_atom.atoms[0].name, SA[0], SA[1], SA[2]))
        xyz_file.write("{} {} {} {}\n".format(single_atom.atoms[0].name, SA[0], SA[1], SA[2]))
        batch_file.write("end\n\n")

        #Emit monomer A with ghost atoms(B)
        batch_file.write("$new_job\n")
        batch_file.write("{}\n".format(DEFAULT_ORCA_HEADER))
        batch_file.write("* xyz 0 1\n")
        for atom in linear_molecule.atoms:
            batch_file.write("{} {} {} {}\n".format(atom.name, atom.coord[0], atom.coord[1], atom.coord[2]))
     
        batch_file.write("{}: {} {} {}\n".format(single_atom.atoms[0].name, SA[0], SA[1], SA[2]))
        batch_file.write("end\n\n")

        #Emit monomer B with ghost atoms(A)
        batch_file.write("$new_job\n")
        batch_file.write("{}\n".format(DEFAULT_ORCA_HEADER))
        batch_file.write("* xyz 0 1\n")
        for atom in linear_molecule.atoms:
            batch_file.write("{}: {} {} {}\n".format(atom.name, atom.coord[0], atom.coord[1], atom.coord[2]))
        batch_file.write("{} {} {} {}\n".format(single_atom.atoms[0].name, SA[0], SA[1], SA[2]))
        batch_file.write("end\n\n")

    pass

def sample_tasks_N2_Ar_pes(N2dist):
    N2_mol = Molecule()
    N2_mol.add_atom("N", [N2dist/2, 0, 0])
    N2_mol.add_atom("N", [-N2dist/2, 0, 0])

    Ar_mol = Molecule()
    Ar_mol.add_atom("Ar", [0, 0, 0])

    cb_linear_mol_vs_single_atom(N2_mol, Ar_mol, "./", 1)

def sample_tasks_N2_Ne_pes(N2dist):
    N2_mol = Molecule()
    N2_mol.add_atom("N", [N2dist/2, 0, 0])
    N2_mol.add_atom("N", [-N2dist/2, 0, 0])

    Ne_mol = Molecule()
    Ne_mol.add_atom("Ne", [0, 0, 0])

    cb_linear_mol_vs_single_atom(N2_mol, Ne_mol, "./", 1)

def app_main():
    sample_tasks_N2_Ne_pes(d_N2)
    pass

if __name__ == "__main__":
    app_main()