bohr2angs = 0.529177
ha2cminv = 220000.00000000003

d_N2 = 1.099998  
d_N2_SolidState = 1.06  

DEFAULT_ORCA_HEADER = """! RHF  DLPNO-CCSD(T) aug-cc-pVQZ aug-cc-pVQZ/C RIJCOSX def2/J VERYTIGHTSCF  DIRECT PMODEL
%maxcore 3000
%pal nprocs 3 
end"""

DLPNO_CCSDT_CBS_ORCH_HEADER = """! RHF DLPNO-CCSD(T) Extrapolate(4/5,aug-cc) aug-cc-pVQZ/C RIJCOSX def2/J VERYTIGHTSCF PMODEL
%maxcore 3000
%pal nprocs 4 
end
"""

CCSDT_CBS_ORCA_HEADER = """! RHF CCSD(T) Extrapolate(4/5,aug-cc) VERYTIGHTSCF PMODEL
%maxcore 19000
%pal nprocs 4 
end"""

CCSDT_CBS_ORCA_HEADER_TZP = """! RHF CCSD(T) aug-cc-pVTZ VERYTIGHTSCF PMODEL
%maxcore 3000
%pal nprocs 4 
end"""

CCSDT_CBS_ORCA_HEADER_QZP = """! RHF CCSD(T) aug-cc-pVQZ VERYTIGHTSCF PMODEL
%maxcore 3000
%pal nprocs 4 
end"""

CCSDT_CBS_RI_ORCA_HEADER = """! RHF CCSD(T) Extrapolate(4/5,aug-cc) VERYTIGHTSCF PMODEL
%maxcore 19000
%pal nprocs 4 
end"""

AUG_CCPVXZ_BASIS_NAMES = {
    2 : "aug-cc-pVDZ",
    3 : "aug-cc-pVTZ",
    4 : "aug-cc-pVQZ",
    5 : "aug-cc-pV5Z"
}

CCPVXZ_BASIS_NAMES = {
    2 : "cc-pVDZ",
    3 : "cc-pVTZ",
    4 : "cc-pVQZ",
    5 : "cc-pV5Z"
}

N2_Ar_DISTANCES_BOHR = [
    4.75, 5.00, 5.25, 5.50, 5.75, 6.00, 
    6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 
    7.75, 8.00, 8.25, 8.50, 8.75, 9.00,
    9.25, 9.50, 9.75, 10.0, 11.0, 12.0,
    14.0, 16.0, 18.0, 20.0, 25.0, 30.0
]

N2_Ar_DISTANCES = [r*bohr2angs for r in N2_Ar_DISTANCES_BOHR ]

LEGENDRE_ROOTS_ANGLE_ARCCOS = [
    2.963498195463516,
    2.732792416247772,
    2.500726327116406,
    2.268356043649630,
    2.035874587341303,
    1.803344977489339,
    1.570796326794897,
    1.338247676100454,
    1.105718066248490,
    0.8732366099401631,
    0.6408663264733868,
    0.4088002373420213,
    0.1780944581262768
]

LEGENDRE_ROOTS_ANGLE_WITH_ZERO = [
    2.963498195463516,
    2.732792416247772,
    2.500726327116406,
    2.268356043649630,
    2.035874587341303,
    1.803344977489339,
    1.570796326794897,
    1.338247676100454,
    1.105718066248490,
    0.8732366099401631,
    0.6408663264733868,
    0.4088002373420213,
    0.1780944581262768,
    0.0000000000000000
]

LEGENDRE_ROOTS_ANGLE_WITH_0_HALFPI = [
    1.570796326794897,
    1.338247676100454,
    1.105718066248490,
    0.8732366099401631,
    0.6408663264733868,
    0.4088002373420213,
    0.1780944581262768,
    0.0000000000000000
]