import math
import numpy as np


def transfer_masses(ff_params):
    """Transfer masses from atomtypes to atoms in the force field parameters."""
    
    # Create a lookup dictionary for atomtype masses
    atomtype_masses = {at.name: at.mass for at in ff_params.atomtypes}
    
    # Update masses in all molecules
    for i, atomtypes in enumerate(ff_params.atomtypes):
        found_match = False
        # for each molecules
        for j, params in enumerate(ff_params.molparams):
            # go through each atom, looking at the name and comparing with atomtypes[i] name
            for k, atoms in enumerate(params.atoms):
                # Get the mass from atomtypes based on the atom's type
                if atomtypes.name == atoms.type:
                    ff_params.atomtypes[i].mass = params.atoms_dict[atomtypes.name].mass
                    ff_params.atomtypes[i].charge = params.atoms_dict[atomtypes.name].charge
                    found_match = True
        if not found_match:
            raise ValueError(f"Atom type '{atomtypes.name}' not found in atoms")

    return ff_params
    
class FFParameters:
    def __init__(self, defaults, atomtypes, molparams, system, molecules):
        self.defaults=defaults
        self.atomtypes = atomtypes
        self.molparams = molparams
        self.system = system
        self.molecules = molecules
        
        
class Atomtypes:
    def __init__(self, name, mass, charge, ptype, sigma, epsilon):
        self.name = name
        #self.atomicnr = atomicnr
        self.mass = mass
        self.charge = charge
        self.ptype = ptype
        self.sigma = sigma
        self.epsilon = epsilon

class Defaults:
    def __init__(self, nbfunc, comb_rule, gen_pairs, fudgeLJ, fudgeQQ):
        self.nbfunc = nbfunc
        self.comb_rule = comb_rule
        self.gen_pairs = gen_pairs
        self.fudgeLJ = fudgeLJ
        self.fudgeQQ = fudgeQQ

class MolParam:
    def __init__(self, name, nrexcl, atoms, atoms_dict=None, bonds=None, angles=None, dihedrals=None):
        self.name = name
        self.nrexcl = nrexcl
        self.atoms = atoms
        self.atoms_dict = atoms_dict
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals

class Atoms:
    def __init__(self, nr, type, resnr, residue, atom, cgnr, charge, mass):
        self.nr = nr
        self.type = type
        self.resnr = resnr
        self.residue = residue
        self.atom = atom
        self.cgnr = cgnr
        self.charge = charge
        self.mass = mass

class Bonds:
    def __init__(self, ai, aj, funct, r, k):
        self.ai = ai
        self.aj = aj
        self.funct = funct
        self.r = r
        self.k = k

class Pairs:
    def __init__(self, ai, aj, funct):
        self.ai = ai
        self.aj = aj
        self.funct = funct

class Angles:
    def __init__(self, ai, aj, ak, funct, theta, cth):
        self.ai = ai
        self.aj = aj
        self.ak = ak
        self.funct = funct
        self.theta = theta
        self.cth = cth

class Dihedrals:
    def __init__(self, i, j, k, l, func, coeffs):
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.func = func
        self.coeffs = coeffs

def read_top_file(top_file: str):
    # Read forcefield and topology file
    # ffParameters = FFParameters()
    system, defaults = None, None  # these are temporary objects (structs) for this function
    atomtypes = []
    moleculetypes = []  # 1 object for each molecule type
    top_name = "?"
    top_nrexcl = 0
    molparams = []
    topatoms = []
    topatoms_dict = {}
    topbonds = []
    toppairs = []
    topangles = []
    topdihedrals = []
    zone = "?"
    molecules = {}  # equivalent to Dict{String,Int64}()
    mol_number = 0
    mol_number_check = 0
    mol_type = 0  # dummy counter

    # First, scan to the bottom and get the number of molecules in the system
    # double check is implemented, we count number of molecules defined along the way
    # Actually, we do the reverse of these two things :)
    with open(top_file) as file:
        for lin in file:
            line = lin.split(";")[0].strip()
            if "[ moleculetype ]" in line:
                mol_number += 1
            if "#include" in line:
                filename = line.split('"')[1]
                with open(filename) as file:
                    for lins in file:
                        lines = lins.split(";")[0].strip()
                        if "[ moleculetype ]" in lines:
                            mol_number += 1
            # Implement a double check on the number of molecule types in
            # topology file by checking the [ molecule ] section
            if "[ molecules ]" in line:
                zone = line[1:-1].strip()
            if zone == "molecules" and len(line.split()) == 2:
                mol_number_check += 1

    if mol_number != mol_number_check:
        raise ValueError(
            f"The number of moleculetypes in topology file {top_file} {mol_number} does not match the number in section [ molecules ] of said file {mol_number_check}."
        )

    ################################################################################
    #                  Read parameters from Topology File
    ################################################################################

    with open(top_file) as file:
        for lin in file:
            line = lin.split(";")[0].strip()  # I only care about what is on the LHS
            if len(line) == 0:  # of the first semi-colon
                continue
            if line.startswith('[') and line.endswith(']'):
                if (zone in ["moleculetype", "system"]) and (len(topatoms) > 0):
                    molparams.append(
                        MolParam(
                            top_name,
                            top_nrexcl,
                            topatoms[:],
                            topatoms_dict,
                            topbonds[:],
                            topangles[:],
                            topdihedrals[:],
                        )
                    )
                    topatoms = []
                    topatoms_dict = {}
                    topbonds = []
                    toppairs = []
                    topangles = []
                    topdihedrals = []
                    mol_type += 1
                elif (zone in ["moleculetype", "system"]) and len(topatoms) == 0:
                    mol_type += 1
                zone = line[1:-1].strip()
                continue
            if zone == "defaults" and len(line.split()) > 2:
                s = line.split()
                defaults = Defaults(
                    int(s[0]),
                    int(s[1]),
                    s[2],
                    float(s[3]),
                    float(s[4]),
                )
            elif zone == "atomtypes" and len(line.split()) > 2:
                inline = line.split()
                atomtypes.append(
                    Atomtypes(
                        inline[0],
                        #inline[1],
                        float(inline[1]),
                        float(inline[2]),
                        inline[3],
                        float(inline[4]),
                        float(inline[5]),
                    )
                )
            elif zone == "molecules":  # last item in topology
                inline = line.split()
                molecules[inline[0]] = int(inline[1])
            elif "#include" in line:
                zone = "itp"
                # open itp file and get atoms/bonds/angles/pairs/dihedrals
                filename = line.split('"')[1]  # name is in quotes, which is the 2nd segment ...     #include "name"
                itpzone = "?"  # using double quotation mark as delimiter -> segment:   (1)     (2)
                print(f"Reading Topology information from itp file for {filename}")
                with open(filename) as file:
                    # First store everything in lists, then add to appropriate
                    # objects as vectors once lengths are known
                    itp_name = "?"
                    itp_nrexcl = 0
                    itp_atoms = []
                    itp_atoms_dict = {}
                    itp_bonds = []
                    itp_pairs = []
                    itp_angles = []
                    itp_dihedrals = []
                    for itplin in file:  # read each line of file
                        itpline = itplin.split(";")[0].strip()  # I only care about what is on the LHS
                        if len(itpline) == 0:  # of the first semi-colon
                            continue
                        if itpline.startswith('[') and itpline.endswith(']'):
                            itpzone = itpline[1:-1].strip()  # there is a word between the brackets... get that word...
                            continue
                        if itpzone == "moleculetype" and len(itpline.split()) >= 1:  # Make new molecule object
                            inline = itpline.split()
                            itp_name = inline[0]
                            itp_nrexcl = int(inline[1])
                        elif itpzone == "atoms" and len(itpline.split()) > 2:
                            s = itpline.split()
                            itp_atoms.append(
                                Atoms(
                                    int(s[0]),
                                    s[1],
                                    int(s[2]),
                                    s[3],
                                    s[4],
                                    int(s[5]),
                                    float(s[6]),
                                    float(s[7]),
                                )
                            )
                            itp_atoms_dict[s[1]] = Atoms(
                                    int(s[0]),
                                    s[1],
                                    int(s[2]),
                                    s[3],
                                    s[4],
                                    int(s[5]),
                                    float(s[6]),
                                    float(s[7]),
                                )
                            
                        elif itpzone == "bonds" and len(itpline.split()) > 2:
                            s = itpline.split()
                            itp_bonds.append(
                                Bonds(
                                    int(s[0]),
                                    int(s[1]),
                                    int(s[2]),
                                    float(s[3]),
                                    float(s[4]),
                                )
                            )
                        elif itpzone == "pairs" and len(itpline.split()) > 2:
                            s = itpline.split()
                            itp_pairs.append(
                                Pairs(
                                    int(s[0]),
                                    int(s[1]),
                                    int(s[2]),
                                )
                            )
                        elif itpzone == "angles" and len(itpline.split()) > 2:
                            s = itpline.split()
                            itp_angles.append(
                                Angles(
                                    int(s[0]),
                                    int(s[1]),
                                    int(s[2]),
                                    int(s[3]),
                                    math.radians(float(s[4])),
                                    float(s[5]),
                                )
                            )
                        elif itpzone == "dihedrals" and len(itpline.split()) > 2:
                            s = itpline.split()
                            if len(s) > 9:  # improper dihedral
                                itp_dihedrals.append(
                                    Dihedrals(
                                        int(s[0]),
                                        int(s[1]),
                                        int(s[2]),
                                        int(s[3]),
                                        int(s[4]),
                                        [float(s[i]) for i in range(5, 11)],
                                    )
                                )
                            else:  # proper dihedral
                                itp_dihedrals.append(
                                    Dihedrals(
                                        int(s[0]),
                                        int(s[1]),
                                        int(s[2]),
                                        int(s[3]),
                                        int(s[4]),
                                        float(s[5]),
                                        float(s[6]),
                                        int(s[7]),
                                    )
                                )
                    molparams.append(
                        MolParam(
                            itp_name,
                            itp_nrexcl,
                            itp_atoms[:],
                            itp_atoms_dict,
                            itp_bonds[:],
                            itp_angles[:],
                            itp_dihedrals[:],
                        )
                    )
            elif zone == "moleculetype":
                inline = line.split()
                top_name = inline[0]  # instantiate a molecule object and add to the type list
                top_nrexcl = int(inline[1])
            elif zone == "atoms" and len(line.split()) > 3:
                s = line.split()
                temp_atoms = Atoms(
                    int(s[0]),
                    s[1],
                    int(s[2]),
                    s[3],
                    s[4],
                    int(s[5]),
                    float(s[6]),
                    float(s[7]),
                )
                topatoms.append(temp_atoms)

                # make dictionary too for easier lookup later
                topatoms_dict[s[1]] = temp_atoms
                
            elif zone == "bonds" and len(line.split()) > 2:
                s = line.split()
                topbonds.append(
                    Bonds(
                        int(s[0]),
                        int(s[1]),
                        int(s[2]),
                        float(s[3]),
                        float(s[4]),
                    )
                )
            elif zone == "pairs" and len(line.split()) > 2:
                s = line.split()
                toppairs.append(
                    Pairs(
                        int(s[0]),
                        int(s[1]),
                        int(s[2]),
                    )
                )
            elif zone == "angles" and len(line.split()) > 2:
                s = line.split()
                topangles.append(
                    Angles(
                        int(s[0]),
                        int(s[1]),
                        int(s[2]),
                        int(s[3]),
                        math.radians(float(s[4])),
                        float(s[5]),
                    )
                )
            elif zone == "dihedrals" and len(line.split()) > 2:
                s = line.split()
                if len(s) > 9:  # improper dihedral
                    topdihedrals.append(
                        Dihedrals(
                            int(s[0]),
                            int(s[1]),
                            int(s[2]),
                            int(s[3]),
                            int(s[4]),
                            [float(s[i]) for i in range(5, 11)],
                        )
                    )
                else:  # proper dihedral
                    topdihedrals.append(
                        Dihedrals(
                            int(s[0]),
                            int(s[1]),
                            int(s[2]),
                            int(s[3]),
                            int(s[4]),
                            float(s[5]),
                            float(s[6]),
                            int(s[7]),
                        )
                    )
            elif zone == "system":  # second last item in topology
                system = line.strip()
                print(f"System name is: {system}")

    system_topology = FFParameters(
        defaults,
        atomtypes[:],
        molparams[:],
        system,
        molecules,
    )
    
    return transfer_masses(system_topology) #system_topology
