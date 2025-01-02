import os

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
    def __init__(self, name, nrexcl, atoms, bonds=None, angles=None, dihedrals=None):
        self.name = name
        self.nrexcl = nrexcl
        self.atoms = atoms
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

def read_section(file, section_name):
    """
    Helper function to read a section of the file and parse it.
    """
    zone_data = []
    #with open(filename, 'r') as file
    for line in file:
        line = line.split(';')[0].strip()
        
        if not line:
            continue
        #if line.startswith(f'[{section_name}]'):
        if section_name in line:
            while True:
                line = file.readline().strip()
                if len(line.split(';')[0].strip()) <= 1:   # verify this isnt header of section
                    continue
                if not line or line.startswith('['):  # end of section or new section
                    break
                zone_data.append(line)
    return zone_data


def parse_atomtypes(atomtypes_data):
    atomtypes = []
    for line in atomtypes_data:
        s = line.split()
        atomtypes.append(Atomtypes(s[0], float(s[1]), float(s[2]), s[3], float(s[4]), float(s[5])))
        
    return atomtypes


def parse_molecules(molecules_data):
    molecules = {}
    for line in molecules_data:
        s = line.split()
        molecules[s[0]] = int(s[1])
    return molecules


def parse_bonds(bonds_data):
    bonds = []
    for line in bonds_data:
        s = line.split()
        bonds.append(Bonds(int(s[0]), int(s[1]), int(s[2]), float(s[3]), float(s[4])))
    return bonds


def parse_angles(angles_data):
    angles = []
    for line in angles_data:
        s = line.split()
        angles.append(Angles(int(s[0]), int(s[1]), int(s[2]), int(s[3]), float(s[4]), float(s[5])))
    return angles


def read_top_file(top_file):
    atomtypes = []
    molparams = []
    molecules = {}

    # If there are multiple molecules, we need to first get "defaults" and "atomtypes"
    # We then have to parse each molecule, delimited by "moleculetype"
    # We then end with "system" and "molecules"
    # Also, we can have itp files embedded that we have to look for.

        bonds_data = read_section(file, "bonds")
        angles_data = read_section(file, "angles")



        atomtypes = parse_atomtypes(atomtypes_data)
        molecules = parse_molecules(molecules_data)
        bonds = parse_bonds(bonds_data)
        angles = parse_angles(angles_data)
        # Process other sections if needed
        # Further parsing logic (e.g., [bonds], [angles], etc.)
        
    return molparams, molecules, atomtypes

# Example usage
#top_file = "your_topology_file.top"
#molparams, molecules, atomtypes = read_top_file(top_file)

