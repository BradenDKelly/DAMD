import math
import numpy as np


def transfer_masses(ff_params):
    """Transfer masses from atomtypes to atoms in the force field parameters."""
    
  
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

def parse_molecule_section(lines, start_index=0):
    """Parse molecule section from a list of lines."""
    mol_name = "?"
    mol_nrexcl = 0
    atoms = []
    atoms_dict = {}
    bonds = []
    pairs = []
    angles = []
    dihedrals = []
    
    current_zone = "?"
    
    for line in lines[start_index:]:
        line = line.split(";")[0].strip()
        if not line:
            continue
            
        if line.startswith('[') and line.endswith(']'):
            current_zone = line[1:-1].strip()
            if current_zone not in ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals"]:
                break
            continue
            
        if current_zone == "moleculetype" and len(line.split()) >= 1:
            inline = line.split()
            mol_name = inline[0]
            mol_nrexcl = int(inline[1])
        elif current_zone == "atoms" and len(line.split()) > 2:
            s = line.split()
            atom = Atoms(
                int(s[0]), s[1], int(s[2]), s[3], s[4],
                int(s[5]), float(s[6]), float(s[7])
            )
            atoms.append(atom)
            atoms_dict[s[1]] = atom
        elif current_zone == "bonds" and len(line.split()) > 2:
            s = line.split()
            bonds.append(Bonds(
                int(s[0]), int(s[1]), int(s[2]),
                float(s[3]), float(s[4])
            ))
        elif current_zone == "pairs" and len(line.split()) > 2:
            s = line.split()
            pairs.append(Pairs(
                int(s[0]), int(s[1]), int(s[2])
            ))
        elif current_zone == "angles" and len(line.split()) > 2:
            s = line.split()
            angles.append(Angles(
                int(s[0]), int(s[1]), int(s[2]), int(s[3]),
                math.radians(float(s[4])), float(s[5])
            ))
        elif current_zone == "dihedrals" and len(line.split()) > 2:
            s = line.split()
            if len(s) > 9:  # improper dihedral
                dihedrals.append(Dihedrals(
                    int(s[0]), int(s[1]), int(s[2]), int(s[3]),
                    int(s[4]), [float(s[i]) for i in range(5, 11)]
                ))
            else:  # proper dihedral
                dihedrals.append(Dihedrals(
                    int(s[0]), int(s[1]), int(s[2]), int(s[3]),
                    int(s[4]), float(s[5]), float(s[6]), int(s[7])
                ))
                
    return MolParam(mol_name, mol_nrexcl, atoms, atoms_dict, bonds, angles, dihedrals), current_zone

def count_molecules(filename):
    """Count the number of unique molecule types in a file."""
    count = 0
    seen_molecules = set()  # To track unique molecule names
    with open(filename) as file:
        for line in file:
            line = line.split(";")[0].strip()
            if "[ moleculetype ]" in line:
                # Read the next line to get the molecule name
                next_line = next(file).strip()
                if next_line and not next_line.startswith(";"):
                    molecule_name = next_line.split()[0]  # Get the name of the molecule
                    if molecule_name not in seen_molecules:
                        seen_molecules.add(molecule_name)
                        count += 1
            if "#include" in line:
                included_file = line.split('"')[1]
                count += count_molecules(included_file)  # Recursively count in included files
    return count

def read_top_file(top_file: str):
    # Initialize containers
    system, defaults = None, None
    atomtypes = []
    molparams = []
    molecules = {}
    
    # Verify molecule count consistency
    mol_count = count_molecules(top_file)
    mol_count_check = 0
    
    with open(top_file) as file:
        lines = file.readlines()
        
    current_zone = "?"
    i = 0
    while i < len(lines):
        line = lines[i].split(";")[0].strip()
        if not line:
            i += 1
            continue
            
        if line.startswith('[') and line.endswith(']'):
            current_zone = line[1:-1].strip()
            i += 1
            continue
            
        if "#include" in line:
            filename = line.split('"')[1]
            with open(filename) as included_file:
                included_lines = included_file.readlines()
                # Parse the included file for molecules
                while included_lines:
                    mol_param, new_zone = parse_molecule_section(included_lines)
                    if mol_param.name != "?":
                        molparams.append(mol_param)
                        mol_count_check += 1  # Increment count for each molecule parsed
                    included_lines = included_lines[len(mol_param.atoms) + 1:]  # Move past the parsed section
            i += 1
            continue
            
        if current_zone == "defaults" and len(line.split()) > 2:
            s = line.split()
            defaults = Defaults(
                int(s[0]), int(s[1]), s[2],
                float(s[3]), float(s[4])
            )
        elif current_zone == "atomtypes" and len(line.split()) > 2:
            inline = line.split()
            atomtypes.append(Atomtypes(
                inline[0], float(inline[1]), float(inline[2]),
                inline[3], float(inline[4]), float(inline[5])
            ))
        elif current_zone == "molecules":
            inline = line.split()
            if len(inline) == 2:
                molecules[inline[0]] = int(inline[1])
                mol_count_check += 1  # Increment count for each molecule defined in this section
        elif current_zone == "moleculetype":
            mol_param, new_zone = parse_molecule_section(lines, i)
            if mol_param.name != "?":
                molparams.append(mol_param)
                mol_count_check += 1  # Increment count for each molecule parsed
            current_zone = new_zone
        
        i += 1
            
    if mol_count != mol_count_check:
        raise ValueError(
            f"The number of molecule types in topology file {top_file} ({mol_count}) "
            f"does not match the number in section [ molecules ] of said file ({mol_count_check})."
        )
        
    system_topology = FFParameters(
        defaults, atomtypes, molparams, system, molecules
    )
    
    return transfer_masses(system_topology)
