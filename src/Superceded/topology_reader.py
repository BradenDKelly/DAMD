def read_top_file(top_file):
    """
    This function reads a topology file and extracts relevant topology information.
    The data is stored in appropriate objects for processing.
    """
    import re

    # Initialize data structures
    system, defaults = None, None
    atomtypes = []
    moleculetypes = []
    top_name = "?"
    top_nr_excl = 0
    molparams = []
    topatoms, topbonds, toppairs, topangles, topdihedrals = [], [], [], [], []
    zone = "?"
    molecules = {}
    mol_number, mol_number_check = 0, 0
    mol_type = 0  # Dummy counter

    # First pass: Count molecules and validate the topology file
    with open(top_file, 'r') as file:
        for line in file:
            line = line.split(";")[0].strip()
            if "[ moleculetype ]" in line:
                mol_number += 1
            elif "#include" in line:
                filename = re.search(r'"(.*?)"', line).group(1)
                with open(filename, 'r') as incl_file:
                    for incl_line in incl_file:
                        if "[ moleculetype ]" in incl_line:
                            mol_number += 1
            if "[ molecules ]" in line:
                zone = line[1:-2].strip()
            if zone == "molecules" and len(line.split()) == 2:
                mol_number_check += 1
    print(f"Number of mol_number {mol_number} and mol_number_check {mol_number_check}")
    if mol_number != mol_number_check:
        raise ValueError(
            f"The number of molecule types in {top_file} does not match the count in the [ molecules ] section."
        )

    # Second pass: Extract topology parameters
    with open(top_file, 'r') as file:
        for line in file:
            
            line = line.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                if zone in {"moleculetype", "system"} and len(topatoms) > 0:
                    molparams.append({
                        "name": top_name,
                        "nr_excl": top_nr_excl,
                        "atoms": topatoms[:],
                        "bonds": topbonds[:],
                        "angles": topangles[:],
                        "dihedrals": topdihedrals[:]
                    })
                    topatoms, topbonds, toppairs, topangles, topdihedrals = [], [], [], [], []
                    mol_type += 1
                zone = line[1:-1].strip()
                continue
            if zone == "defaults" and len(line.split()) > 2:
                parts = line.split()
                defaults = {
                    "nbfunc": int(parts[0]),
                    "comb_rule": int(parts[1]),
                    "gen_pairs": parts[2],
                    "fudgeLJ": float(parts[3]),
                    "fudgeQQ": float(parts[4])
                }
            elif zone == "atomtypes" and len(line.split()) > 2:
                parts = line.split()
                atomtypes.append({
                    "name": parts[0],
                    #"atomicnr": parts[1],
                    "mass": float(parts[1]),
                    "charge": float(parts[2]),
                    "ptype": parts[3],
                    "sigma": float(parts[4]),
                    "epsilon": float(parts[5])
                })
            elif zone == "molecules":
                parts = line.split()
                molecules[parts[0]] = int(parts[1])
            elif "#include" in line:
                filename = re.search(r'"(.*?)"', line).group(1)
                print(f"Reading topology information from itp file: {filename}")
                # Additional processing for included files can be added here
                # ...

    return {
        "system": system,
        "defaults": defaults,
        "atomtypes": atomtypes,
        "moleculetypes": moleculetypes,
        "molecules": molecules,
        "molparams": molparams
    }