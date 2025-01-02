import re
import numpy as np
from collections import Counter
from typing import List, Tuple

class Topology:
    def __init__(self, mol_name: str, box: List[float], coords: List[List[float]],
                 atomnm: List[str], resnm: List[str], resnr: List[int], elem: List[str]):
        self.mol_name = mol_name
        self.box = box
        self.coords = coords
        self.atomnm = atomnm
        self.resnm = resnm
        self.resnr = resnr
        self.elem = elem

    def __str__(self):
        return (f"Topology: {self.mol_name}\n"
                f"Box Dimensions: {self.box}\n"
                f"Number of Atoms: {len(self.coords)}\n"
                f"Residues: {Counter(self.resnm)}")

def read_pdb(pdb_name: str) -> Topology:
    mol_name = pdb_name.split(".")[0]
    box = [None] * 3
    coords, atomnm, resnm, resnr, elem = [], [], [], [], []
    df = 0.1

    try:
        with open(pdb_name, 'r') as file:
            for line in file:
                if line.startswith(("ATOM", "HETATM")):
                    coords.append([
                        df * float(line[30:38].strip()),
                        df * float(line[39:46].strip()),
                        df * float(line[47:54].strip())
                    ])
                    atomnm.append(line[12:15].strip())
                    resnm.append(line[17:21].strip())
                    resnr.append(int(line[22:27].strip()))
                    elem.append(line[76:78].strip() if len(line) >= 77 else "")
                elif line.startswith("CRYST1"):
                    box = [
                        df * float(line[6:15].strip()),
                        df * float(line[15:24].strip()),
                        df * float(line[24:33].strip())
                    ]
    except FileNotFoundError:
        raise ValueError(f"File {pdb_name} not found.")
    except Exception as e:
        raise ValueError(f"An error occurred while reading {pdb_name}: {e}")

    residue_counts = Counter(resnm)
    if len(residue_counts) == 1:
        print(f"{list(residue_counts.keys())[0]} is alone in {pdb_name}")
    else:
        print(f"Residues in {pdb_name}: {residue_counts}")

    return Topology(mol_name, box, coords, atomnm, resnm, resnr, elem)

