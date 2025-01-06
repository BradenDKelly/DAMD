




def read_cnf(input="cnf_input.inp"):
    r = []
    e = []
    i = 0
    box1 = 0.0
    #box=9.42953251

    with open(input, 'r') as file:
        for line in file:
            i += 1

            if i == 2:  #length(split(line)) == 1  && split(line) != typeof(Flo)
                box1 = float(line.strip())  # 9.42953251 #parse(Float64, split(line)) # this should be on the 2nd Line
                #print("hardcoded box at line 257 in ReadCNF: ", box1)

            if i >= 3:  #length(split(line)) > 4
                lin = line.split()

                r.append([
                    float(lin[0].strip()),
                    float(lin[1].strip()),
                    float(lin[2].strip())
                ])

                e.append([
                    float(lin[3].strip()),
                    float(lin[4].strip()),
                    float(lin[5].strip()),
                    float(lin[6].strip())
                ])

    return r, e, box1


def read_nist(filename):
    lj_coords = []
    qq_coords = []
    charge = []
    atom_type = []
    atom_name = []
    mol_name = []
    num = 7
    print(filename)
    i = 0
    box1 = 0

    with open(filename, 'r') as file:
        for line in file:
            i += 1
            if i == 1:  # length(split(line)) == 1  && split(line) != typeof(Flo)
                box1 = float(line.split()[0])  # 9.42953251

            if len(line.strip().split()) > 2 and i > 2:
                line = line.strip()
                split_line = line.split()
                coords = [float(split_line[1]), float(split_line[2]), float(split_line[3])]
                qq_coords.append(coords)
                mol_name.append("WAT")

                if split_line[4] == "O":
                    num = 7
                    atom_name.append("O1")
                    atom_type.append(1)
                    charge.append(-2 * 0.42380)
                    lj_coords.append(coords)
                elif split_line[4] == "H":
                    num += 1
                    atom_name.append(f"H{num}")
                    atom_type.append(2)
                    charge.append(0.42380)

    qq_q = [charge[i] for i in range(len(charge))]
    qq_r = [tuple(qq_coords[i]) for i in range(len(qq_coords))]  # Using tuple instead of SVector

    return qq_r, qq_q, rm, qq_r, atom_tracker, box1, atom_name, atom_type


def init_cubic_grid(n: int, rho: float):
    """
    ------------------------------------------------------------------------
    Created by Braden Kelly
    ------------------------------------------------------------------------
    Creates an initial configuration
    ------------------------------------------------------------------------
    input:  n       number of particles
           rho     density
    output: coords  coordinates of n particles
    ------------------------------------------------------------------------
    """
    # Calculate box length (L)
    L = (n / rho) ** (1.0 / 3.0)

    # Calculate the lowest perfect cube that will contain all of the particles
    n_cube = 2
    while (n_cube**3 < n):
        n_cube += 1

    coords = numpy.zeros((3, n))
    # initial position of particle 1
    posit = [0, 0, 0]

    # begin assigning particle positions
    for i in range(n):
        coords[:, i] = numpy.array([posit[0] + 0.01, posit[1] + 0.01, posit[2] + 0.01]) * (L / n_cube)
        # Advancing the index (posit)
        posit[0] += 1
        if posit[0] == n_cube:
            posit[0] = 0
            posit[1] += 1
            if posit[1] == n_cube:
                posit[1] = 0
                posit[2] += 1

    # Convert to list of 3D vectors
    return [numpy.array([coords[0, i], coords[1, i], coords[2, i]]) for i in range(coords.shape[1])]

