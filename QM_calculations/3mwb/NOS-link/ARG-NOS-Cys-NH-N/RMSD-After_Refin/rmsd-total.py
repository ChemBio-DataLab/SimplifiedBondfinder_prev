import numpy as np

def parse_gjf(filename):
    """
    Parses a Gaussian .gjf file to extract atomic symbols and coordinates.
    Returns a list of tuples: (symbol, np.array([x, y, z]))
    """
    atoms = []
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Variables to track the parsing state
    parsing_coords = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Detect the charge and multiplicity line (e.g., "0 1")
        if not parsing_coords:
            if stripped == '':
                continue  # Skip blank lines
            if stripped[0].isdigit() or (stripped[0] == '-' and stripped[1:].isdigit()):
                parts = stripped.split()
                if len(parts) >= 2:
                    charge, multiplicity = parts[:2]
                    # Start parsing coordinates from the next line
                    parsing_coords = True
                continue
        else:
            if stripped == '':
                break  # End of coordinates
            parts = stripped.split()
            if len(parts) < 4:
                continue  # Not enough parts to contain symbol and coordinates
            symbol = parts[0]
            try:
                # Assume that the last three parts are x, y, z
                x, y, z = map(float, parts[-3:])
                atoms.append((symbol, np.array([x, y, z])))
            except ValueError:
                # If conversion fails, skip the line
                continue
    
    if not atoms:
        print(f"Warning: No atoms were parsed from {filename}. Please check the file format.")
    
    return atoms

def kabsch(P, Q):
    """
    The Kabsch algorithm: calculates the optimal rotation matrix U
    that minimizes the RMSD between P and Q.
    P and Q should be numpy arrays of shape (N, 3)
    """
    # Compute the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Compute the optimal rotation using SVD
    V, S, Wt = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    U = np.dot(V, Wt)

    return U

def align_structures(refined_atoms, optimized_atoms):
    """
    Aligns optimized_atoms to refined_atoms using the Kabsch algorithm.
    Returns the aligned optimized_atoms coordinates.
    """
    # Extract coordinates
    P = np.array([atom[1] for atom in refined_atoms])
    Q = np.array([atom[1] for atom in optimized_atoms])

    if P.size == 0 or Q.size == 0:
        raise ValueError("One of the atom lists is empty. Cannot proceed with alignment.")

    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Center the points
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # Compute the optimal rotation matrix
    U = kabsch(Q_centered, P_centered)

    # Rotate and translate Q
    Q_rotated = np.dot(Q_centered, U)
    Q_aligned = Q_rotated + centroid_P

    return Q_aligned

def write_gjf(filename, header_lines, atoms, charge='0', multiplicity='1'):
    """
    Writes a .gjf file with given header and atoms.
    header_lines: list of strings (lines before the coordinates)
    atoms: list of tuples (symbol, np.array([x, y, z]))
    charge: string representing the charge
    multiplicity: string representing the multiplicity
    """
    with open(filename, 'w') as file:
        for line in header_lines:
            file.write(line)
        file.write('\n')  # Blank line
        file.write('Aligned structure\n\n')
        # Write charge and multiplicity
        file.write(f"{charge} {multiplicity}\n")
        for atom in atoms:
            symbol, coord = atom
            file.write(f"{symbol}  {coord[0]:>15.6f} {coord[1]:>15.6f} {coord[2]:>15.6f}\n")

def write_xyz(filename, atoms, comment=''):
    """
    Writes an .xyz file with given atoms.
    atoms: list of tuples (symbol, np.array([x, y, z]))
    comment: string for the comment line
    """
    with open(filename, 'w') as file:
        file.write(f"{len(atoms)}\n")
        file.write(f"{comment}\n")
        for atom in atoms:
            symbol, coord = atom
            file.write(f"{symbol}  {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}\n")

def extract_header(filename):
    """
    Extracts the header lines from a .gjf file (everything before the charge and multiplicity line).
    """
    header = []
    charge = '0'
    multiplicity = '1'
    with open(filename, 'r') as file:
        for line in file:
            stripped = line.strip()
            if stripped == '':
                continue  # Skip blank lines
            parts = stripped.split()
            if len(parts) >= 2 and (parts[0].isdigit() or (parts[0].startswith('-') and parts[0][1:].isdigit())):
                # Detected charge and multiplicity line; extract them and stop
                charge = parts[0]
                multiplicity = parts[1]
                break
            header.append(line)
    return header, charge, multiplicity

def calculate_rmsd(coords1, coords2):
    """
    Calculates the Root-Mean-Square Deviation (RMSD) between two sets of coordinates.
    coords1 and coords2 should be numpy arrays of shape (N, 3)
    """
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

def main():
    Refined_stru = 'After_Refine-NoH.gjf'
    optimized_gjf = '3mwb-NH2NOS-QM-optimized-NoH.gjf'
    output_gjf = 'aligned_QM-optimized.gjf'
    output_xyz = 'aligned_QM-optimized.xyz'
    Refined_xyz = 'After_Refine-NoH.xyz'

    # Parse the .gjf files
    refined_atoms = parse_gjf(Refined_stru)
    optimized_atoms = parse_gjf(optimized_gjf)

    print(f"Refined structure parsed: {len(refined_atoms)} atoms.")
    print(f"Optimized structure parsed: {len(optimized_atoms)} atoms.")

    # Check if the number of atoms and their symbols match
    if len(refined_atoms) != len(optimized_atoms):
        print("Error: The two structures have different numbers of atoms.")
        return
    for i, (atom1, atom2) in enumerate(zip(refined_atoms, optimized_atoms)):
        if atom1[0] != atom2[0]:
            print(f"Error: Atom type mismatch at position {i+1}: {atom1[0]} vs {atom2[0]}")
            return

    # Align the structures
    try:
        aligned_coords = align_structures(refined_atoms, optimized_atoms)
    except ValueError as ve:
        print(f"Error during alignment: {ve}")
        return

    # Prepare the aligned atoms list
    aligned_atoms = []
    for i, atom in enumerate(refined_atoms):
        symbol = atom[0]
        coord = aligned_coords[i]
        aligned_atoms.append((symbol, coord))

    # Extract header and charge/multiplicity from the refined .gjf file
    header, charge, multiplicity = extract_header(Refined_stru)

    # Write the aligned structure to a new .gjf file
    write_gjf(output_gjf, header, aligned_atoms, charge, multiplicity)
    print(f"Aligned structure has been written to {output_gjf}")

    # Write the aligned structure to an .xyz file
    write_xyz(output_xyz, aligned_atoms, comment='Aligned Optimized Structure')
    print(f"Aligned structure has been written to {output_xyz}")

    # Optionally, write the refined structure to an .xyz file for comparison
    write_xyz(Refined_xyz, refined_atoms, comment='Refined Structure')
    print(f"Refined structure has been written to {Refined_xyz}")

    # Calculate overall RMSD
    aligned_coords_array = np.array([atom[1] for atom in aligned_atoms])
    refined_coords_array = np.array([atom[1] for atom in refined_atoms])
    overall_rmsd = calculate_rmsd(refined_coords_array, aligned_coords_array)

    # Print only the total (overall) RMSD
    print(f"\nOverall RMSD: {overall_rmsd:.4f} Å")

if __name__ == "__main__":
    main()

