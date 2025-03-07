def gjf_to_xyz(gjf_file_path, xyz_file_path):
    """
    Converts a Gaussian .gjf file to an .xyz file.
    
    Parameters
    ----------
    gjf_file_path : str
        Path to the input .gjf file.
    xyz_file_path : str
        Path to the output .xyz file.
    """
    geometry_data = []
    with open(gjf_file_path, 'r') as gjf_file:
        for line in gjf_file:
            line = line.strip()
            # Skip empty lines or lines not starting with an atom symbol
            if not line or not line[0].isalpha():
                continue
            
            parts = line.split()
            # Expecting Atom, FixFlag (or Charge), X, Y, Z
            if len(parts) >= 5:
                atom_symbol = parts[0]
                try:
                    x, y, z = map(float, parts[-3:])
                    geometry_data.append((atom_symbol, x, y, z))
                except ValueError:
                    continue  # Skip lines that don't have valid coordinates

    if not geometry_data:
        raise ValueError(f"No valid atomic data found in {gjf_file_path}")

    # Write to XYZ file
    with open(xyz_file_path, 'w') as xyz_file:
        xyz_file.write(f"{len(geometry_data)}\n")
        xyz_file.write(f"Converted from {gjf_file_path}\n")
        for atom_symbol, x, y, z in geometry_data:
            xyz_file.write(f"{atom_symbol:2s} {x:15.6f} {y:15.6f} {z:15.6f}\n")


# Example usage
gjf_file_path = "6pgd-NHOS_QM-optimized.gjf"  # Replace with your .gjf file path
xyz_file_path = "6pgd-NHOS_QM-optimized.xyz"  # Replace with your desired .xyz file path

gjf_to_xyz(gjf_file_path, xyz_file_path)
print(f"Converted {gjf_file_path} to {xyz_file_path}.")

