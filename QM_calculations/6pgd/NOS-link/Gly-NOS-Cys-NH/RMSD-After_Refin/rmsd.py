import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calculate_rmsd(xyz_file1, xyz_file2, visualize=True):
    """
    Aligns two XYZ molecular structures, calculates RMSD, and visualizes them.
    
    Parameters:
    - xyz_file1: Path to the first XYZ file (reference structure)
    - xyz_file2: Path to the second XYZ file (target structure)
    - visualize: Whether to show a 3D plot of the aligned molecules
    
    Returns:
    - RMSD value (float)
    """

    # Load the XYZ files
    u1 = mda.Universe(xyz_file1, format="XYZ")
    u2 = mda.Universe(xyz_file2, format="XYZ")

    # Align the second structure to the first
    align.alignto(u2, u1, select="all")

    # Get aligned coordinates
    coords1 = u1.atoms.positions
    coords2 = u2.atoms.positions

    # Compute RMSD
    rmsd_value = np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))

    print(f"RMSD between the two structures: {rmsd_value:.3f} Å")

    # Visualization
    if visualize:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        
        ax.scatter(coords1[:, 0], coords1[:, 1], coords1[:, 2], c="red", label="Reference (XYZ1)")
        ax.scatter(coords2[:, 0], coords2[:, 1], coords2[:, 2], c="blue", alpha=0.6, label="Aligned (XYZ2)")

        ax.set_xlabel("X-axis (Å)")
        ax.set_ylabel("Y-axis (Å)")
        ax.set_zlabel("Z-axis (Å)")
        ax.set_title("Aligned Molecular Structures")
        ax.legend()
        plt.show()

    return rmsd_value

# Example usage
xyz_file1 = "6pgd-NHOS_QM-optimized.xyz"  # Reference structure
xyz_file2 = "After_Refine.xyz"  # Target structure

rmsd = calculate_rmsd(xyz_file1, xyz_file2)

