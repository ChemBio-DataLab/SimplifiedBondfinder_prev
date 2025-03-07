import os
import re

def extract_gibbs_free_energy_from_log(log_file_path):
    """
    Extracts the "Sum of electronic and thermal Free Energies" value from a Gaussian log file.

    Args:
        log_file_path (str): Path to the Gaussian log file.

    Returns:
        tuple: File name and the extracted Gibbs free energy value in Hartree.
    """
    # Check if the log file exists
    if not os.path.exists(log_file_path):
        raise FileNotFoundError(f"The file {log_file_path} does not exist.")

    # Initialize variable to store the free energy value
    gibbs_free_energy = None

    # Read the log file and search for the relevant line
    with open(log_file_path, 'r') as log_file:
        for line in log_file:
            if "Sum of electronic and thermal Free Energies" in line:
                # Extract the energy value using regex
                match = re.search(r"-\d+\.\d+", line)
                if match:
                    gibbs_free_energy = float(match.group(0))
                break

    # Check if the value was found
    if gibbs_free_energy is None:
        raise ValueError("Could not find the 'Sum of electronic and thermal Free Energies' value in the log file.")

    return os.path.basename(log_file_path), gibbs_free_energy

def process_all_log_files_in_directory(directory):
    """
    Processes all *.log files in a directory, extracts Gibbs free energy values, and saves them to an output file.

    Args:
        directory (str): Path to the directory containing *.log files.

    Returns:
        None
    """
    output_file_path = os.path.join(directory, "extracted_gibbs_energies.txt")

    with open(output_file_path, 'w') as output_file:
        output_file.write("File Name\tGibbs Free Energy (Hartree)\n")
        output_file.write("---------------------------------------\n")

        for file_name in os.listdir(directory):
            if file_name.endswith(".log"):
                log_file_path = os.path.join(directory, file_name)
                try:
                    file_name, gibbs_free_energy = extract_gibbs_free_energy_from_log(log_file_path)
                    output_file.write(f"{file_name}\t{gibbs_free_energy}\n")
                except (FileNotFoundError, ValueError) as e:
                    print(f"Skipping {file_name}: {e}")

    print(f"Extracted Gibbs free energies saved to {output_file_path}")

# Example usage
if __name__ == "__main__":
    current_directory = os.getcwd()
    process_all_log_files_in_directory(current_directory)

