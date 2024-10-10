import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

def generate_xyz_from_smiles(smiles_code, output_xyz_file="generated_structure.xyz"):
    try:
        subprocess.run(['obabel', f'-:{smiles_code}', '--gen3D', '-O', output_xyz_file], check=True)
        print(f"XYZ structure generated from SMILES and saved as {output_xyz_file}.")
        return output_xyz_file
    except subprocess.CalledProcessError as e:
        print(f"Error generating XYZ file from SMILES: {e}")
        return None

def optimize_geometry(input_xyz_file, output_xyz_file="optimized_structure.xyz"):
    xtb_path = 'xtb'
    try:
        subprocess.run([xtb_path, input_xyz_file, '--opt', '--gfn2', '--disp'], check=True)
        subprocess.run(['cp', 'xtbopt.xyz', output_xyz_file], check=True)
        print(f"Geometry optimization completed successfully. Optimized structure saved as {output_xyz_file}.")
        return output_xyz_file
    except subprocess.CalledProcessError as e:
        print(f"Error during geometry optimization: {e}")
        return None

def calculate_single_molecule_energy(optimized_xyz_file, conda_env_path="/path/to/conda/env/bin"):
    try:
        env = {"PATH": conda_env_path + ":" + os.environ["PATH"]}
        result = subprocess.run(['xtb', optimized_xyz_file, '--gfn2', '--sp', '--disp'], check=True, env=env, capture_output=True, text=True)
        energy_line = next(line for line in result.stdout.splitlines() if 'SCC energy' in line)
        single_molecule_energy = float(energy_line.split()[3])
        print(f"Single molecule SCC energy: {single_molecule_energy:.6f} Eh")
        return single_molecule_energy
    except subprocess.CalledProcessError as e:
        print(f"Error calculating single molecule energy: {e}")
        return None



def stack_and_optimize_dimer(optimized_xyz_file, stack_distance, conda_env_path="/path/to/conda/env/bin", output_file="optimized_dimer.xyz"):
    """
    Stacks two helicenes at a specified distance, optimizes the dimer geometry, and saves the optimized structure.

    Parameters:
    - optimized_xyz_file (str): Path to the optimized XYZ file of the single molecule.
    - stack_distance (float): Distance to stack the second molecule above the first one.
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.
    - output_file (str): Path where the optimized dimer structure will be saved.

    Returns:
    - str: The path to the optimized dimer structure file.
    """
    try:
        # Step 1: Stack the molecules at the specified distance
        stacked_xyz_file = "stacked_dimer.xyz"
        with open(optimized_xyz_file, 'r') as f:
            lines = f.readlines()

        num_atoms = int(lines[0].strip())
        atom_data = lines[2:]

        # Extract coordinates of the single molecule
        coords = []
        for line in atom_data:
            parts = line.split()
            element = parts[0]
            x, y, z = map(float, parts[1:4])
            coords.append((element, x, y, z))

        # Stack the second molecule at the specified distance
        stacked_coords = coords.copy()
        for element, x, y, z in coords:
            stacked_coords.append((element, x, y, z + stack_distance))

        # Write the stacked dimer structure to a new XYZ file
        with open(stacked_xyz_file, 'w') as f_out:
            f_out.write(f"{len(stacked_coords)}\n")
            f_out.write("Stacked dimer of helicenes\n")
            for atom in stacked_coords:
                f_out.write(f"{atom[0]} {atom[1]:.6f} {atom[2]:.6f} {atom[3]:.6f}\n")

        # Step 2: Optimize the dimer geometry using xTB
        print(f"Starting dimer geometry optimization for {stacked_xyz_file}...")
        env = {"PATH": conda_env_path + ":" + os.environ["PATH"]}
        subprocess.run(['xtb', stacked_xyz_file, '--opt', '--gfn2', '--disp'], check=True, env=env)

        # Save the optimized dimer geometry
        subprocess.run(['cp', 'xtbopt.xyz', output_file], check=True)
        print(f"Dimer optimization completed successfully. Optimized dimer saved as {output_file}.")

        return output_file

    except subprocess.CalledProcessError as e:
        print(f"Error during dimer optimization: {e}")
        return None

def extract_scc_energy_from_log(log_file):
    """
    Extracts the SCC energy from an xTB log file.

    Parameters:
    - log_file (str): The path to the log file.

    Returns:
    - float: The extracted SCC energy value.
    """
    try:
        with open(log_file, 'r') as file:
            for line in file:
                if 'SCC energy' in line:
                    return float(line.split()[3])  # Extract the SCC energy value from the line
    except (IOError, ValueError, IndexError) as e:
        print(f"Error reading SCC energy from {log_file}: {e}")
        return None


def generate_neutral_structure_from_radical(radical_smiles, radical_xyz_file="radical_structure.xyz", neutral_xyz_file="neutral_structure.xyz"):
    """
    Generates the neutral form of a molecule from its radical form using Open Babel.

    Parameters:
    - radical_smiles (str): The SMILES code of the radical form of the molecule.
    - radical_xyz_file (str): Path where the radical XYZ file will be saved.
    - neutral_xyz_file (str): Path where the neutral XYZ file will be saved.

    Returns:
    - tuple: Paths to the generated XYZ files of the radical and neutral forms.
    """
    try:
        # Step 1: Generate the radical form using Open Babel
        subprocess.run(['obabel', f'-:{radical_smiles}', '--gen3D', '-O', radical_xyz_file], check=True)
        print(f"Radical XYZ structure generated from SMILES and saved as {radical_xyz_file}.")

        # Step 2: Clean the SMILES string to generate the neutral form without explicit hydrogens
        cleaned_smiles = clean_smiles(radical_smiles)

        # Step 3: Convert the cleaned SMILES to generate the neutral structure
        subprocess.run(['obabel', f'-:{cleaned_smiles}', '--gen3D', '-O', neutral_xyz_file], check=True)
        print(f"Neutral XYZ structure generated from cleaned SMILES and saved as {neutral_xyz_file}.")

        return radical_xyz_file, neutral_xyz_file
    except subprocess.CalledProcessError as e:
        print(f"Error generating structures from SMILES: {e}")
        return None, None

def calculate_reorganization_energy(radical_smiles, neutral_smiles):
    radical_xyz = generate_xyz_from_smiles(radical_smiles, "radical_structure.xyz")
    neutral_xyz = generate_xyz_from_smiles(neutral_smiles, "neutral_structure.xyz")

    if radical_xyz and neutral_xyz:
        optimized_radical = optimize_geometry(radical_xyz, "optimized_radical.xyz")
        optimized_neutral = optimize_geometry(neutral_xyz, "optimized_neutral.xyz")

        if optimized_radical and optimized_neutral:
            E1 = calculate_single_molecule_energy(optimized_neutral)
            neutral_in_charged_state_log = run_single_point_calculation(optimized_neutral, "neutral_in_charged_state.log", charge=1, spin=2)
            E2 = extract_scc_energy_from_log(neutral_in_charged_state_log)

            E3 = calculate_single_molecule_energy(optimized_radical)
            charged_in_neutral_state_log = run_single_point_calculation(optimized_radical, "charged_in_neutral_state.log", charge=0, spin=1)
            E4 = extract_scc_energy_from_log(charged_in_neutral_state_log)

            if None not in [E1, E2, E3, E4]:
                reorganization_energy = (E2 - E1) + (E4 - E3)
                print(f"Reorganization Energy (λ): {reorganization_energy:.6f} eV")
                return reorganization_energy
            else:
                print("Failed to extract all energy values for reorganization energy calculation.")

    print("Failed to calculate reorganization energy.")
    return None

def clean_smiles(smiles):
    """
    Removes explicit hydrogens from a SMILES string for neutral molecule generation.

    Parameters:
    - smiles (str): The SMILES string to be cleaned.

    Returns:
    - str: The cleaned SMILES string without explicit hydrogens.
    """
    # Remove explicit hydrogens denoted as [H] and standard single hydrogen characters
    cleaned_smiles = smiles.replace('[H]', '').replace('H', '')
    return cleaned_smiles


def extract_orbital_energies(xtb_output_file):
    """
    Extracts the HOMO and LUMO orbital energies from an xTB output file.

    Parameters:
    - xtb_output_file (str): The path to the xTB output file.

    Returns:
    - tuple: HOMO and LUMO energy levels (in eV), or None if extraction fails.
    """
    homo_energy = None
    lumo_energy = None

    with open(xtb_output_file, 'r') as file:
        for line in file:
            # Check for the HOMO and LUMO lines using specific labels
            if '(HOMO)' in line:
                try:
                    parts = line.split()
                    homo_energy = float(parts[-2])  # Extract the HOMO energy in eV
                except (ValueError, IndexError) as e:
                    print(f"Error extracting HOMO energy from line: {line.strip()}")
            elif '(LUMO)' in line:
                try:
                    parts = line.split()
                    lumo_energy = float(parts[-2])  # Extract the LUMO energy in eV
                except (ValueError, IndexError) as e:
                    print(f"Error extracting LUMO energy from line: {line.strip()}")

    if homo_energy is None or lumo_energy is None:
        print(f"Failed to extract HOMO/LUMO energies from {xtb_output_file}. Please check the file format.")
    else:
        print(f"HOMO energy: {homo_energy:.6f} eV, LUMO energy: {lumo_energy:.6f} eV")

    return homo_energy, lumo_energy


def calculate_transfer_integral(homo_dimer, lumo_dimer, homo_monomer, lumo_monomer):
    delta_homo = abs(homo_dimer - homo_monomer)
    delta_lumo = abs(lumo_dimer - lumo_monomer)
    transfer_integral = max(delta_homo, delta_lumo) / 2
    print(f"Estimated Transfer Integral: {transfer_integral:.6f} eV")
    return transfer_integral

def run_single_point_calculation(input_xyz, output_file, charge=0, spin=1):
    """
    Runs a single-point energy calculation using xTB for a given XYZ file.
    
    Parameters:
    - input_xyz (str): The path to the input XYZ file.
    - output_file (str): The name of the output log file.
    - charge (int): The charge state of the molecule (default is 0 for neutral).
    - spin (int): The spin multiplicity of the molecule (default is 1 for a closed-shell singlet).
    
    Returns:
    - str: The path to the output log file.
    """
    try:
        # Run the xTB single-point calculation with the specified charge and spin state
        subprocess.run(['xtb', input_xyz, '--gfn2', '--sp', '--chrg', str(charge), '--uhf', str(spin - 1), '--disp'], 
                       check=True, stdout=open(output_file, 'w'))
        print(f"Single-point calculation completed for {input_xyz} with charge={charge} and spin={spin}. Results saved to {output_file}.")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"Error during single-point calculation for {input_xyz}: {e}")
        return None

def calculate_charge_carrier_mobility(reorganization_energy, transfer_integral, temperature=300):
    e = 1.602e-19
    k_B = 8.6173e-5
    hbar = 6.582e-16

    prefactor = (transfer_integral**2 / hbar**2) * (np.pi * hbar / reorganization_energy)
    exponential_term = np.exp(-reorganization_energy / (4 * k_B * temperature))
    D = prefactor * np.sqrt(np.pi / (4 * reorganization_energy * k_B * temperature)) * exponential_term
    mobility = (e / (k_B * temperature)) * D
    mobility_in_cm2 = mobility * 1e4

    print(f"Estimated charge carrier mobility: {mobility_in_cm2:.6f} cm^2/V·s")
    return mobility_in_cm2


def calculate_bipolar_mobility(hole_mobility, electron_mobility):
    """
    Calculate the overall bipolar mobility using the harmonic mean of hole and electron mobilities.

    Parameters:
    - hole_mobility (float): The mobility of the holes in cm^2/V·s.
    - electron_mobility (float): The mobility of the electrons in cm^2/V·s.

    Returns:
    - float: The overall bipolar mobility in cm^2/V·s.
    """
    if hole_mobility > 0 and electron_mobility > 0:
        bipolar_mobility = (2 * hole_mobility * electron_mobility) / (hole_mobility + electron_mobility)
        print(f"Overall Bipolar Mobility: {bipolar_mobility:.6f} cm^2/V·s")
        return bipolar_mobility
    else:
        print("One or both mobilities are zero or undefined. Cannot calculate bipolar mobility.")
        return None

def main(radical_smiles, neutral_for_hole_smiles, neutral_for_electron_smiles):
    reorganization_energy_hole = calculate_reorganization_energy(radical_smiles, neutral_for_hole_smiles)
    reorganization_energy_electron = calculate_reorganization_energy(radical_smiles, neutral_for_electron_smiles)

    if reorganization_energy_hole and reorganization_energy_electron:
        optimized_neutral = optimize_geometry("neutral_structure.xyz", "optimized_neutral.xyz")
        optimized_dimer_file = stack_and_optimize_dimer(optimized_neutral, stack_distance=3.5, output_file="final_optimized_dimer.xyz")

        monomer_output = run_single_point_calculation(optimized_neutral, "monomer_output.log", charge=0, spin=1)
        dimer_output = run_single_point_calculation(optimized_dimer_file, "dimer_output.log", charge=0, spin=1)

        homo_monomer, lumo_monomer = extract_orbital_energies(monomer_output)
        homo_dimer, lumo_dimer = extract_orbital_energies(dimer_output)

        transfer_integral = calculate_transfer_integral(homo_dimer, lumo_dimer, homo_monomer, lumo_monomer)

        hole_mobility = calculate_charge_carrier_mobility(reorganization_energy_hole, transfer_integral)
        electron_mobility = calculate_charge_carrier_mobility(reorganization_energy_electron, transfer_integral)

        print(f"Hole Mobility: {hole_mobility:.6f} cm^2/V·s")
        print(f"Electron Mobility: {electron_mobility:.6f} cm^2/V·s")

        # Calculate the overall bipolar mobility
        bipolar_mobility = calculate_bipolar_mobility(hole_mobility, electron_mobility)

if __name__ == "__main__":
    radical_smiles = '[H]C1=C([H])C2=C(C3=C([H])C(N4C5=C([N.]C(C)=N4)C=CC=C5O6)=C6C([H])=C3C([H])=C2[H])C7=C1C([H])=C([H])C8=C7C9=C([H])C([H])=C([H])C([H])=C9C([H])=C8[H]'  # Example radical SMILES
    neutral_for_hole_smiles = '[H]C1=C([H])C2=C(C3=C([H])C(N4C5=C([N]C(C)=N4)C=CC=C5O6)=C6C([H])=C3C([H])=C2[H])C7=C1C([H])=C([H])C8=C7C9=C([H])C([H])=C([H])C([H])=C9C([H])=C8[H]'  # Neutral SMILES for hole transport
    neutral_for_electron_smiles = '[H]C1=C([H])C2=C(C3=C([H])C(N4C5=C([N..]C(C)=N4)C=CC=C5O6)=C6C([H])=C3C([H])=C2[H])C7=C1C([H])=C([H])C8=C7C9=C([H])C([H])=C([H])C([H])=C9C([H])=C8[H]'  # Neutral SMILES for electron transport
    main(radical_smiles, neutral_for_hole_smiles, neutral_for_electron_smiles)

