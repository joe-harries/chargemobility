import subprocess  # Allows you to spawn new processes, connect to input/output/error pipes, and obtain return codes
import numpy as np  # Fundamental package for numerical computations
import matplotlib.pyplot as plt  # Used for plotting and visualizing results
import os  # Provides a way of interacting with the operating system, including file and environment management

def generate_xyz_from_smiles(smiles_code, output_xyz_file="generated_structure.xyz"):
    """
    Generates a 3D structure in XYZ format from a SMILES code using Open Babel.

    Parameters:
    - smiles_code (str): The SMILES code of the molecule.
    - output_xyz_file (str): The path where the generated XYZ file will be saved.

    Returns:
    - str: The path to the generated XYZ file.
    """
    try:
        # Use Open Babel to convert SMILES to a 3D optimized XYZ file
        subprocess.run(['obabel', f'-:{smiles_code}', '--gen3D', '-O', output_xyz_file], check=True)
        print(f"XYZ structure generated from SMILES and saved as {output_xyz_file}.")
        return output_xyz_file
    except subprocess.CalledProcessError as e:
        # Handle any error that occurred during the process
        print(f"Error generating XYZ file from SMILES: {e}")
        return None

def optimize_geometry(input_xyz_file, output_xyz_file, conda_env_path="/path/to/conda/env/bin"):
    """
    Optimize the geometry of a molecule using xTB. If optimization fails, it retries from the xtblast.xyz file.
    
    Parameters:
    - input_xyz_file (str): The path to the input XYZ file containing the molecule's structure.
    - output_xyz_file (str): The path where the optimized structure will be saved.
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.
    
    Returns:
    - str: The path to the optimized structure file.
    """
    xtb_path = 'xtb'  # Define the path to the xTB executable

    try:
        # Run xTB geometry optimization with full path to the input file
        print(f"Starting geometry optimization for {input_xyz_file}...")
        env = {"PATH": conda_env_path + ":" + os.environ["PATH"]}
        
        subprocess.run([xtb_path, input_xyz_file, '--opt', '--gfn2', '--disp'], check=True, env=env)

        # Save the optimized geometry to the specified output file
        subprocess.run(['cp', 'xtbopt.xyz', output_xyz_file], check=True)

        print(f"Geometry optimization completed successfully. Optimized structure saved as {output_xyz_file}.")
        return output_xyz_file

    except subprocess.CalledProcessError as e:
        print(f"Error during geometry optimization for {input_xyz_file}: {e}")

        # Check if xtblast.xyz exists to retry optimization
        if os.path.exists("xtblast.xyz"):
            print("Retrying optimization from the last saved geometry (xtblast.xyz)...")

            try:
                subprocess.run([xtb_path, 'xtblast.xyz', '--opt', '--gfn2', '--disp'], check=True, env=env)

                # Save the new optimized geometry to the output file
                subprocess.run(['cp', 'xtbopt.xyz', output_xyz_file], check=True)

                print(f"Optimization completed successfully from xtblast.xyz. Optimized structure saved as {output_xyz_file}.")
                return output_xyz_file

            except subprocess.CalledProcessError as retry_error:
                print(f"Failed to optimize from xtblast.xyz: {retry_error}")
                return None
        else:
            print("xtblast.xyz not found. Unable to retry optimization.")
            return None


def calculate_single_molecule_energy(optimized_xyz_file, conda_env_path="/path/to/conda/env/bin"):
    """
    Calculate the energy of a single optimized molecule using xTB.

    Parameters:
    - optimized_xyz_file (str): The path to the optimized XYZ file of the molecule.
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.

    Returns:
    - float: The energy of the single molecule.
    """
    try:
        # Set up the environment path and run the single-point energy calculation
        env = {"PATH": conda_env_path + ":" + os.environ["PATH"]}
        result = subprocess.run(['xtb', optimized_xyz_file, '--gfn2', '--sp', '--disp'], check=True, env=env, capture_output=True, text=True)

        # Extract the energy value from the xTB output
        energy_line = next(line for line in result.stdout.splitlines() if 'SCC energy' in line)
        single_molecule_energy = float(energy_line.split()[3])  # Extract energy value
        print(f"Single molecule SCC energy: {single_molecule_energy:.6f} Eh")
        return single_molecule_energy
    except subprocess.CalledProcessError as e:
        # Handle any error that occurred during the energy calculation
        print(f"Error calculating single molecule energy: {e}")
        return None

def change_charge_and_calculate_energy(optimized_xyz_file, charge_change, conda_env_path="/path/to/conda/env/bin"):
    """
    Modify the charge state of a molecule by adding or removing an electron and calculate the new single-point energy.

    Parameters:
    - optimized_xyz_file (str): Path to the optimized XYZ file of the neutral or radical molecule.
    - charge_change (int): The charge change (e.g., +1 for removing an electron, -1 for adding an electron).
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.

    Returns:
    - float: The single-point energy of the modified molecule.
    """
    try:
        # Step 1: Define the charge parameter for xTB
        if charge_change == -1:
            charge = "-1"  # Adding an electron results in a negative charge
        elif charge_change == +1:
            charge = "+1"  # Removing an electron results in a positive charge
        else:
            raise ValueError("charge_change must be either +1 or -1.")

        # Step 2: Run the single-point energy calculation with the modified charge
        env = {"PATH": conda_env_path + ":" + os.environ["PATH"]}
        result = subprocess.run(['xtb', optimized_xyz_file, '--gfn2', '--sp', '--chrg', charge, '--disp'], 
                                check=True, env=env, capture_output=True, text=True)
        
        # Extract energy from the output using your existing single-point energy calculation method
        energy_line = next(line for line in result.stdout.splitlines() if 'SCC energy' in line)
        single_molecule_energy = float(energy_line.split()[3])  # Extract energy value
        
        print(f"Single molecule energy with charge {charge}: {single_molecule_energy:.6f} Eh")
        return single_molecule_energy

    except subprocess.CalledProcessError as e:
        print(f"Error calculating energy with charge {charge}: {e}")
        return None
    except Exception as e:
        print(f"General error: {e}")
        return None

    
def calculate_reorganization_energy(neutral_xyz_file, radical_xyz_file, conda_env_path="/path/to/conda/env/bin"):
    """
    Calculate the reorganization energy using the four-point method.

    Parameters:
    - neutral_xyz_file (str): Path to the optimized XYZ file for the neutral molecule.
    - radical_xyz_file (str): Path to the optimized XYZ file for the radical molecule.
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.

    Returns:
    - float: The calculated reorganization energy in eV.
    """
    # Step 1: Calculate the single-point energy of the optimized neutral molecule
    E_neutral_optimized = calculate_single_molecule_energy(neutral_xyz_file, conda_env_path)
    
    # Step 2: Calculate the single-point energy of the optimized radical molecule
    E_radical_optimized = calculate_single_molecule_energy(radical_xyz_file, conda_env_path)
    
    # Step 3: Calculate the energy of the neutral molecule in the radical's geometry
    E_neutral_in_radical_geom = change_charge_and_calculate_energy(radical_xyz_file, -1, conda_env_path)
    
    # Step 4: Calculate the energy of the radical molecule in the neutral's geometry
    E_radical_in_neutral_geom = change_charge_and_calculate_energy(neutral_xyz_file, +1, conda_env_path)
    
    # Step 5: Calculate the reorganization energy using the four-point method
    reorganization_energy = (E_neutral_in_radical_geom - E_neutral_optimized) + (E_radical_in_neutral_geom - E_radical_optimized)
    
    print(f"Reorganization Energy (λ): {reorganization_energy:.6f} eV")
    return reorganization_energy


def stack_and_optimize_dimer(optimized_xyz_file, stack_distance, conda_env_path="/opt/anaconda3/envs/xtb_env/bin", output_file="optimized_dimer.xyz"):
    """
    Stacks two helicenes at a specified distance and optimizes the dimer geometry using xTB.
    
    Parameters:
    - optimized_xyz_file (str): Path to the optimized XYZ file of the single molecule.
    - stack_distance (float): Distance to stack the second molecule above the first one.
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.
    - output_file (str): Path to save the optimized dimer structure.
    
    Returns:ß
    - str: Path to the optimized dimer structure file.
    """
    try:
        # Read coordinates from the optimized XYZ file
        with open(optimized_xyz_file, 'r') as file:
            lines = file.readlines()

        num_atoms = int(lines[0].strip())
        atom_lines = lines[2:]

        # Extract atom coordinates
        coords = [line.split() for line in atom_lines]

        # Stack a second molecule by adjusting the z-coordinate
        stacked_coords = coords[:]
        for i in range(num_atoms):
            element, x, y, z = coords[i]
            new_z = float(z) + stack_distance
            stacked_coords.append([element, x, y, str(new_z)])

        # Write the stacked dimer to a new XYZ file
        stacked_xyz_file = "stacked_dimer.xyz"
        with open(stacked_xyz_file, 'w') as f:
            f.write(f"{2 * num_atoms}\n\n")  # Total atoms (twice the monomer)
            for atom in stacked_coords:
                f.write(" ".join(atom) + "\n")

        # Optimize the stacked dimer using xTB
        env = {"PATH": conda_env_path + ":" + os.environ["PATH"]}
        subprocess.run(['xtb', stacked_xyz_file, '--opt', '--gfn2', '--disp'], check=True, env=env)

        # Check for optimized output and copy to the final output file
        if not os.path.exists("xtbopt.xyz"):
            raise FileNotFoundError("Optimized dimer file (xtbopt.xyz) not found.")
        subprocess.run(['cp', 'xtbopt.xyz', output_file], check=True)

        print(f"Optimized dimer saved as {output_file}.")
        return output_file

    except Exception as e:
        print(f"Error during dimer optimization: {e}")
        return None




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
    """
    Calculate the transfer integral using the energy splitting method.
    
    Parameters:
    - homo_dimer (float): HOMO energy of the dimer in eV.
    - lumo_dimer (float): LUMO energy of the dimer in eV.
    - homo_monomer (float): HOMO energy of the monomer in eV.
    - lumo_monomer (float): LUMO energy of the monomer in eV.
    
    Returns:
    - float: Estimated transfer integral in eV.
    """
    delta_homo = abs(homo_dimer - homo_monomer)
    delta_lumo = abs(lumo_dimer - lumo_monomer)
    
    transfer_integral = max(delta_homo, delta_lumo) / 2  # Choose the larger splitting as a conservative estimate

    print(f"Estimated Transfer Integral: {transfer_integral:.6f} eV")
    return transfer_integral


def calculate_charge_carrier_mobility(reorganization_energy, transfer_integral, temperature=300):
    """
    Calculate the charge carrier mobility using Marcus theory.

    Parameters:
    - reorganization_energy (float): Reorganization energy (λ) in eV.
    - transfer_integral (float): Transfer integral (V) in eV.
    - temperature (float): Temperature in Kelvin (default is 300 K).

    Returns:
    - float: Estimated charge carrier mobility (μ) in cm^2/V·s.
    """
    # Constants
    e = 1.602e-19  # Elementary charge in coulombs
    k_B = 8.6173e-5  # Boltzmann constant in eV/K
    hbar = 6.582e-16  # Reduced Planck's constant in eV·s

    # Calculate the diffusion coefficient (D) using Marcus theory
    prefactor = (transfer_integral**2 / hbar**2) * (np.pi * hbar / reorganization_energy)
    exponential_term = np.exp(-reorganization_energy / (4 * k_B * temperature))
    D = prefactor * np.sqrt(np.pi / (4 * reorganization_energy * k_B * temperature)) * exponential_term

    # Calculate the charge carrier mobility (μ)
    mobility = (e / (k_B * temperature)) * D

    # Convert mobility to cm^2/V·s (1 m^2/V·s = 10,000 cm^2/V·s)
    mobility_in_cm2 = mobility * 1e4

    print(f"Estimated charge carrier mobility: {mobility_in_cm2:.6f} cm^2/V·s")
    return mobility_in_cm2




def test_charge_carrier_mobility_workflow(radical_smiles, neutral_smiles, stack_distance, conda_env_path="/path/to/conda/env/bin"):
    """
    Test function that runs the entire workflow: 
    1. Generates 3D structures from SMILES.
    2. Optimizes the geometries.
    3. Calculates reorganization energy using the four-point method.
    4. Stacks two radicals to estimate the transfer integral.
    5. Calculates charge carrier mobility.

    Parameters:
    - radical_smiles (str): The SMILES code of the radical molecule.
    - neutral_smiles (str): The SMILES code of the neutral molecule (radical removed).
    - stack_distance (float): Distance to stack the two radical molecules for transfer integral estimation.
    - conda_env_path (str): Path to the Conda environment's bin directory where xTB is located.
    
    Returns:
    - float: Estimated charge carrier mobility (in cm²/V·s).
    """
    
    try:
        # Step 1: Generate XYZ files from SMILES codes
        print("Generating XYZ files from SMILES...")
        radical_xyz = generate_xyz_from_smiles(radical_smiles, "radical_structure.xyz")
        neutral_xyz = generate_xyz_from_smiles(neutral_smiles, "neutral_structure.xyz")
        
        if radical_xyz is None or neutral_xyz is None:
            raise Exception("Failed to generate XYZ files from SMILES.")
        
        # Step 2: Optimize the geometries for radical and neutral molecules
        print("Optimizing geometries...")
        optimized_radical_xyz = optimize_geometry(radical_xyz, "optimized_radical.xyz", conda_env_path)
        optimized_neutral_xyz = optimize_geometry(neutral_xyz, "optimized_neutral.xyz", conda_env_path)
        
        if optimized_radical_xyz is None or optimized_neutral_xyz is None:
            raise Exception("Failed to optimize geometries.")
        
        # Step 3: Calculate the reorganization energy using the four-point method
        print("Calculating reorganization energy using the four-point method...")
        reorg_energy = calculate_reorganization_energy(optimized_neutral_xyz, optimized_radical_xyz, conda_env_path)
        
        if reorg_energy is None:
            raise Exception("Failed to calculate reorganization energy.")
        
        # Step 4: Stack two radicals and optimize the dimer for transfer integral calculation
        print("Stacking radicals and optimizing the dimer...")
        optimized_dimer_file = stack_and_optimize_dimer(optimized_radical_xyz, stack_distance, conda_env_path, output_file="optimized_dimer.xyz")
        
        if optimized_dimer_file is None:
            raise Exception("Failed to optimize the dimer.")
        
        # Step 5: Run single-point energy calculations for both dimer and monomer
        print("Running single-point energy calculations for dimer and monomer...")
        dimer_output = calculate_single_molecule_energy(optimized_dimer_file, conda_env_path)
        monomer_output = calculate_single_molecule_energy(optimized_radical_xyz, conda_env_path)
        
        if dimer_output is None or monomer_output is None:
            raise Exception("Failed to perform single-point energy calculations.")
        
        # Step 6: Extract HOMO and LUMO energies for both dimer and monomer
        print("Extracting orbital energies...")
        homo_dimer, lumo_dimer = extract_orbital_energies("dimer_output.log")
        homo_monomer, lumo_monomer = extract_orbital_energies("monomer_output.log")
        
        if homo_dimer is None or lumo_dimer is None or homo_monomer is None or lumo_monomer is None:
            raise Exception("Failed to extract orbital energies.")
        
        # Step 7: Calculate the transfer integral from the energy splittings
        print("Calculating transfer integral...")
        transfer_integral = calculate_transfer_integral(homo_dimer, lumo_dimer, homo_monomer, lumo_monomer)
        
        if transfer_integral is None:
            raise Exception("Failed to calculate transfer integral.")
        
        # Step 8: Calculate charge carrier mobility using Marcus theory
        print("Calculating charge carrier mobility...")
        mobility = calculate_charge_carrier_mobility(reorg_energy, transfer_integral)
        
        if mobility is None:
            raise Exception("Failed to calculate charge carrier mobility.")
        
        print(f"Final estimated charge carrier mobility: {mobility:.6f} cm²/V·s")
        print(reorg_energy)
        
        return mobility
    
    except Exception as e:
        print(f"Error during workflow execution: {e}")
        return None


# Example usage:

# SMILES code for the radical and neutral forms of the molecule
radical_smiles = "[H]C1=C([H])C2=C(C3=C(C4=C([H])C=C3N2COCC[Si](C)(C)C)C5=C(C([H])=C4[H])C([H])=C([H])C(C6N(C(C)(C)C(C)(C)N67=O.)[O])=C5[H])C8=C1C([H])=C([H])C9=C7([H])C([H])=C(Br)C([H])=C89" # Example radical (Nitrobenzene)
neutral_smiles = "[H]C1=C([H])C2=C(C3=C(C4=C([H])C=C3N2COCC[Si](C)(C)C)C5=C(C([H])=C4[H])C([H])=C([H])C(C6N(C(C)(C)C(C)(C)N67=O)[O])=C5[H])C8=C1C([H])=C([H])C9=C7([H])C([H])=C(Br)C([H])=C89" 

# Distance to stack two radicals for transfer integral calculation
stack_distance = 3.5  # In Ångstroms

# Run the test function to calculate charge carrier mobility
mobility = test_charge_carrier_mobility_workflow(radical_smiles, neutral_smiles, stack_distance, "/path/to/conda/env/bin")
