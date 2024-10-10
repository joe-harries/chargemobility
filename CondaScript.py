import subprocess

def optimize_geometry(input_xyz_file, output_xyz_file="optimized_6helicene.xyz"):
    """
    Optimize the geometry of a molecule using xTB.

    Parameters:
    - input_xyz_file (str): The path to the input XYZ file containing the molecule's structure.
    - output_xyz_file (str): The path where the optimized structure will be saved.

    Returns:
    - str: The path to the optimized structure file.
    """
    try:
        # Run xTB geometry optimization
        print(f"Starting geometry optimization for {input_xyz_file}...")
        subprocess.run(['xtb', input_xyz_file, '--opt', '--gfn2'], check=True)

        # Save the optimized geometry to the specified output file
        subprocess.run(['cp', 'xtbopt.xyz', output_xyz_file], check=True)

        print(f"Geometry optimization completed successfully. Optimized structure saved as {output_xyz_file}.")
        return output_xyz_file

    except subprocess.CalledProcessError as e:
        print(f"Error during geometry optimization: {e}")

        return None

# Example usage of the function
if __name__ == "__main__":
    input_file = '6helicene.xyz'
    output_file = 'optimized_6helicene.xyz'

    # Call the optimization function
    optimized_file = optimize_geometry(input_file, output_file)

    # Check the results
    if optimized_file:
        print(f"Optimized geometry saved to: {optimized_file}")
    else:
        print("Optimization failed.")
