import os
from rdkit import Chem
from rdkit.Chem import AllChem

# Set the input directory containing mol2 files and the charge file path
input_dir = './sampl7'
output_base_dir = './sampl7'  # Base output directory
charge_file_path = './sampl7/charge.txt'

# Read the charge file and create a dictionary mapping microstate ID to charge
charge_dict = {}
with open(charge_file_path, 'r') as f:
    for line in f:
        line = line.strip()  # Remove any leading/trailing whitespace
        if not line or line.startswith('#'):
            continue  # Skip empty lines or comment lines
        
        parts = line.split()
        if len(parts) < 2:  # Ensure there are at least two elements (ID and charge)
            print(f"Skipping malformed line: {line}")
            continue
        
        microstate_id = parts[0]
        try:
            charge = int(parts[-1])  # Charge is the last item in the line
            charge_dict[microstate_id] = charge
        except ValueError:
            print(f"Skipping line with invalid charge value: {line}")

# Loop through each file in the directory
for filename in os.listdir(input_dir):
    if filename.endswith('.mol2'):
        # Construct the full path to the file
        file_path = os.path.join(input_dir, filename)
        print(f'Processing {filename}')
        
        # Read the molecule from the mol2 file
        m = Chem.MolFromMol2File(file_path)
        if m is None:
            print(f"Could not read molecule from {file_path}. Skipping.")
            continue
        
        # Add hydrogens and generate conformers
        m2 = Chem.AddHs(m)
        confids = AllChem.EmbedMultipleConfs(m2, numConfs=100)
        print('Number of conformers: %d' % len(confids))
        
        # Filter unique conformers based on energy
        uniqueEs = []
        uniqueIDs = []
        
        for confid in confids:
            AllChem.MMFFOptimizeMolecule(m2, confId=confid)
            ff = AllChem.MMFFGetMoleculeForceField(m2, AllChem.MMFFGetMoleculeProperties(m2), confId=confid)
            E = ff.CalcEnergy()
            keep = 1
            for Eold in uniqueEs:
                if (E - Eold) ** 2 < 0.000001:
                    keep = 0
            if keep > 0:
                uniqueEs.append(E)
                uniqueIDs.append(confid)

        # Sort conformers by energy
        sortedEs = [(x, y) for x, y in sorted(zip(uniqueEs, uniqueIDs))]
        Emin = sortedEs[0][0]
        print('Lowest energy: %.4f' % Emin)

        # Extract microstate ID from filename (assuming format like 'SM25_micro000.mol2')
        microstate_id = filename.split('.')[0]  # Extract base name without extension
        if microstate_id in charge_dict:
            charge = charge_dict[microstate_id]
        else:
            print(f"Charge not found for {microstate_id}. Assuming charge 0.")
            charge = 0  # Default to 0 if charge is not found

        # Create the output directory based on the filename
        output_dir = os.path.join(output_base_dir, microstate_id)  # Create output path
        os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist

        # Write the Gaussian files for the 5 lowest energy conformers
        for i in range(min(20, len(sortedEs))):  # Limit to 5 conformers
            E, confid = sortedEs[i]
            gjf_filename = f'0{i+1}.gjf'  # Naming the file as 01.gjf, 02.gjf, etc.
            gjf_file_path = os.path.join(output_dir, gjf_filename)
            chk_filename = f'0{i+1}.chk'  # Naming the checkpoint file as 01.chk, 02.chk, etc.

            with open(gjf_file_path, 'w') as gjf_file:
                # Write header with charge info
                gjf_file.write(f"%mem=20GB\n%nprocshared=32\n%chk={chk_filename}\n")
                gjf_file.write("#P nosymm wB97XD scf(tight,fermi) integral(grid=ultrafine) opt 6-31G(d) scrf(smd,solvent=water)\n\n")
                gjf_file.write(f'Optimized Conformer with charge {charge}\n\n')
                gjf_file.write(f'{charge} 1\n')  # Write charge and multiplicity
                
                # Write atom coordinates
                for atom_id in range(m2.GetNumAtoms()):
                    atom = m2.GetAtomWithIdx(atom_id)
                    pos = m2.GetConformer(confid).GetAtomPosition(atom_id)
                    gjf_file.write(f'{atom.GetSymbol()}  {pos[0]:.8f}  {pos[1]:.8f}  {pos[2]:.8f}\n')
                
                # Write footer for the Gaussian job
                gjf_file.write("\n--link1--\n")
                gjf_file.write(f"%mem=20GB\n%nprocshared=128\n%chk={chk_filename}\n")
                gjf_file.write("#P nosymm wB97XD scf(tight,fermi) integral(grid=ultrafine) freq=noraman chkbasis geom=allcheck guess=read scrf=checkpoint\n\n")
                gjf_file.write("\n--link1--\n")
                gjf_file.write(f"%mem=20GB\n%nprocshared=128\n%chk={chk_filename}\n")
                gjf_file.write("#P nosymm wB97XD scf(tight,fermi) integral(grid=ultrafine) 6-311++G(2d,p) geom=allcheck guess=read scrf=checkpoint\n\n")

            print(f'Wrote Gaussian file for {filename}: {gjf_filename}')

print('Processing complete.')
