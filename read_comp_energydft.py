#!/usr/bin/env python3
# This script read energies from output files and calculate composite energies 

import os,sys
import subprocess
import re
N=1
if(len(sys.argv)>1):
    N=int(sys.argv[1])
#print('Reading %s th DFT energy '%(N))

def extract_energy(line):
    match = re.search(r'(-?\d+\.\d+)', line)
    return float(match.group(1)) if match else None

def extract_sm_number(dirname):
   match = re.match(r'SM(\d+)_micro\d+', dirname)
   return int(match.group(1)) if match else 0

# Get the current working directory
parent = os.getcwd()

results = []

# Iterate through all items in the current directory
for d in sorted(os.listdir(parent), key=extract_sm_number):
    # Check if the item is a directory
    if os.path.isdir(d):
        newd = os.path.join(parent, d)
        
        # Change to the new directory
        os.chdir(newd)
        
        # Iterate through all files ending with 'out' in the new directory
        for f in [file for file in os.listdir() if file.endswith('out')]:
            # Check if "hermal Free" exists in the file
            result = subprocess.run(['grep', '-c', 'hermal Free', f], capture_output=True, text=True)
            c = int(result.stdout.strip())

            # Check if the DFT directory and file are present 
            dftd=d.lower()
            newddft=os.path.join(newd,dftd)
            fdft='%s/%s'%(newddft,f)
            
            if c > 0 and os.path.isdir(dftd) and os.path.isfile(fdft):
                # Get the "hermal Free" line
                result = subprocess.run(['grep', 'hermal Free', f], capture_output=True, text=True)
                g = result.stdout.strip()
                
                # Get the second last "SCF Don" line (DZ energy)
                result = subprocess.run(f"tac {f} | grep -m2 'SCF Don' | tail -n1", shell=True, capture_output=True, text=True)
                e1 = result.stdout.strip()
                
                # Get the last "SCF Don" line (TZ energy)
                result = subprocess.run(f"tac {f} | grep -m1 'SCF Don'", shell=True, capture_output=True, text=True)
                e2 = result.stdout.strip()

                # Get the Nth DFT energy from the DFT single-point subdirectory , as well as the name 
                result = subprocess.run(['grep', '-c', 'SCF Don', fdft], capture_output=True, text=True)
                N0 = int(result.stdout.strip())
                if(N<=0 or N>N0):
                    raise TypeError('DFT energy request %d does not match number of DFT calcs %d\n'%(N,N0))
                result = subprocess.run(f"cat {fdft} | grep -m{N} 'SCF Don'|tail -n1", shell=True, capture_output=True, text=True)
                e3 = result.stdout.strip()
                
                # Extract energy values
                g_energy = extract_energy(g)
                e1_energy = extract_energy(e1)
                e2_energy = extract_energy(e2)
                e3_energy = extract_energy(e3)
                #print('Reading file %s with %d DFT methods energy %12.6f\n'%(fdft,N0,e3_energy))
                
                if g_energy is None or e3_energy is None or e1_energy is None:
                    print("Error: One or more energy values are None")
                    print(f"g_energy: {g_energy}, e3_energy: {e3_energy}, e1_energy: {e1_energy}")
                    # You can choose to skip this calculation, set a default value, or raise an exception
                    calculated_energy = None  # or some default value
                else:
                    calculated_energy = g_energy + e3_energy - e1_energy
                

                # Store the result
                results.append((d, calculated_energy))
        
        # Change back to the parent directory
        os.chdir(parent)

# Sort results based on directory names
results.sort(key=lambda x: extract_sm_number(x[0]))

# Print the sorted results
for d, energy in results:
    if energy is None:
        print(f"{d} None")
    else:
       #print(f"{d} composite energy: {energy:.8f}")
       print(f"{d} {energy:.8f}")
   
