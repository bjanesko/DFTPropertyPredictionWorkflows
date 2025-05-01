#!/usr/bin/env python3 
# This script generates relative energy between reference microstates and with pH 0 correction 

import csv
from collections import defaultdict

def read_energy_file(filename):
    data = defaultdict(lambda: defaultdict(list))
    with open(filename, 'r') as f:
        for line in f:
            microstate_id, energy = line.strip().split()
            sm_group, micro = microstate_id.split('_')
            data[sm_group][micro].append(float(energy))
    return data

def read_smiles_file(filename):
    data = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 3:
                microstate_id, smiles, charge = row[:3]
                try:
                    charge = int(charge)
                except ValueError:
                    charge = 0  # Default to 0 if charge is not a valid integer
            else:
                microstate_id = row[0]
                smiles = row[1] if len(row) > 1 else ''
                charge = 0
            data[microstate_id] = (smiles, charge)
    return data

def calculate_relative_energy(energy, reference_energy, charge, reference_charge):
    charge_difference = charge - reference_charge
    return (energy - reference_energy) * 627.509 + charge_difference * 270.297

def process_data(energy_file, smiles_file, output_file):
    energy_data = read_energy_file(energy_file)
    smiles_data = read_smiles_file(smiles_file)

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        for sm_group, microstates in energy_data.items():
            micro000_id = f"{sm_group}_micro000"
            reference_energy = min(microstates['micro000'])
            _, reference_charge = smiles_data.get(micro000_id, ('', 0))
            
            for micro, energies in microstates.items():
                if micro == 'micro000':
                    continue
                
                microstate_id = f"{sm_group}_{micro}"
                lowest_energy = min(energies)
                smiles, charge = smiles_data.get(microstate_id, ('', 0))
                
                relative_energy = calculate_relative_energy(lowest_energy, reference_energy, charge, reference_charge)
                
                writer.writerow([
                    micro000_id,
                    microstate_id,
                    charge,
                    f"{relative_energy:.2f}",
                    "0.00",  # Placeholder for standard error
                    "0.00"   # Placeholder for uncertainty
                ])

# Usage
energy_file = 'energy3.txt'
smiles_file = 'smiles_data.csv'
output_file = 'energy3.csv'

process_data(energy_file, smiles_file, output_file)

