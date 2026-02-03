#!/usr/bin/env python3
import os
import glob
import numpy as np
import csv
from typing import Dict, List, Tuple

class TitrationCalculator:
    """
    Calculate titration curves and LogD values from Gaussian output files.
    """
    
    def __init__(self):
        self.a2kJ = 2625.5  # Hartree to kJ/mol conversion
        self.RT = 298.15 * 0.008314  # RT in kJ/mol at T=298 K
        self.H_plus = -0.48  # Proton energy in Hartree
        
    def get_last_line_with_pattern(self, filename: str, pattern: str, occurrence: int = 1) -> str:
        """Get the last occurrence of a pattern in a file."""
        matches = []
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if pattern in line:
                        matches.append(line.strip())
        except FileNotFoundError:
            return ""
        
        if len(matches) >= occurrence:
            return matches[-occurrence]
        return ""
    
    def extract_scf_energy(self, line: str) -> float:
        """Extract SCF energy value from 'SCF Done' line."""
        if not line:
            return 0.0
        parts = line.split()
        try:
            if '=' in parts:
                idx = parts.index('=')
                if idx + 1 < len(parts):
                    return float(parts[idx + 1])
        except (ValueError, IndexError):
            pass
        return 0.0
    
    def extract_free_energy(self, filename: str) -> float:
        """Extract the lowest thermal free energy correction."""
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            free_energies = []
            for line in lines:
                if "hermal Free" in line:
                    parts = line.split()
                    if parts:
                        try:
                            free_energies.append(float(parts[-1]))
                        except ValueError:
                            continue
            
            if free_energies:
                return min(free_energies)
        except FileNotFoundError:
            pass
        
        return 0.0
    
    def getG(self, filename: str) -> float:
        """
        Calculate Gibbs free energy from Gaussian output file.
        G = E(large basis) + G0(thermal correction) - E0(small basis)
        """
        E_line = self.get_last_line_with_pattern(filename, "SCF Done", 1)
        E = self.extract_scf_energy(E_line)
        
        E0_line = self.get_last_line_with_pattern(filename, "SCF Done", 2)
        E0 = self.extract_scf_energy(E0_line)
        
        G0 = self.extract_free_energy(filename)
        
        G = E + G0 - E0
        
        return G
    
    def find_lowest_energy(self, directory: str) -> Tuple[float, str]:
        """Find the lowest energy conformer in a directory."""
        pattern = os.path.join(directory, '*.out')
        files = glob.glob(pattern)
        
        if not files:
            return 0.0, ""
        
        lowest_G = None
        lowest_file = ""
        
        for f in files:
            if '/sp/' in f or '\\sp\\' in f:
                continue
            
            G = self.getG(f)
            
            if lowest_G is None or G < lowest_G:
                lowest_G = G
                lowest_file = f
        
        return lowest_G if lowest_G is not None else 0.0, lowest_file
    
    def collect_energies(self, base_dirs: Dict[str, str]) -> Dict[str, float]:
        """
        Collect lowest energies from directories.
        
        Args:
            base_dirs: Dictionary with keys like 'neutral_water', 'neutral_octanol', etc.
                      and values as directory paths
        
        Returns:
            Dictionary of energies
        """
        energies = {}
        
        for key, directory in base_dirs.items():
            if os.path.exists(directory):
                G, filename = self.find_lowest_energy(directory)
                energies[key] = G
                print(f"{key}: G = {G:.6f} Hartree from {filename}")
            else:
                print(f"Warning: Directory {directory} not found")
                energies[key] = 0.0
        
        return energies
    
    def calculate_Ka_values(self, Gwat0: float, Gwat1: float, Gwat2: float) -> Dict[str, float]:
        """
        Calculate acid dissociation constants.
        """
        # Energy differences in kJ/mol
        DGw21 = self.a2kJ * (Gwat2 - Gwat1 - self.H_plus)
        DGw10 = self.a2kJ * (Gwat1 - Gwat0 - self.H_plus)
        
        # Ka values: Ka = exp(-DG/(RT)) 
        Ka1 = np.exp(-DGw21 / self.RT)
        Ka2 = np.exp(-DGw10 / self.RT)
        
        return {
            'DGw21': DGw21,
            'DGw10': DGw10,
            'Ka1': Ka1,
            'Ka2': Ka2
        }
    
    def calculate_populations_water(self, pH_values: List[float], Ka1: float, Ka2: float) -> np.ndarray:
        """
        Calculate populations of each species in water as a function of pH.
        
        Returns array with columns: [pH, cH, cA0, cA1, cA2]
        """
        results = []
        
        for pH in pH_values:
            cH = 10**(-pH)
            
            # Denominator for all species
            denom = 1 + Ka1/cH + Ka1*Ka2/(cH**2)
            
            # Populations
            cA2 = 1 / denom  # Doubly protonated (most acidic)
            cA1 = (Ka1/cH) / denom  # Singly protonated
            cA0 = (Ka1*Ka2/(cH**2)) / denom  # Neutral
            
            results.append([pH, cH, cA0, cA1, cA2])
        
        return np.array(results)
    
    def calculate_partition_coefficients(self, energies: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate partition coefficients for each charge state.
        
        DGow = Goct - Gwat (in kJ/mol)
        """
        Gwat0 = energies.get('neutral_water', 0.0)
        Gwat1 = energies.get('plusone_water', 0.0)
        Gwat2 = energies.get('twoplus_water', 0.0)
        
        Goct0 = energies.get('neutral_octanol', 0.0)
        Goct1 = energies.get('plusone_octanol', 0.0)
        Goct2 = energies.get('twoplus_octanol', 0.0)
        
        DGow00 = self.a2kJ * (Goct0 - Gwat0)
        DGow11 = self.a2kJ * (Goct1 - Gwat1)
        DGow22 = self.a2kJ * (Goct2 - Gwat2)
        
        # Partition coefficients: P = exp(-DGow/RT)
        P0 = np.exp(-DGow00 / self.RT)  
        P1 = np.exp(-DGow11 / self.RT)
        P2 = np.exp(-DGow22 / self.RT)
        
        return {
            'DGow00': DGow00,
            'DGow11': DGow11,
            'DGow22': DGow22,
            'P0': P0,
            'P1': P1,
            'P2': P2
        }
    
    def calculate_logD(self, pH_values: List[float], Ka1: float, Ka2: float, 
                      P0: float, P1: float, P2: float) -> np.ndarray:
        """
        Calculate distribution coefficient (LogD) as a function of pH.
        
        LogD = log10((cA0_oct + cA1_oct + cA2_oct) / (cA0_wat + cA1_wat + cA2_wat))
        """
        results = []
        
        for pH in pH_values:
            cH = 10**(-pH)
            
            # Populations in water
            denom_wat = 1 + Ka1/cH + Ka1*Ka2/(cH**2)
            cA2_wat = 1 / denom_wat
            cA1_wat = (Ka1/cH) / denom_wat
            cA0_wat = (Ka1*Ka2/(cH**2)) / denom_wat
            
            # Populations in octanol (using partition coefficients)
            cA2_oct = cA2_wat * P2
            cA1_oct = cA1_wat * P1
            cA0_oct = cA0_wat * P0
            
            # Distribution coefficient
            total_oct = cA0_oct + cA1_oct + cA2_oct
            total_wat = cA0_wat + cA1_wat + cA2_wat
            
            if total_wat > 0:
                D = total_oct / total_wat
                logD = np.log10(D) if D > 0 else -10
            else:
                logD = -10
            
            results.append([pH, cH, cA0_wat, cA1_wat, cA2_wat, 
                           cA0_oct, cA1_oct, cA2_oct, total_wat, logD])
        
        return np.array(results)
    
    def save_results(self, output_file: str, energies: Dict[str, float], 
                    Ka_values: Dict[str, float], partition_data: Dict[str, float],
                    populations: np.ndarray, logD_data: np.ndarray):
        """Save all results to CSV file."""
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow(['Predicted titration curve and LogD'])
            writer.writerow([])
            
            # Step 1: Energies
            writer.writerow(['Step 1: Computed Gibbs free energies (Hartree)'])
            writer.writerow([])
            writer.writerow(['Gwat0', energies.get('neutral_water', 0.0)])
            writer.writerow(['Gwat1', energies.get('plusone_water', 0.0)])
            writer.writerow(['Gwat2', energies.get('twoplus_water', 0.0)])
            writer.writerow(['Goct0', energies.get('neutral_octanol', 0.0)])
            writer.writerow(['Goct1', energies.get('plusone_octanol', 0.0)])
            writer.writerow(['Goct2', energies.get('twoplus_octanol', 0.0)])
            writer.writerow([])
            writer.writerow(['H(+)', self.H_plus])
            writer.writerow([])
            
            # Step 2: Ka values and partition coefficients
            writer.writerow(['Step 2: Equilibrium constants and partition coefficients'])
            writer.writerow([])
            writer.writerow(['', 'Energy (kJ/mol)', 'Ka value'])
            writer.writerow(['DGw21', Ka_values['DGw21'], Ka_values['Ka1']])
            writer.writerow(['DGw10', Ka_values['DGw10'], Ka_values['Ka2']])
            writer.writerow([])
            writer.writerow(['', 'Energy (kJ/mol)', 'Partition Coeff'])
            writer.writerow(['DGow00', partition_data['DGow00'], partition_data['P0']])
            writer.writerow(['DGow11', partition_data['DGow11'], partition_data['P1']])
            writer.writerow(['DGow22', partition_data['DGow22'], partition_data['P2']])
            writer.writerow([])
            
            # Step 3: Populations and LogD
            writer.writerow(['Step 3: Populations in water as a function of pH'])
            writer.writerow([])
            writer.writerow(['pH', 'cH', 'cAwat (neutral)', 'cHAwat (+1)', 'cH2Awat (+2)'])
            for row in populations:
                writer.writerow([f"{row[0]:.1f}", f"{row[1]:.2e}", 
                               f"{row[2]:.2e}", f"{row[3]:.2e}", f"{row[4]:.2e}"])
            writer.writerow([])
            
            # Step 4: LogD values
            writer.writerow(['Step 4: Distribution coefficient (LogD) as a function of pH'])
            writer.writerow([])
            writer.writerow(['pH', 'cH', 'cAwat', 'cHAwat', 'cH2Awat', 
                           'cAoct', 'cHAoct', 'cH2Aoct', 'Total', 'LogD'])
            for row in logD_data:
                writer.writerow([f"{row[0]:.1f}", f"{row[1]:.2e}",
                               f"{row[2]:.2e}", f"{row[3]:.2e}", f"{row[4]:.2e}",
                               f"{row[5]:.2e}", f"{row[6]:.2e}", f"{row[7]:.2e}",
                               f"{row[8]:.2e}", f"{row[9]:.4f}"])
        
        print(f"\nResults saved to: {output_file}")


def main():
    """
    Main function to run titration curve calculations.
    """
    import sys
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate titration curves and LogD values from Gaussian output files')
    parser.add_argument('base_path', nargs='?', default='.', 
                       help='Base directory containing neutral/plusOne/twoplus subdirectories (default: current directory)')
    parser.add_argument('--ph', nargs='+', type=float,
                       help='Specific pH values to calculate LogD (e.g., --ph 2.0 5.5 7.4 9.0)')
    parser.add_argument('--ph-range', nargs=3, type=float, metavar=('START', 'STOP', 'STEP'),
                       help='pH range to calculate (start, stop, step) (e.g., --ph-range 1.0 12.0 0.5)')
    parser.add_argument('--no-csv', action='store_true',
                       help='Do not save results to CSV file')
    
    args = parser.parse_args()
    
    calculator = TitrationCalculator()
    
    # Define directory structure
    base_path = args.base_path
    base_dirs = {
        'neutral_water': f'{base_path}/neutral/water',
        'neutral_octanol': f'{base_path}/neutral/octanol',
        'plusone_water': f'{base_path}/plusOne/water',
        'plusone_octanol': f'{base_path}/plusOne/octanol',
        'twoplus_water': f'{base_path}/twoplus/water',
        'twoplus_octanol': f'{base_path}/twoplus/octanol'
    }
    
    print("="*80)
    print("Titration Curve and LogD Calculator")
    print("="*80)
    print("\nStep 1: Collecting energies from Gaussian output files")
    print("-"*80)
    
    # Collect energies
    energies = calculator.collect_energies(base_dirs)
    
    # Calculate Ka values
    print("\nStep 2: Calculating equilibrium constants")
    print("-"*80)
    Ka_values = calculator.calculate_Ka_values(
        energies.get('neutral_water', 0.0),
        energies.get('plusone_water', 0.0),
        energies.get('twoplus_water', 0.0)
    )
    
    print(f"DGw21 = {Ka_values['DGw21']:.2f} kJ/mol, Ka1 = {Ka_values['Ka1']:.2e}")
    print(f"DGw10 = {Ka_values['DGw10']:.2f} kJ/mol, Ka2 = {Ka_values['Ka2']:.2e}")
    
    # Calculate partition coefficients
    print("\nStep 3: Calculating partition coefficients")
    print("-"*80)
    partition_data = calculator.calculate_partition_coefficients(energies)
    
    print(f"DGow00 = {partition_data['DGow00']:.2f} kJ/mol, P0 = {partition_data['P0']:.2e}")
    print(f"DGow11 = {partition_data['DGow11']:.2f} kJ/mol, P1 = {partition_data['P1']:.2e}")
    print(f"DGow22 = {partition_data['DGow22']:.2f} kJ/mol, P2 = {partition_data['P2']:.2e}")
    
    # Determine pH values to calculate
    if args.ph:
        # Use specific pH values provided
        pH_values = sorted(args.ph)
        print(f"\nUsing specified pH values: {pH_values}")
    elif args.ph_range:
        # Use pH range
        start, stop, step = args.ph_range
        pH_values = np.arange(start, stop + step/2, step).tolist()
        print(f"\nUsing pH range: {start} to {stop} with step {step}")
    else:
        # Default: pH 1.0 to 12.5 with step 0.5
        pH_values = np.arange(1.0, 12.51, 0.5).tolist()
    
    # Calculate populations in water
    print("\nStep 4: Calculating species populations")
    print("-"*80)
    populations = calculator.calculate_populations_water(
        pH_values, Ka_values['Ka1'], Ka_values['Ka2']
    )
    
    # Calculate LogD
    print("\nStep 5: Calculating LogD values")
    print("-"*80)
    logD_data = calculator.calculate_logD(
        pH_values,
        Ka_values['Ka1'], Ka_values['Ka2'],
        partition_data['P0'], partition_data['P1'], partition_data['P2']
    )
    
    # Print LogD values on screen
    print("\n" + "="*80)
    print("LogD Values at Different pH")
    print("="*80)
    print(f"{'pH':<8} {'LogD':<12} {'cAwat':<12} {'cHAwat':<12} {'cH2Awat':<12}")
    print("-"*80)
    for row in logD_data:
        pH, cH, cA0_wat, cA1_wat, cA2_wat, cA0_oct, cA1_oct, cA2_oct, total, logD = row
        print(f"{pH:<8.2f} {logD:<12.4f} {cA0_wat:<12.2e} {cA1_wat:<12.2e} {cA2_wat:<12.2e}")
    print("="*80)
    
    # Save results to CSV unless --no-csv flag is set
    if not args.no_csv:
        output_file = 'titration_curve_results.csv'
        calculator.save_results(
            output_file, energies, Ka_values, partition_data,
            populations, logD_data
        )
    
    print("\n" + "="*80)
    print("Calculation complete!")
    print("="*80)


if __name__ == "__main__":
    main()
