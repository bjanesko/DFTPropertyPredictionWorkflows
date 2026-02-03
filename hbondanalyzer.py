import numpy as np
import math
import glob
from typing import List, Tuple, Dict

class HydrogenBondAnalyzer:
    """
    A class to analyze hydrogen bonds in molecular structures from XYZ files.
    
    Hydrogen bond criteria:
    - Distance between donor heavy atom (D) and acceptor atom (A): typically < 3.5 Å
    - Distance between hydrogen (H) and acceptor atom (A): typically < 2.5 Å  
    - Angle D-H...A: typically > 120°
    """
    
    def __init__(self, max_da_distance=3.5, max_ha_distance=2.5, min_angle=120):
        self.max_da_distance = max_da_distance  # Max distance between donor and acceptor
        self.max_ha_distance = max_ha_distance  # Max distance between hydrogen and acceptor
        self.min_angle = min_angle  # Minimum D-H...A angle in degrees
        
        # Common hydrogen bond donors and acceptors
        self.donors = {'N', 'O', 'F'}  # Atoms that can donate hydrogen
        self.acceptors = {'N', 'O', 'F'}  # Atoms that can accept hydrogen
    
    def parse_xyz_file(self, filename: str) -> Tuple[List[str], np.ndarray]:
        """
        Parse an XYZ file and return atom symbols and coordinates.
        
        Args:
            filename: Path to the XYZ file
            
        Returns:
            Tuple of (atom_symbols, coordinates_array)
        """
        atoms = []
        coords = []
        
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Skip the first two lines (number of atoms and comment)
        for line in lines[2:]:
            parts = line.strip().split()
            if len(parts) >= 4:
                symbol = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                atoms.append(symbol)
                coords.append([x, y, z])
        
        return atoms, np.array(coords)
    
    def parse_gaussian_output(self, filename: str) -> Tuple[List[str], np.ndarray]:
        """
        Parse a Gaussian output file and extract the final optimized geometry.
        
        Args:
            filename: Path to the Gaussian output file (.out or .log)
            
        Returns:
            Tuple of (atom_symbols, coordinates_array)
        """
        atoms = []
        coords = []
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the last occurrence of "Standard orientation:"
        # This contains the final optimized geometry
        last_coord_start = -1
        
        for i, line in enumerate(lines):
            if "Standard orientation:" in line or "Input orientation:" in line:
                last_coord_start = i
        
        if last_coord_start == -1:
            raise ValueError("No geometry found in Gaussian output file")
        
        # Parse the coordinate table
        # Skip 4 lines of header after "Standard orientation:"
        coord_start = last_coord_start + 5
        
        atomic_numbers = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
            9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
            16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 35: 'Br', 53: 'I'
        }
        
        for line in lines[coord_start:]:
            # Stop at the separator line
            if '-----' in line:
                break
            
            parts = line.strip().split()
            if len(parts) >= 6:
                try:
                    atomic_num = int(parts[1])
                    x = float(parts[3])
                    y = float(parts[4])
                    z = float(parts[5])
                    
                    symbol = atomic_numbers.get(atomic_num, f'X{atomic_num}')
                    atoms.append(symbol)
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        
        if not atoms:
            raise ValueError("Failed to parse coordinates from Gaussian output file")
        
        return atoms, np.array(coords)
    
    def calculate_distance(self, coord1: np.ndarray, coord2: np.ndarray) -> float:
        """Calculate Euclidean distance between two points."""
        return np.linalg.norm(coord1 - coord2)
    
    def calculate_angle(self, coord1: np.ndarray, coord2: np.ndarray, coord3: np.ndarray) -> float:
        """
        Calculate angle between three points (in degrees).
        coord2 is the vertex of the angle.
        """
        vec1 = coord1 - coord2
        vec2 = coord3 - coord2
        
        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        # Clamp to avoid numerical errors
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_rad = np.arccos(cos_angle)
        return np.degrees(angle_rad)
    
    def find_bonded_hydrogens(self, atoms: List[str], coords: np.ndarray, 
                             atom_idx: int, max_bond_length: float = 1.2) -> List[int]:
        """
        Find hydrogen atoms bonded to a given atom.
        
        Args:
            atoms: List of atom symbols
            coords: Array of atomic coordinates
            atom_idx: Index of the atom to check
            max_bond_length: Maximum bond length to consider atoms as bonded
            
        Returns:
            List of indices of bonded hydrogen atoms
        """
        bonded_hydrogens = []
        atom_coord = coords[atom_idx]
        
        for i, atom in enumerate(atoms):
            if i != atom_idx and atom == 'H':
                distance = self.calculate_distance(atom_coord, coords[i])
                if distance <= max_bond_length:
                    bonded_hydrogens.append(i)
        
        return bonded_hydrogens
    
    def identify_hydrogen_bonds(self, atoms: List[str], coords: np.ndarray) -> List[Dict]:
        """
        Identify hydrogen bonds in the molecular structure.
        
        Args:
            atoms: List of atom symbols
            coords: Array of atomic coordinates
            
        Returns:
            List of dictionaries containing hydrogen bond information
        """
        hydrogen_bonds = []
        
        # Find all potential donors (heavy atoms bonded to hydrogen)
        for i, atom in enumerate(atoms):
            if atom in self.donors:
                bonded_hydrogens = self.find_bonded_hydrogens(atoms, coords, i)
                
                # For each hydrogen bonded to this donor
                for h_idx in bonded_hydrogens:
                    h_coord = coords[h_idx]
                    donor_coord = coords[i]
                    
                    # Look for acceptor atoms
                    for j, acceptor_atom in enumerate(atoms):
                        if j != i and j != h_idx and acceptor_atom in self.acceptors:
                            acceptor_coord = coords[j]
                            
                            # Check distance criteria
                            da_distance = self.calculate_distance(donor_coord, acceptor_coord)
                            ha_distance = self.calculate_distance(h_coord, acceptor_coord)
                            
                            if (da_distance <= self.max_da_distance and 
                                ha_distance <= self.max_ha_distance):
                                
                                # Check angle criterion
                                angle = self.calculate_angle(donor_coord, h_coord, acceptor_coord)
                                
                                if angle >= self.min_angle:
                                    hbond = {
                                        'donor_idx': i,
                                        'donor_atom': atom,
                                        'hydrogen_idx': h_idx,
                                        'acceptor_idx': j,
                                        'acceptor_atom': acceptor_atom,
                                        'da_distance': round(da_distance, 3),
                                        'ha_distance': round(ha_distance, 3),
                                        'dha_angle': round(angle, 1)
                                    }
                                    hydrogen_bonds.append(hbond)
        
        return hydrogen_bonds
    
    def print_results(self, hydrogen_bonds: List[Dict], atoms: List[str]):
        """Print hydrogen bond analysis results."""
        if not hydrogen_bonds:
            print("No hydrogen bonds found with the given criteria.")
            return
        
        print(f"Found {len(hydrogen_bonds)} hydrogen bond(s):")
        print("-" * 80)
        print(f"{'#':<3} {'Donor':<15} {'H':<5} {'Acceptor':<15} {'D...A(Å)':<10} {'H...A(Å)':<10} {'D-H...A(°)':<10}")
        print("-" * 80)
        
        for i, hb in enumerate(hydrogen_bonds, 1):
            donor_info = f"{atoms[hb['donor_idx']]}{hb['donor_idx']+1}"
            h_info = f"H{hb['hydrogen_idx']+1}"
            acceptor_info = f"{atoms[hb['acceptor_idx']]}{hb['acceptor_idx']+1}"
            
            print(f"{i:<3} {donor_info:<15} {h_info:<5} {acceptor_info:<15} "
                  f"{hb['da_distance']:<10} {hb['ha_distance']:<10} {hb['dha_angle']:<10}")
    
    def save_results(self, hydrogen_bonds: List[Dict], atoms: List[str], output_file: str):
        """Save hydrogen bond analysis results to a file."""
        with open(output_file, 'w') as f:
            f.write("Hydrogen Bond Analysis Results\n")
            f.write("=" * 50 + "\n\n")
            
            if not hydrogen_bonds:
                f.write("No hydrogen bonds found with the given criteria.\n")
                return
            
            f.write(f"Found {len(hydrogen_bonds)} hydrogen bond(s):\n\n")
            f.write(f"Criteria used:\n")
            f.write(f"  - Maximum D...A distance: {self.max_da_distance} Å\n")
            f.write(f"  - Maximum H...A distance: {self.max_ha_distance} Å\n")
            f.write(f"  - Minimum D-H...A angle: {self.min_angle}°\n\n")
            
            f.write("-" * 80 + "\n")
            f.write(f"{'#':<3} {'Donor':<15} {'H':<5} {'Acceptor':<15} {'D...A(Å)':<10} {'H...A(Å)':<10} {'D-H...A(°)':<10}\n")
            f.write("-" * 80 + "\n")
            
            for i, hb in enumerate(hydrogen_bonds, 1):
                donor_info = f"{atoms[hb['donor_idx']]}{hb['donor_idx']+1}"
                h_info = f"H{hb['hydrogen_idx']+1}"
                acceptor_info = f"{atoms[hb['acceptor_idx']]}{hb['acceptor_idx']+1}"
                
                f.write(f"{i:<3} {donor_info:<15} {h_info:<5} {acceptor_info:<15} "
                       f"{hb['da_distance']:<10} {hb['ha_distance']:<10} {hb['dha_angle']:<10}\n")


def main():
    """Example usage of the HydrogenBondAnalyzer."""
    # Initialize the analyzer with custom criteria if needed
    analyzer = HydrogenBondAnalyzer(
        max_da_distance=3.5,  # Maximum donor-acceptor distance in Å
        max_ha_distance=2.5,  # Maximum hydrogen-acceptor distance in Å
        min_angle=120         # Minimum D-H...A angle in degrees
    )
    
    # Example usage
    try:
        # Choose your file type and provide the path
        # For Gaussian output files (.out or .log):
        #filename = 'molecule.out'  # or 'molecule.log'
        #atoms, coords = analyzer.parse_gaussian_output(filename)
        for filename in glob.glob("*.out"):
            atoms, coords = analyzer.parse_gaussian_output(filename)
            hbonds = analyzer.identify_hydrogen_bonds(atoms, coords)
            output_name = filename.replace('.out', '_hbonds.txt')
            analyzer.save_results(hbonds, atoms, output_name)   
        
        # For XYZ files:
        # filename = 'molecule.xyz'
        # atoms, coords = analyzer.parse_xyz_file(filename)
        
        print(f"Loaded {len(atoms)} atoms from {filename}")
        
        # Identify hydrogen bonds
        hydrogen_bonds = analyzer.identify_hydrogen_bonds(atoms, coords)
        
        # Print results
        analyzer.print_results(hydrogen_bonds, atoms)
        
        # Save results to file
        analyzer.save_results(hydrogen_bonds, atoms, 'hydrogen_bonds_results.txt')
        print(f"\nResults saved to 'hydrogen_bonds_results.txt'")
        
    except FileNotFoundError:
        print("File not found. Please provide a valid file path.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
