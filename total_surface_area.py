import sys
import re
import numpy as np

def parse_gaussian_geometry(logfile):
    """Extract final geometry from Gaussian output"""
    with open(logfile, 'r') as f:
        lines = f.readlines()
    
    # Find the last Standard orientation section
    last_geom_start = -1
    for i, line in enumerate(lines):
        if 'Standard orientation:' in line or 'Input orientation:' in line:
            last_geom_start = i
    
    if last_geom_start == -1:
        raise ValueError("No geometry found in file")
    
    # Parse geometry (starts 5 lines after header)
    atoms = []
    for i in range(last_geom_start + 5, len(lines)):
        line = lines[i]
        if '-------' in line:
            break
        
        parts = line.split()
        if len(parts) >= 6:
            atom_num = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            atoms.append((atom_num, x, y, z))
    
    return atoms

def calculate_vdw_surface(atoms):
    """Calculate approximate van der Waals surface area"""
    
    # Van der Waals radii (Angstroms)
    vdw_radii = {
        1: 1.20,   # H
        6: 1.70,   # C
        7: 1.55,   # N
        8: 1.52,   # O
        9: 1.47,   # F
        15: 1.80,  # P
        16: 1.80,  # S
        17: 1.75,  # Cl
        35: 1.85,  # Br
        53: 1.98,  # I
    }
    
    atom_names = {
        1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
        15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'
    }
    
    total_area = 0
    nonpolar_area = 0
    atom_areas = []
    
    for atom_num, x, y, z in atoms:
        radius = vdw_radii.get(atom_num, 2.0)
        atom_area = 4 * np.pi * radius**2
        total_area += atom_area
        
        atom_name = atom_names.get(atom_num, 'X')
        atom_areas.append((atom_name, atom_area))
        
        # Carbon and hydrogen are nonpolar
        if atom_num in [6, 1]:  # C and H
            nonpolar_area += atom_area
    
    return {
        'total_area': total_area,
        'nonpolar_area': nonpolar_area,
        'polar_area': total_area - nonpolar_area,
        'nonpolar_percentage': (nonpolar_area/total_area)*100 if total_area > 0 else 0,
        'atom_areas': atom_areas
    }

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 surface_simple.py <gaussian_output_file>")
        sys.exit(1)
    
    gaussian_file = sys.argv[1]
    
    try:
        print(f"\nParsing {gaussian_file}...")
        atoms = parse_gaussian_geometry(gaussian_file)
        print(f"Found {len(atoms)} atoms")
        
        results = calculate_vdw_surface(atoms)
        
        print(f"\n{'='*60}")
        print(f"Van der Waals Surface Area Analysis")
        print(f"{'='*60}")
        print(f"Total Surface Area:     {results['total_area']:10.2f} Å²")
        print(f"Nonpolar Surface Area:  {results['nonpolar_area']:10.2f} Å²")
        print(f"Polar Surface Area:     {results['polar_area']:10.2f} Å²")
        print(f"Nonpolar Percentage:    {results['nonpolar_percentage']:10.1f} %")
        print(f"{'='*60}\n")
        
        # Optional: Show breakdown by atom type
        atom_counts = {}
        for atom_name, area in results['atom_areas']:
            if atom_name not in atom_counts:
                atom_counts[atom_name] = {'count': 0, 'area': 0}
            atom_counts[atom_name]['count'] += 1
            atom_counts[atom_name]['area'] += area
        
        print("Breakdown by atom type:")
        print(f"{'Atom':<6} {'Count':<8} {'Total Area (Ų)':<15}")
        print("-" * 30)
        for atom, data in sorted(atom_counts.items()):
            print(f"{atom:<6} {data['count']:<8} {data['area']:>10.2f}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
