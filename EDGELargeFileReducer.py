# This code is part of my Final Year Project and is used to reduce the filesize of the HiRes simulation data
# This is due to computational limitations, and to allow for full analysis to be completed on the HiRes galaxy simulations

import numpy as np

def process_large_file(input_file):
    radii = []
    masses = []

    with open(input_file, 'r') as file:
        next(file)  # Skip the first row
        next(file)  # Skip the second row
        for line in file:
            values = line.strip().split(',')
            # Extract x, y, z positions (columns 1, 2, 3) and mass (column 7)
            x, y, z = float(values[0]), float(values[1]), float(values[2])
            mass = float(values[6])
            radius = np.sqrt(x**2 + y**2 + z**2)
            radii.append(radius)
            masses.append(mass)

            # Delete x, y, z to free up memory
            del x, y, z

    return radii, masses

galaxies = ['605', '624', '1445', '1459']
file_types = ['star', 'gas', 'dm']

for galaxy in galaxies:
    for file_type in file_types:
        input_file = f'{galaxy}/{galaxy}_fiducial_hires_{file_type}.data'
        output_file = f'{galaxy}/{galaxy}_{file_type}_mass_HIRES.txt'

        sorted_radii, cumulative_mass = process_large_file(input_file)
        
        # Write radii and associated masses to file
        with open(output_file, 'w') as outfile:
            for radius, mass in zip(sorted_radii, cumulative_mass):
                outfile.write(f"{radius} {mass}\n")
                 
        print(f"File {output_file} has been created.")
