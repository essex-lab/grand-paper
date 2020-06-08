
import argparse
import numpy as np
from matplotlib import pyplot as plt
from simtk.unit import *


def calc_density(N, V):
    """
    Calculate the density for a given volume and number of water molecules

    Parameters
    ----------
    N : int
        Number of water molecules
    V : simtk.unit.Quantity
        Volume, in appropriate units

    Returns
    -------
    density : simtk.unit.Quantity
        Calculated density, in units of g/mL
    """
    water_mass = (2*(1.008) + 15.999) * grams / mole
    mass = N * water_mass / AVOGADRO_CONSTANT_NA
    density = mass / V
    return density.in_units_of(grams / (centimeter**3))


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', nargs='+', default=['gcmc-md.log'], help='GCMC/MD log file.')
parser.add_argument('-n', '--nwats', default=2094, type=int, help='Initial number of waters in the simulation')
parser.add_argument('-v', '--volume', default=None, type=float, help='Volume of the system (nm^3)')
parser.add_argument('-t', '--time', default=1, type=int, help='Number of timesteps per move')
parser.add_argument('-o', '--output', default='npt.csv', help='Name of the CSV file containing the output')
args = parser.parse_args()


moves = []
move_dict = {}  # Store the densities for a given number of moves (proxy for time)

# Calculate the initial density
if args.nwats is not None:
   moves.append(0)
   dens = calc_density(args.nwats, args.volume * nanometers ** 3)._value
   move_dict[0] = [dens]

# Read in the density at each point for each simulation
for filename in args.data:
    with open(filename, 'r') as f:
        for line in f.readlines():
            # Check that this is a reporting line
            if 'Current N' in line:
                for i in range(len(line.split())):
                    # Read the number of moves executed at this point
                    if line.split()[i] == 'move(s)':
                        n_moves = int(line.split()[i-1])
                        if n_moves not in moves:
                            moves.append(n_moves)
                    # Read the value of N at this point
                    elif line.split()[i] == 'Current':
                        n_wats = int(line.split()[i+3].strip('.'))
                # Calculate density
                dens = calc_density(n_wats, args.volume * nanometers ** 3)._value
                # Add this density to the appropriate list (or create a new one)
                try:
                    move_dict[n_moves].append(dens)
                except:
                    move_dict[n_moves] = [dens]

# Write data to a CSV file
with open(args.output, 'w') as f:
    f.write('Time (ns),Density (g/mL)\n')
    for i in range(len(moves)):
        # Convert this number of moves to a point in time
        time = moves[i] * args.time * 0.002 / 1000
        f.write('{}'.format(time))
        for d in move_dict[moves[i]]:
            f.write(',{}'.format(d))
        f.write('\n')

