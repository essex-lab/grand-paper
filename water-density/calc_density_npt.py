
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
parser.add_argument('-d', '--data', default=['npt.out'], nargs='+', help='Stdout output from NPT MD simulation')
parser.add_argument('-n', '--nwats', default=2094, type=int, help='Number of waters in the simulation')
parser.add_argument('-v', '--volume', default=None, type=float, help='Initial volume of the system (nm^3)')
parser.add_argument('-o', '--output', default='npt.csv', help='Name of the CSV file containing the output')
args = parser.parse_args()


times = []
densities = []

# Calculate initial density
if args.volume is not None:
   times.append(0.0)
   densities.append([calc_density(args.nwats, args.volume*nanometers**3)._value])

# Calculate the density for each frame of each simulation
for filename in args.data:
    with open(filename, 'r') as f:
        # Columns to read the time and volume from
        n_cols = False
        time_col = False
        volume_col = False
    
        for line in f.readlines():
            # Read the header to figure out which columns contain the time and volume
            if line.startswith('#'):
                n_cols = len(line.split(','))
                for i in range(n_cols):
                    if 'Time' in line.split(',')[i]:
                        time_col = i
                    elif 'Box Volume' in line.split(',')[i]:
                        volume_col = i
            elif n_cols is False:
                continue
            # Read data
            elif len(line.split(',')) == n_cols:
                # Read time and volume
                time = float(line.split(',')[time_col]) / 1000  # Convert from ps to ns
                volume = float(line.split(',')[volume_col]) * nanometers ** 3
                # Calculate density (remove units)
                density = calc_density(args.nwats, volume)._value
                # Check if a value of the density has already been calculated for this time value
                # If not, create a new list for this point in time
                if not any([np.isclose(time, t) for t in times]):
                    times.append(time)
                    densities.append([density])
                else:
                    for i in range(len(times)):
                        if np.isclose(time, times[i]):
                            densities[i].append(density)

# Write data out to a CSV file
with open(args.output, 'w') as f:
    f.write('Time (ns),Density (g/mL)\n')
    for i in range(len(times)):
        f.write('{}'.format(times[i]))
        for d in densities[i]:
            f.write(',{}'.format(d))
        f.write('\n')

