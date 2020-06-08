"""
run.py
Marley Samways

Script to run a single calculation of the standard volume of water. This involves
calculating the average volume per molecule of a water box from an NPT simulation
"""


import argparse

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import grand


# Load system
pdb = PDBFile('../water_box-eq.pdb')

# Setup system
ff = ForceField("tip3p.xml")
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=12.0*angstroms,
                         constraints=HBonds,
                         switchDistance=10*angstroms)

# Make sure the LJ interactions are being switched
for f in range(system.getNumForces()):
    force = system.getForce(f)
    if 'NonbondedForce' == force.__class__.__name__:
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(1.0*nanometer)

# Define barostat - very important
system.addForce(MonteCarloBarostat(1*bar, 298*kelvin, 25))

# Calculate the volume
volume = grand.potential.calc_std_volume(system=system,
                                         topology=pdb.topology,
                                         positions=pdb.positions,
                                         box_vectors=pdb.topology.getPeriodicBoxVectors(),
                                         temperature=298*kelvin,
                                         n_samples=10000,
                                         n_equil=2500)

# Print result
print('Standard volume = {}'.format(volume.in_units_of(angstroms**3)))

