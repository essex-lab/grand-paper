"""
run.py
Marley Samways

This script is to run GCMC/MD using a slightly different value of the excess
chemical potential, to test its effect on the water density.
"""

import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from openmmtools.integrators import BAOABIntegrator

import grand


# Check which file to start from
parser = argparse.ArgumentParser()
parser.add_argument('-pdb', '--pdb', default='water.pdb',
                    help='PDB file to start from')
args = parser.parse_args()

# Load in a water box PDB...
pdb = PDBFile(args.pdb)

# Add ghost waters,
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology,
                                                             pdb.positions,
                                                             n=100,
                                                             pdb='water-ghosts.pdb')

# Load force field and create system
ff = ForceField('tip3p.xml')
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=12.0*angstroms,
                         switchDistance=10.0*angstroms,
                         constraints=HBonds)

# Make sure the LJ interactions are being switched
for f in range(system.getNumForces()):
    force = system.getForce(f)
    if 'NonbondedForce' == force.__class__.__name__:
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(1.0*nanometer)

# Create GCMC sampler object
gcmc_mover = grand.samplers.StandardGCMCSystemSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=298*kelvin,
                                                      excessChemicalPotential=-6.25*kilocalorie_per_mole,
                                                      standardVolume=30.345*angstroms**3,
                                                      boxVectors=np.array(pdb.topology.getPeriodicBoxVectors()),
                                                      log='density-1.log',
                                                      ghostFile='ghosts-1.txt',
                                                      rst='restart-1.rst7',
                                                      overwrite=False)

# Langevin integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picoseconds)

# Define platform
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Set up system
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Initialise the Sampler
gcmc_mover.initialise(simulation.context, ghosts)

# Run simulation - want to run 12.5M GCMC moves total, walltime may limit this, so we write checkpoints
while gcmc_mover.n_moves < 12500000:
    # Carry out 125 GCMC moves per 250 fs of MD
    simulation.step(125)
    gcmc_mover.move(simulation.context, 125)
    
    # Write data out every 0.1 ns
    if gcmc_mover.n_moves % 50000 == 0:
        gcmc_mover.report(simulation)

