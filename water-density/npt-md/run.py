"""
run.py
Marley Samways

This script runs 100 ns of NPT MD on a box of pure water, in order to see the density
fluctuations observed using a Monte Carlo barostat
"""

import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

from openmmtools.integrators import BAOABIntegrator

import grand


# Load in a water box PDB...
pdb = PDBFile('../water_box-eq.pdb')

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
        # Also need to switch off the long-range LJ, as this isn't used in GCMC/MD
        force.setUseDispersionCorrection(False)

# Langevin integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picoseconds)

# Define the barostat
system.addForce(MonteCarloBarostat(1*bar, 298*kelvin, 25))

# Define platform
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Set up system
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Add a reporter - write out every 0.5 ns
simulation.reporters.append(StateDataReporter(stdout, 250000, step=True, time=True, potentialEnergy=True,
                                              temperature=True, volume=True))

# Run for 100 ns
simulation.step(50000000)

