from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import grand


# Load PDB
pdb = PDBFile('bpti-npt.pdb')

# Create system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=12*angstrom,
                                 switchDistance=10*angstrom, constraints=HBonds)

# Make sure the LJ interactions are being switched
for f in range(system.getNumForces()):
    force = system.getForce(f)
    if 'NonbondedForce' == force.__class__.__name__:
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(1.0*nanometer)
        # Also need to switch off the long-range LJ, as this isn't used in GCMC/MD
        force.setUseDispersionCorrection(False)

# Define integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picosecond)

# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, integrator, platform)

# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298*kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Run 500 ps
simulation.step(250000)

# Write out a PDB
positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
with open('bpti-nvt2.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, positions, f)

