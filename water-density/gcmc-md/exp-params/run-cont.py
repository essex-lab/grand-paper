"""
run-cont.py
Marley Samways

This script is to run GCMC/MD on a simulation box of pure water, sampling the entire system.
This is not how GCMC/MD would normally be run, but this is done in order to assess whether
the system will sample the correct density, where fluctuations in density arise from changes
in the number of particles as the volume is held constant.

This script is intended to continue a stopped simulation
"""

import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from openmmtools.integrators import BAOABIntegrator

import grand


def read_moves(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    n_completed = int(lines[-1].split()[4])
    n_accepted = int(lines[-1].split()[7].strip('('))

    return n_completed, n_accepted


# Check which run this is
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--run', type=int, default=2,
                    help='Which leg this represents in the full simulation')
args = parser.parse_args()

# Loading some old variables
n_moves, n_accepted = read_moves('density-{}.log'.format(args.run-1))
ghosts = grand.utils.read_ghosts_from_file('ghosts-{}.txt'.format(args.run-1))[-1]

# Load in the .pdb water box (including ghosts) to get the topology
pdb = PDBFile('water-ghosts.pdb')

# Load in the .rst7 to get the checkpointed positions and velocities
rst7 = AmberInpcrdFile('restart-{}.rst7'.format(args.run - 1))

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
                                                      excessChemicalPotential=-6.3*kilocalorie_per_mole,
                                                      standardVolume=30.0*angstroms**3,
                                                      boxVectors=np.array(pdb.topology.getPeriodicBoxVectors()),
                                                      log='density-{}.log'.format(args.run),
                                                      ghostFile='ghosts-{}.txt'.format(args.run),
                                                      rst='restart-{}.rst7'.format(args.run),
                                                      overwrite=False)

# Langevin integrator
integrator = BAOABIntegrator(298*kelvin, 1.0/picosecond, 0.002*picoseconds)

# Define platform
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Set up system
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(rst7.getPositions())  # Load positions from checkpoint
simulation.context.setVelocities(rst7.getVelocities())  # Load velocities from checkpoint
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Initialise the Sampler
gcmc_mover.initialise(simulation.context, ghosts)
# Set the number of moves to that left off at
gcmc_mover.n_moves = n_moves
gcmc_mover.n_accepted = n_accepted

# Run simulation - want to run 50M GCMC moves total, walltime may limit this, so we write checkpoints
while gcmc_mover.n_moves < 50000000:
    # Carry out 125 GCMC moves per 250 fs of MD
    simulation.step(125)
    gcmc_mover.move(simulation.context, 125)
    
    # Write data out every 0.5 ns
    if gcmc_mover.n_moves % 250000 == 0:
        gcmc_mover.report(simulation)
    

