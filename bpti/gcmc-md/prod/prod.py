from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
from sys import stdout
import numpy as np
import mdtraj
import grand


# Load PDB
pdb = PDBFile('../equil/bpti-uvt2.pdb')

# Add ghost waters
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology, pdb.positions,
                                                             n=15, pdb='bpti-ghosts.pdb')

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

# Define reference atoms for the GCMC sphere
ref_atoms = [{'name': 'CA', 'resname': 'TYR', 'resid': '10', 'chain': 0},
             {'name': 'CA', 'resname': 'ASN', 'resid': '43', 'chain': 0}]

# Define GCMC Sampler
gcmc_mover = grand.samplers.StandardGCMCSphereSampler(system=system,
                                                      topology=pdb.topology,
                                                      temperature=298*kelvin,
                                                      referenceAtoms=ref_atoms,
                                                      sphereRadius=4.2*angstroms,
                                                      log='bpti-prod.log',
                                                      excessChemicalPotential=-6.09*kilocalorie_per_mole,
                                                      standardVolume=30.345*angstroms**3,
                                                      ghostFile='bpti-prod.txt',
                                                      dcd='bpti-prod.dcd',
                                                      rst='bpti-prod.rst7',
                                                      overwrite=True)

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

# Prepare the GCMC sphere
gcmc_mover.initialise(simulation.context, ghosts)

# Run GCMC/MD equilibration (500k GCMC moves over 10 ns - 50 moves every 1 ps)
for i in range(10000):
    simulation.step(500)
    gcmc_mover.move(simulation.context, 50)
    # Write out a frame every 20 ps
    if (i+1) % 20 == 0:
        gcmc_mover.report(simulation)

#
# Format trajectory for visualisation
#

# Remove ghost waters from GCMC region
trj = grand.utils.shift_ghost_waters(ghost_file='bpti-prod.txt',
                                     topology='bpti-ghosts.pdb',
                                     trajectory='bpti-prod.dcd')

# Centre the trajectory on a particular residue
trj = grand.utils.recentre_traj(t=trj, resname='TYR', name='CA', resid=10)

# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output='bpti-gcmc-md.dcd')

# Write out a PDB trajectory of the GCMC sphere
grand.utils.write_sphere_traj(radius=4.2,
                              ref_atoms=ref_atoms,
                              topology='bpti-ghosts.pdb',
                              trajectory='bpti-gcmc-md.dcd',
                              output='gcmc_sphere.pdb',
                              initial_frame=True)

# Cluster water sites
grand.utils.cluster_waters(topology='bpti-ghosts.pdb',
                           trajectory='bpti-gcmc-md.dcd',
                           sphere_radius=4.2,
                           ref_atoms=ref_atoms,
                           cutoff=2.4,
                           output='watclusts.pdb')

