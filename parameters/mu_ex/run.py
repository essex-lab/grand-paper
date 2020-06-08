"""
run.py
Marley Samways

Script to run a single hydration free energy calculation for water, used to
calculate the excess chemical potential of water
"""

import argparse

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import grand


parser = argparse.ArgumentParser()
parser.add_argument('-l', '--log', default='dG.log', help='Free energy log file')
args = parser.parse_args()

# Load PDB file
pdb = PDBFile('../water_box-eq.pdb')

# Set up system
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

# Run free energy calculation using grand
free_energy = grand.potential.calc_mu_ex(system=system,
                                         topology=pdb.topology,
                                         positions=pdb.positions,
                                         box_vectors=pdb.topology.getPeriodicBoxVectors(),
                                         temperature=298*kelvin,
                                         n_lambdas=30,
                                         n_samples=1000,
                                         n_equil=5000,
                                         log_file=args.log)


