# Simulation of Bovine Pancreatic Trypsin Inhibitor

This directory contains the scripts and files necessary to carry out the
simulations GCMC/MD and canonical MD simulations reported for bovine pancreatic
trypsin inhibitor (BPTI).

The GCMC/MD simulations can be found in the `gcmc-md/` sub-directory, and the NVT
simulations in `nvt-md`.
Each simulation involves three equilibration stages, followed by a production run.
For GCMC/MD, this is:
```commandline
cd gcmc-md/equil/
python equil-uvt1.py
python equil-npt.py
python equil-uvt2.py
cd ../prod/
python prod.py
```
and for NVT MD:
```commandline
cd nvt-md/equil/
python equil-nvt1.py
python equil-npt.py
python equil-nvt2.py
cd ../prod/
python prod.py
```

The trajectory files for visualisation are `gcmc-md/prod/bpti-gcmc-md.dcd` and 
`nvt-md/prod/bpti-nvt.dcd`.
These can optionally be visualised using a PyMOL script included in _grand_,
for example:
```commandline
cd gcmc-md/prod
pymol ~/software/grand/grand/scripts/gcmc_pymol.py -- -top bpti-ghosts.pdb -trj bpti-gcmc-md.dcd -s gcmc_sphere.pdb
```

The clustered water locations are saved as `gcmc-md/prod/watclusts.pdb` and 
`nvt-md/prod/watclusts.pdb`.
In these PDB files, the occupancy and temperature factor columns both contain
the fraction of frames in which the cluster is present.
