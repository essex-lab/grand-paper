# Water Density Analysis

This directory contains the scripts needed to carry out the density analysis
of water using both GCMC/MD and NPT, as reported.
This directory contains the PDB file of the initial water structure 
(`water_box-eq2.pdb`), with a density of 1.004 g/mL.
The structure used for the simulation shown in the Supporting Information is
`water_box-eq.pdb` - the density of this structure is 0.978 g/mL

## Directory Structure

The directories here are arranged as follows:

- `npt-md` : Directory for the NPT simulations
- `gcmc-md` : Directory for the GCMC/MD simulations
    - `sim-params` : GCMC/MD simulations using the calculated values for the excess
    chemical potential and standard state volume (-6.09 kcal/mol and 30.345 A^3,
    respectively)
    - `exp-params` : GCMC/MD simulations using the experimental values for the excess
    chemical potential and standard state volume (-6.30 kcal/mol and 30.000 A^3,
    respectively)

## Simulation

### NPT MD

For each repeat, simply run:
```commandline
cd npt-md
python run.py > npt-01.out
```

### GCMC/MD

To run GCMC/MD, using the calculated parameters, run the following:
```commandline
cd gcmc-md/sim-params
mkdir run1
cd run1
python ../run.py -pdb ../../../water_box-eq2.pdb
```
The same can also be done with `../../../water_box-eq.pdb`, as presented in the SI.
This simulation is very slow (given the unusually high frequency of GCMC
moves, owing to the large system size), so the simulation may need to be 
restarted several times:
```commandline
python ../run-cont.py -r 2
python ../run-cont.py -r 3
```
and so on.
The same is also carried out for sub-directories `run2` and `run3`.

Similarly, for the simulation using the experimental parameters:
```commandline
cd gcmc-md/exp-params
mkdir run1
cd run1
python ../run.py -pdb ../../../water_box-eq2.pdb
python ../run-cont.py -r 2
python ../run-cont.py -r 3
```

## Analysis

The density values at each point in time for the NPT simulations are calculated
as:
```commandline
python calc_density_npt.py -d npt-md/*.out -n 2094 -v 62.391 -o npt-densities.csv
```
where the `-n` flag indicates the number of waters in the system, and the `-v`
flag indicates the initial volume (cubic nanometres).
Similarly, the densities can be calculated for each of the GCMC/MD simulations:
```commandline
python calc_density_uvt.py -d gcmc-md/sim-params/run?/density-?.log -n 2094 -v 62.391 -o uvt-densities-new.csv
python calc_density_uvt.py -d gcmc-md/exp-params/run?/density-?.log -n 2094 -v 62.391 -o uvt-densities-old.csv
```
where `-v` gives the volume of the simulation, and the `-n` flag takes the
initial number of water molecules.

The graphical comparisons can then be plotted using the commands below.
```commandline
python plot_densities.py --npt npt-densities.csv --uvt uvt-densities-new.csv --xlim 100.5 -o density-new.pdf
python plot_densities.py --npt npt-densities.csv --uvt uvt-densities-old.csv --xlim 100.5 -o density-old.pdf
```
