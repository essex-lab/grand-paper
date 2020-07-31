# Water Density Analysis

This directory contains the scripts needed to carry out the density analysis
of water using both GCMC/MD and NPT, as reported in the Supporting Information.

## Directory Structure

The directories here are arranged as follows:

- `npt-md` : Directory for shorter NPT simulations (25 ns, reporting every 0.1 ns)
- `gcmc-md` : Directory for the shorter GCMC/MD simulations
    - `equil` : GCMC/MD simulations starting from densities which are approximately 20 %
    higher or lower than the equilibrium density.
        - `low` : Simulations starting from a density around 20 % too low.
        - `high` : Simulations starting from a density around 20 % too high.
    - `sensitivity` : GCMC/MD simulations using different values for the excess
    chemical potential to assess the sensitivity of the results to this parameter.
        - `mu_-6.15` : Excess chemical potential of -6.15 kcal/mol.
        - `mu_-6.20` : Excess chemical potential of -6.20 kcal/mol.
        - `mu_-6.25` : Excess chemical potential of -6.25 kcal/mol.

## Simulation

### NPT MD

As before, for each repeat, simply run:
```commandline
cd npt-md
python run.py > npt-01.out
```

### GCMC/MD

To run GCMC/MD, starting from a low density, we run
```commandline
cd gcmc-md/equil/low/
mkdir run1
cd run1
python ../run.py -pdb ../../../../../water_box-eq.pdb
```
Here, 20 % of the water molecules are selected at random, and 
deleted prior to the simulation starting.
The simulations can be restarted as described previously.

Similarly, starting from a high density, we can run the following:
```commandline
cd gcmc-md/equil/high/
mkdir run1
cd run1
python ../run.py -pdb ../water_box-high.pdb
```
where `water_box-high.pdb` has been prepared with a density which is
around 20 % higher than the equilibrium value.

For the sensitivity analysis, we run simulations starting from a density 
of 1.004 g/mL, with excess chemical potential values of -6.15, -6.20 and -6.25
kcal/mol, as follows:
```commandline
cd gcmc-md/sensitivity/mu_-6.15
mkdir run1
cd run1
python ../run.py -pdb ../../../../../water_box-eq.pdb
```
This can then be repeated for directories `mu_-6.20` and `mu_-6.25`.

## Analysis

As before, the densities can be calculated as:
```commandline
python ../calc_density_npt.py -d npt-md/*.out -n 2094 -v 62.391 -o npt-25ns.csv
python ../calc_density_uvt.py -d gcmc-md/equil/low/run?/density-?.log -n 1676 -v 62.391 -o uvt-low.csv
python ../calc_density_uvt.py -d gcmc-md/equil/high/run?/density-?.log -n 2500 -v 62.391 -o uvt-high.csv
python ../calc_density_uvt.py -d gcmc-md/sensitivity/mu_-6.15/run?/density-?.log -n 2094 -v 62.391 -o uvt-mu_-6.15.csv
python ../calc_density_uvt.py -d gcmc-md/sensitivity/mu_-6.20/run?/density-?.log -n 2094 -v 62.391 -o uvt-mu_-6.20.csv
python ../calc_density_uvt.py -d gcmc-md/sensitivity/mu_-6.25/run?/density-?.log -n 2094 -v 62.391 -o uvt-mu_-6.25.csv
```

And these data can be plotted as:
```commandline
python ../plot_densities.py --npt npt-25ns.csv --uvt uvt-low.csv --xlim 25.5 -o density-low.pdf
python ../plot_densities.py --npt npt-25ns.csv --uvt uvt-high.csv --xlim 25.5 -o density-high.pdf
python ../plot_densities.py --npt npt-25ns.csv --uvt uvt-mu_-6.15.csv --xlim 25.5 -o density-mu_-6_15.pdf
python ../plot_densities.py --npt npt-25ns.csv --uvt uvt-mu_-6.20.csv --xlim 25.5 -o density-mu_-6_20.pdf
python ../plot_densities.py --npt npt-25ns.csv --uvt uvt-mu_-6.25.csv --xlim 25.5 -o density-mu_-6_25.pdf
```

