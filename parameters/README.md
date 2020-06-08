# Calculation of GCMC Parameters

This directory is used for the calculation of the excess chemical potential
and standard state volume of water.
These values are calculated in the `mu_ex` and `v_std` sub-directories.

Each repeat for the excess chemical potential is carried out as:
```commandline
cd mu_ex
python run.py -l dG-run01.log
python run.py -l dG-run02.log
```
The results can be simply extracted from the log files:
```commandline
grep 'Excess chemical potential' dG-run*log
```

The standard state volume is calculated similarly:
```commandline
cd v_std
python run.py > std_vol_01.out
python run.py > std_vol_02.out
grep 'Standard volume' std_vol*out
```
