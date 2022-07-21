# TBmodel_CrI3
Tight binding model for CrI3 (monolayer, bilayer, twisted bilayer) written in Python.

How to use:
1. Create the npz files (numpy savez) of the Hamiltonian parameters of monolayers 
 (spin up and down, VLR and VRL types)  
 For each `TBband.*.py` file, bands of single spin sepecies without spin-orbit coupling
 corresponds to each file will be shown. Just close the figure to proceed.  

```
cd mono_singlespin_VRL
python TBband.0.py 
python TBband.1.py
cd ..
cd mono_singlespin_VLR
python TBband.0.py 
python TBband.1.py
cd ..
```

2. Run each code

  2.1 A quick example
  
Move to the monolayer directory `cd monolayer`.  
Generate a MOKE data `python monolayer_MOKE.py` (this can take a few seconds).  
The data file 'DATA.mono.npz' will be generated.  
Run a plotting code `python plot_monolayer.py DATA.mono.npz`.  
