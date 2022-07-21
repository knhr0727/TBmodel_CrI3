# TBmodel_CrI3
Tight binding model for CrI3 (monolayer, bilayer, twisted bilayer) written in Python.

How to use:
1. Create the npz files (numpy savez) of the Hamiltonian parameters of monolayers 
 (spin up and down, VLR and VRL types)

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

2. Run each codes 
1) a quick example

move to the monolayer directory `cd monolayer`.
generate a MOKE data `python monolayer_MOKE.py` (this can take a few seconds).
the data file 'DATA.mono.npz' will be generated.
run a plotting code `python plot_monolayer.py DATA.mono.npz`. 
```


