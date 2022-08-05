# TBmodel_CrI3
Tight binding model for CrI3 (monolayer, bilayer, twisted bilayer) written in Python.

How to use:
1. Create the npz files (numpy savez) of the Hamiltonian parameters of monolayers 
 (spin up and down, VLR and VRL types).  
 For each `TBband.*.py` file, the bands of single spin sepecies without spin-orbit coupling
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
The data file `DATA.mono.npz` will be generated.  
Run a figure plotting code `python plot_monolayer.py DATA.mono.npz`.  

  2.2 List of the codes
  
In `./TBmodel_CrI3/`,  
`plot_bilayer.py`: reads the MOKE and conductivity tensor output files of bilayer cases (untwisted and twisted), `DATA.*.npz`, and plots the graphs of them. 
`plot_STTB_untwisted.py`: reads the spin texture output file of the untwisted bilayer and plots the spin texture figure.  
`plot_STTB_TBL.py`: reads the spin texture output file of the twisted bilayer and plots the spin texture figure.  

In `./TBmodel_CrI3/mono_singlespin_VLR/`,  
`TBband.0.py`: creates spin-up part of the Hamiltonian parameters of monolayer (VLR type)  
`TBband.1.py`: creates spin-down part of the Hamiltonian parameters of monolayer (VLR type)  

In `./TBmodel_CrI3/mono_singlespin_VRL/`,  
`TBband.0.py`: creates spin-up part of the Hamiltonian parameters of monolayer (VRL type)  
`TBband.1.py`: creates spin-down part of the Hamiltonian parameters of monolayer (VRL type)  

In `./TBmodel_CrI3/monolayer/`,  
`monolayer_MOKE.py`: calculates the conductivity tensor and MOKE spectrum of the FM monolayer CrI3. The output file is `DATA.mono.npz`. 
By uncommenting lines 257-265, the bands will be shown.  
`plot_monolayer.py`: reads the conductivity tensor and MOKE output files of the monolayer, `DATA.mono.npz`, and plots the graphs of them.  
`monolayer_mag.py`: calculates the total magnetic moment of the FM monolayer CrI3. Outputs are printed as texts.  
`monolayer_pol.py`: calculates the electric polarization density of the FM monolayer CrI3. Outputs are printed as texts.  
`monolayer_spintexture.py`: calculates the spin texture of the FM monolayer CrI3. The output file is `spintexture_mono.npz`. 
The figure is plotted by commend `python ../plot_STTB_untwisted.py spintexture_mono.npz`  

In `./TBmodel_CrI3/bilayer/`,  
`bilayer_MOKE.py`: calculates the conductivity tensor and MOKE spectrum of the AFM AB-stacked bilayer. The output file is `DATA.biAFM.npz`. 
The figure is plotted by commend `python ../plot_bilayer.py DATA.biAFM.npz`.  
`bilayer_mag.py`: calculates the total magnetic moment of the AFM bilayer. Outputs are printed as texts.  
`bilayer_pol.py`: calculates the electric polarization density of the AFM bilayer. Outputs are printed as texts.  
`bilayer_spintexture.py`: calculates the spin texture of the AFM bilayer. The output file is `spintexture_bilayer.npz`. 
The figure is plotted by commend `python ../plot_STTB_untwisted.py spintexture_bilayer.npz`  

In `./TBmodel_CrI3/TBL_p3/`,  
`TBL_p3_MOKE.py`: calculates the conductivity tensor and MOKE spectrum of the AFM p3 twisted bilayer (TBL) CrI3. The output file is `DATA.p3.npz`. 
The figure is plotted by commend `python ../plot_bilayer.py DATA.p3.npz`.  
By modifying the integer pair (n,m) at lines 31 and 32, one can define the twist angle. The defalt is (n,m) = (2,1) corresponding to 21.79 degrees.  
`TBL_p3_mag.py`: calculates the total magnetic moment of the AFM p3 TBL CrI3. Outputs are printed as texts.  
`TBL_p3_pol.py`: calculates the electric polarization density of the AFM p3 TBL CrI3. Outputs are printed as texts.  
`TBL_p3_spintexture.py`: calculates the spin texture of the AFM p3 CrI3. The output file is `spintexture_TBL_p3.npz`. 
The figure is plotted by commend `python ../plot_STTB_TBL.py spintexture_TBL_p3.npz`  

In `./TBmodel_CrI3/TBL_p321/`,  
`TBL_p321_MOKE.py`: calculates the conductivity tensor and MOKE spectrum of the AFM p321 TBL CrI3. The output file is `DATA.p321.npz`.  
`TBL_p321_mag.py`: calculates the total magnetic moment of the AFM p321 TBL CrI3. Outputs are printed as texts.  
`TBL_p321_pol.py`: calculates the electric polarization density of the AFM p321 TBL CrI3. Outputs are printed as texts.  
`TBL_p321_spintexture.py`: calculates the spin texture of the AFM p321 TBL CrI3. The output file is `spintexture_TBL_p321.npz`.  

