# Tutorial
For this tutorial to work, you need a proper xyz file (one is provided in this directory).

Furthermore, the installation must have been done properly (your PATH should contain the directory where the
python programs ims3d.py, ims3d_harv.py and ims3d_view.py are located) and the dependancies fullfiled.

## Naphtalene with default grid and 1 angström radius surface

1. Generation of the com files

```
cd tutorial
python3 ${IMS3D_PATH}/ims3d.py -r 1 naphtalene.xyz 
```

2. Gaussian calculation

2.1. Transfert of the _com_ files to the cluster (if needed)

```
scp *.com REMOTE_MACHINE:INPUT_PATH
```

2.2. Run gaussian on all _com_ files 

```
subg16 input_batch_00000.com
```
#  ☕☕☕

2.3. Transfert of the _log_ files back from the cluster (if needed)

```
scp REMOTE_MACHINE:INPUT_PATH/*.log .
```

3. Gathering of the computed IMS

```
python3 ${IMS3D_PATH}/ims3d_harv.py input_batch_00000.log
```
4. View of the results
```
python3 ${IMS3D_PATH}/ims3d_view.py
```
When a window pops u, press "q" to exit.
Viewing can be done with different rendering (run ims3d_view -h for available options)
