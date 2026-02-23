Detail of updates [here](update_notes.md)

# Read before use
This program is the fruit of a scientific collaboration of the consortium [Twistar](https://twistar.home.blog/).

This is a suite of programs to generate quantum chemistry program inputs to compute IMS maps, gather results of the user performed calculations and preview the results.

Any publication resulting from the use of this program requires the citation of this article:
["Visualizing electron delocalization in contorted polycyclic aromatic hydrocarbons"](https://doi.org/10.1039/D1SC03368A)
A. Artigas, D. Hagebaum-Reignier, Y. Carissan and Y. Coquerel, Chem. Sci., 2021, DOI: 10.1039/D1SC03368A.

The user is strongly advised to read this publication before using the programs.

# Installation
## Getting the code
```
git clone https://github.com/ycarissan/ims3d.py.git
```

After this is done, you need to add the directory containing the programs to your PATH:

In bash
```
cd ims3d.py
export IMS3D_PATH=$PWD
export PATH=${IMS3D_PATH}:${PATH}
```
To make these programs permanently available you should add this line to your .bashrc file

## Requirements
### Using Conda
Install conda : [https://docs.conda.io/projects/conda/en/latest/index.html](https://docs.conda.io/projects/conda/en/latest/index.html)
Generate an environment using the conda/ims3d_conda_env.yml:
```
conda env create -f ${IMS3D_PATH}/conda/ims3d_conda_env.yml
conda activate ims3d_env
```
__IMPORTANT__: This last command activates a python environment with all necessary dependencies. This will work __only__ with the appropriate
python binary.
Once the environment is activated, the proper python3 binary is used (located at ${CONDA_INSTALLATION_DIRECTORY}/envs/ims3d_env/bin).

### Otherwise
#### Compulsory
- Python 3
- [openbabel](http://openbabel.org/wiki/Main_Page)
- [numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [pymatgen](https://pymatgen.org/)
- [Rdkit](http://rdkit.org/)
- Eugene Eeo's geode library is embedded in the distribution but can be available [here](https://github.com/eugene-eeo/spheres-from-triangles)
#### Optional
- [pyvista](https://www.pyvista.org/) (for the viewer only)

# Testing the graphical tools
All graphics in the code are done using [pyvista](https://www.pyvista.org/).
You can check if pyvista works correctly by:
```
python3 ${IMS3D_PATH}/tools/test_pyvista.py
```
If the pyvista logo appears, you are good to go with graphics.

# Testing the installation
From the `tests/` directory, run the unit tests using the conda environment's Python interpreter:
```
cd ${IMS3D_PATH}/tests
python3 -m unittest -v
```

# Usage
It is highly recommended that you reproduce the tutorial before starting your own calculations

## Tutorial
Read the tutorial [here](https://github.com/ycarissan/ims3d.py/blob/master/tutorial)

## ims3d.py

```
$ python3 ./ims3d.py -h
usage: ims3d.py [-h] [-v] [-d] [-r RADIUS] [-n NPTS] [--batch BATCH]
                [--depth DEPTH] [-o] [-i] [-p] [-a] [-c CYCLE_MAX_SIZE]
                [-f {com,dal,orca}] [--usesymmetry]
                geomfile

Generate gaussian inputs for IMS calculations.

positional arguments:
  geomfile              Geometry file in xyz format. default: geom.xyz

options:
  -h, --help            show this help message and exit
  -v, --verbose         More info
  -d, --debug           Debug info
  -r RADIUS, --radius RADIUS
                        Set the radius to 1 angstrom
  -n NPTS, --npts NPTS  Number of angular points by half circle. default: 12
  --batch BATCH, -b BATCH
                        Change the number of bq per batch. default: infinity
  --depth DEPTH         Change the depth for geodesic grid generation: 3
  -o, --orient          Reorient the molecule along its principal symmetry axis
  -i, --ignoreH         Ignore hydrogen atoms for the generation of the surface
  -p, --preview         Preview the grid and the resulting surface
  -a, --angular         Activate the deprecated angular grid
  -c CYCLE_MAX_SIZE, --cycle-max-size CYCLE_MAX_SIZE
                        Auto detect cycles of max size: 7
  -f {com,dal,orca}, --format {com,dal,orca}
                        output format: com
  --usesymmetry         Use symmetry operations (experimental)

Make sure you use the python3 interpreter that comes with the conda environment.
If this sentence makes no sense and you get error messages for missing librairies, please read the documentation.
```

## ims3d_harv.py

```
$ python3 ./ims3d_harv.py -h
usage: ims3d_harv.py [-h] [--verbose] [--angular] [--npts NPTS]
                     [-f {com,dal,orca}]
                     logfile

Harvest the calcuated data of IMS calculations.

positional arguments:
  logfile               Log filename of a series of calculations. default:
                        input_cycle_01_batch_01.log

options:
  -h, --help            show this help message and exit
  --verbose, -v         More info
  --angular, -a         Handle angular grid
  --npts NPTS, -n NPTS  Number of grid points in all 3 directions. default: 30
  -f {com,dal,orca}, --format {com,dal,orca}
                        output format: com
```

## Using molecular symmetry

When the molecule has a non-trivial point group (i.e. not C1), the `--usesymmetry`
flag reduces the number of Bq grid points to compute by exploiting symmetry.
Only symmetry-unique atoms contribute grid points; the full surface is recovered
during harvesting.

### Preparation step (`ims3d.py`)

```
python3 ims3d.py --usesymmetry geom.xyz
```

This generates:
- a reduced grid of Bq points (only around symmetry-unique atoms)
- `symmetry_operations.bin` — the symmetry operations, used later by `ims3d_harv.py`

### Harvesting step (`ims3d_harv.py`)

```
python3 ims3d_harv.py input_batch_00000.log
```

`ims3d_harv.py` automatically detects `symmetry_operations.bin` in the working
directory. If found, it expands the harvested data to the full surface by applying
all symmetry operations before writing `ims.dat`.

![Alt text](img/diff.gif)
