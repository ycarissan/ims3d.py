To setup a TURBOMOLE calulation:
```
tp -g geom.xyz -d 'b3-lyp'
```
creates a control file.

For a periodic calulation, you need to perfomr a TURBOMOLE calculation with the following keywords :
```
$periodic 1
$cell
    4.6843902891
$kpoints
nkpoints 10
```
To optimize the cell add:
```
$optcell
```
Then launch:
```
get_turbofile
```
creates a submit.job file.
Edit submit.job to run
```
jobex
```
Finally:
```
sbatch submit.job
```
