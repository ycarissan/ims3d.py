To setup a TURBOMOLE calulation:
```
tp -g geom.xyz -n -d 'b3-lyp'
```
creates a control file.

For a periodic calulation, you need to perform a TURBOMOLE calculation with the following keywords __before the $end statement__:
```
$periodic 1
$cell
    4.6843902891
$kpoints
nkpoints 9
```
To optimize the cell add __before the $end statement__:
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
jobex -riper -c 100
```
Finally:
```
sbatch submit.job
```
