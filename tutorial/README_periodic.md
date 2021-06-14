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
and run 
```
jobex
```
