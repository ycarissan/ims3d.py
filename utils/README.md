# Utils
## Usage of all.job (adapt to your needs): 
```
sbatch --array 1-20%4 all.job
```
runs the first 20 files named \*.com in the current directory
and ensures that no more than 4 jobs run simultaneously.
