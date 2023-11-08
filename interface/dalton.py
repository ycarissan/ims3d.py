import pymatgen

def generate_daltonFile(geom, grid, logger, outdir="./", igrid=0, maxbq=50, all_atoms=None, atomtypes=None):
    grid_size = len(grid)
    batch=[]
    for i in range(grid_size//maxbq):
        batch.append(grid[i*maxbq:(i+1)*maxbq])
    remainder = grid_size%maxbq
    if (remainder!=0):
        batch.append(grid[len(grid)-remainder:])
    for ibatch in range(len(batch)):
        print("doing {}/{}".format(ibatch*maxbq, len(grid)))
        daltonfile = outdir + \
                "input_batch_{:05d}.dal".format(ibatch*maxbq)
        f = open(daltonfile, "w")
        try:
            fhead = open("dalton.head", "r")
            for l in fhead.readlines():
                f.write(l)
        except:
            f.write("please provide a dalton.head file\n")
        f.write("Charge=0 Atoms={} Basis=pointcharge\n".format(min(len(grid[ibatch:]),maxbq)))
        for bq in batch[ibatch]:
            f.write(
                "Bq     {0[0]:16.10f} {0[1]:16.10f} {0[2]:16.10f}\n".format(bq))
        try:
            ffoot = open("dalton.foot", "r")
            for l in ffoot.readlines():
                f.write(l)
        except:
            f.write("please provide a dalton.foot file\n")
        f.close()
    return

