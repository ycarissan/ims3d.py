import pymatgen

def generate_daltonFile(geom, grid, logger, outdir="./", igrid=0, maxbq=50, all_atoms=None, atomtypes=None):
    print("doing {}/{}".format(igrid, len(grid)))
    daltonfile = outdir + \
        "input_batch_{:05d}.dal".format(igrid)
    f = open(daltonfile, "w")

    try:
        fhead = open("dalton.head", "r")
        for l in fhead.readlines():
            f.write(l)
    except:
        f.write("please provide a dalton.head file\n")

    f.write("Charge=0 Atoms={} Basis=pointcharge\n".format(min(len(grid[igrid:]),maxbq)))
    nbq = 0
    for at in grid[igrid:]:
        f.write(
            "Bq     {0[0]:16.10f} {0[1]:16.10f} {0[2]:16.10f}\n".format(at))
        nbq = nbq + 1
        igrid = igrid + 1
        if (nbq == maxbq):
            f.close()
            logger.info("Batch generation : {}".format(igrid))
            generate_daltonFile(
                geom, grid, logger, outdir=outdir, igrid=igrid, maxbq = maxbq, all_atoms=all_atoms, atomtypes=  atomtypes)
            break
    f = open(daltonfile, "a")
    try:
        ffoot = open("dalton.foot", "r")
        for l in ffoot.readlines():
            f.write(l)
    except:
        f.write("please provide a dalton.foot file\n")
    f.close()
    return

