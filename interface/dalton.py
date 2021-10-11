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

    if all_atoms==None:
        all_atoms = geom.getUniqueElementsByLabels()
        print("generate all atoms")
    if atomtypes==None:
        atomtypes = geom.getUniqueLabels()
        print("generate all labels")

    for attype in atomtypes:
        if attype.capitalize()=="BQ":
            continue
        elif attype.capitalize()=="C":
            charge = 6.0
        elif attype.capitalize()=="O":
            charge = 8.0
        elif attype.capitalize()=="H":
            charge = 1.0
        nat = len(all_atoms[attype])
        f.write("Charge={} Atoms={} Basis=6-31g\n".format(charge, nat))
        for at in all_atoms[attype]:
            f.write("{:4s} {:16.10f} {:16.10f} {:16.10f}\n".format(at['label'], at['x'], at['y'], at['z']))

#    print(len(grid[igrid:]),maxbq)
    f.write("Charge=0 Atoms={} Basis=pointcharge\n".format(min(len(grid[igrid:]),maxbq)))
    nbq = 0
    for at in grid[igrid:]:
        f.write(
            "Bq     {0[0]:16.10f} {0[1]:16.10f} {0[2]:16.10f}\n".format(at))
        nbq = nbq + 1
        igrid = igrid + 1
        if (nbq == maxbq):
            logger.info("Batch generation : {}".format(igrid))
            generate_daltonFile(
                geom, grid, logger, outdir=outdir, igrid=igrid, maxbq = maxbq, all_atoms=all_atoms, atomtypes=  atomtypes)
            break
    try:
        ffoot = open("dalton.foot", "r")
        for l in ffoot.readlines():
            f.write(l)
    except:
        f.write("please provide a dalton.foot file\n")
    return

