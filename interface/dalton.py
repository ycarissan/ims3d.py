import pymatgen

from interface.generic import split_geom_and_grid

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
        f.write("Charge=0 Atoms={} Basis=pointcharge\n".format(min(len(batch[ibatch]),maxbq)))
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

def readdalfile(logfile):
    """
    Read a dalton output file and store the geometry and the ims values (if any)
    """
    f = open(logfile, "r")
    store_geom = False
    store_ims = False
    index = 0
    geom = []
    for l in f.readlines():
        if store_geom:
            if countdown>0:
                countdown -= 1
            elif len(l)>1:
                atmp = l.split()
                if (len(atmp)==5):
                    geom.append({'label': str(atmp[0]),
                        'x': float(atmp[2])*.529177210,
                        'y': float(atmp[3])*.529177210,
                        'z': float(atmp[4])*.529177210
                        })
                else:
                    geom.append({'label': str(atmp[0]),
                        'x': float(atmp[1])*.529177210,
                        'y': float(atmp[2])*.529177210,
                        'z': float(atmp[3])*.529177210
                        })
            else:
                store_geom=False
        if ("Molecular geometry (au)" in l):
            store_geom = True
            countdown=2
        if ("@2" in l):
            atmp = l.split()
            if not "shielding" in l:
                if len(atmp)==10:
                    geom[index]['ims'] = float(atmp[3])
                    index = index + 1
                elif len(atmp)==9:
                    geom[index]['ims'] = float(atmp[2])
                    index = index + 1
    f.close()
    return split_geom_and_grid(geom)