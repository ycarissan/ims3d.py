from interface.generic import split_geom_and_grid


def generate_orcaFile(geom, grid, logger, outdir="./", igrid=0, maxbq=200):
    orcafile = outdir + \
        "input_batch_{:05d}.inp".format(igrid)
    f = open(orcafile, "w")
    f.write("! B3LYP 6-311++G(d,p) TightSCF NMR\n")
    f.write("%pal nprocs 8 end\n")
    f.write("%maxcore 2000\n")
    nat = 0
    # Si un jour la geometrie devait avoir une charge ou une multiplicité
    charge = geom.charge if hasattr(geom, 'charge') else 0
    multiplicity = geom.multiplicity if hasattr(geom, 'multiplicity') else 1
    f.write("* xyz {} {}\n".format(charge, multiplicity))
    for at in geom.atoms:
        f.write("{:4s} {:16.10f} {:16.10f} {:16.10f}\n".format(at['label'], at['x'], at['y'], at['z']))
        nat = nat + 1
    nbq = 0
    for at in grid[igrid:]:
        f.write("{:4s} {:16.10f} {:16.10f} {:16.10f} NewGTO S 1 1 1e6 1 end NewAuxJGTO S 1 1 2e6 1 end\n".format("H:", at[0], at[1], at[2]))
        nbq = nbq + 1
        nat = nat + 1
        igrid = igrid + 1
        if (nbq == maxbq):
            logger.info("Batch generation : {}".format(igrid))
            generate_orcaFile(
                geom, grid, logger, outdir=outdir, igrid=igrid, maxbq = maxbq)
            break
    f.write("*\n")
    f.close()
    return

def readorcafile(logfile):
    f = open(logfile, "r")
    store_geom = False
    store_ims = False
    index = 0
    geom = []
    nat = 0
    nbq = 0
    for l in f.readlines():
        if store_geom:
            if countdown>0:
                countdown -= 1
            elif len(l)>1:
                atmp = l.split()
#                print(atmp)
                if float(atmp[2]) == 0.0:
                    nbq = nbq + 1
                    lbl = "Bq"
                    if nat==0:
                        nat = len(geom)
                else:
                    lbl = str(atmp[1])
                geom.append({'label': lbl,
                             'x': float(atmp[5])*.529177210,
                             'y': float(atmp[6])*.529177210,
                             'z': float(atmp[7])*.529177210
                            })
            else:
                store_geom=False
        if ("CARTESIAN COORDINATES (A.U.)" in l):
            store_geom = True
            countdown=2
        if ("Total" in l and "iso=" in l):
            atmp = l.split()
            geom[index]['ims'] = float(atmp[5])
            index = index + 1
#    print(geom)
    f.close()
    return split_geom_and_grid(geom)


