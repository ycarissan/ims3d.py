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


