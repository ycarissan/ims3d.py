from interface.generic import split_geom_and_grid


def generate_pyscfFile(geom, grid, logger, outdir="./", igrid=0, maxbq=200):
    batch_start = igrid

    # Collect Bq positions for this batch
    bq_this_batch = []
    for at in grid[igrid:]:
        bq_this_batch.append(at)
        igrid += 1
        if len(bq_this_batch) == maxbq:
            break

    # Write Bq probe positions to a separate file (x y z, Angstrom, one per line)
    bq_filename = "input_batch_{:05d}_bq.txt".format(batch_start)
    with open(outdir + bq_filename, "w") as bf:
        for at in bq_this_batch:
            bf.write("{:16.10f}  {:16.10f}  {:16.10f}\n".format(at[0], at[1], at[2]))

    # Write PySCF script
    scriptfile = outdir + "input_batch_{:05d}.py".format(batch_start)
    unique_labels = sorted(set(at['label'] for at in geom.atoms))

    f = open(scriptfile, "w")
    f.write("from pyscf import gto, scf\n")
    f.write("import numpy as np\n")
    f.write("import copy, sys\n\n")

    # Real molecule (no Bq in basis)
    f.write("# ── Real molecule ─────────────────────────────────────────────────────\n")
    f.write("mol = gto.Mole()\n")
    f.write("mol.atom = '''\n")
    for at in geom.atoms:
        f.write("{:s}  {:16.10f}  {:16.10f}  {:16.10f}\n".format(
            at['label'], at['x'], at['y'], at['z']))
    f.write("'''\n")
    f.write("mol.basis = {\n")
    for lbl in unique_labels:
        f.write("    '{}': '6-311++g**',\n".format(lbl))
    f.write("}\n")
    f.write("mol.charge = 0\n")
    f.write("mol.spin = 0\n")
    f.write("mol.build()\n\n")

    # SCF on real molecule
    f.write("# ── DFT B3LYP ─────────────────────────────────────────────────────────\n")
    f.write("mf = scf.RKS(mol)\n")
    f.write("mf.xc = 'b3lyp'\n")
    f.write("mf.run()\n\n")

    # Read Bq probe positions from companion file
    f.write("# ── Bq probe positions ────────────────────────────────────────────────\n")
    f.write("bq_coords = []\n")
    f.write("with open('{}') as bq_file:\n".format(bq_filename))
    f.write("    for line in bq_file:\n")
    f.write("        line = line.strip()\n")
    f.write("        if line:\n")
    f.write("            x, y, z = map(float, line.split())\n")
    f.write("            bq_coords.append([x, y, z])\n\n")

    # Build mol_nmr = real atoms + X probes (X has no basis → same nao as mol)
    f.write("# ── Extended molecule: real atoms + X probes (no basis on X) ─────────\n")
    f.write("BOHR2ANG = 0.529177210\n")
    f.write("atom_list = []\n")
    f.write("for sym, coord_bohr in mol._atom:\n")
    f.write("    c = [v * BOHR2ANG for v in coord_bohr]\n")
    f.write("    atom_list.append('{} {} {} {}'.format(sym, c[0], c[1], c[2]))\n")
    f.write("for coord in bq_coords:\n")
    f.write("    atom_list.append('X {} {} {}'.format(coord[0], coord[1], coord[2]))\n\n")
    f.write("mol_nmr = gto.Mole()\n")
    f.write("mol_nmr.atom = '; '.join(atom_list)\n")
    f.write("mol_nmr.basis = mol._basis.copy()  # X not in basis dict → no basis functions\n")
    f.write("mol_nmr.charge = mol.charge\n")
    f.write("mol_nmr.spin = mol.spin\n")
    f.write("mol_nmr.build()\n\n")

    # NMR: copy mf, swap molecule (same nao → same mo_coeff dimensions)
    f.write("# ── NMR shielding (dia + para) at all positions including X probes ────\n")
    f.write("mf_nmr = copy.copy(mf)\n")
    f.write("mf_nmr.mol = mol_nmr\n")
    f.write("nmr_calc = mf_nmr.NMR()\n")
    f.write("shielding = nmr_calc.kernel()\n\n")

    # Output
    f.write("# ── Output ───────────────────────────────────────────────────────────\n")
    f.write("for i in range(mol_nmr.natm):\n")
    f.write("    lbl = mol_nmr.atom_symbol(i)\n")
    f.write("    coord = mol_nmr.atom_coord(i) * BOHR2ANG\n")
    f.write("    iso = np.trace(shielding[i]) / 3.0\n")
    f.write("    print('IMS3D_RESULT {:s} {:16.10f} {:16.10f} {:16.10f} {:16.10f}'.format(\n")
    f.write("        lbl, coord[0], coord[1], coord[2], iso))\n")
    f.write("sys.stdout.flush()\n")
    f.close()

    # Recurse for next batch
    if len(bq_this_batch) == maxbq and igrid < len(grid):
        if logger:
            logger.info("Batch generation : {}".format(igrid))
        generate_pyscfFile(geom, grid, logger, outdir=outdir, igrid=igrid, maxbq=maxbq)
    return


def readpyscffile(outfile):
    geom = []
    with open(outfile) as f:
        for l in f:
            if l.startswith("IMS3D_RESULT"):
                parts = l.split()
                lbl_raw = parts[1].lower()
                lbl = "Bq" if (lbl_raw == 'x' or lbl_raw.startswith('ghost')) else parts[1]
                geom.append({'label': lbl,
                             'x': float(parts[2]),
                             'y': float(parts[3]),
                             'z': float(parts[4]),
                             'ims': float(parts[5])})
    return split_geom_and_grid(geom)
