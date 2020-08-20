#!/usr/bin/python3

import sys
import numpy as np
from scipy.spatial.distance import cdist
import re
import os
import glob
import argparse
import logging
import jsonUtils
import jmol_interface
from constants import Bohr2Angstrom

# Create logger
logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)

# create console handler and set level to error
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)

# create file handler and set level to info
fh = logging.FileHandler("log_nics_harv")
fh.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


def generate_cubefile(geom, grid, grid_values, dx, dy, dz, nptx, npty, nptz):
    """
    Generate a cube file for the given geometry and grid
    """
    cubefile = "nics.cube"
    fio = open(cubefile, "w+")
    nat = len(geom)
    fio.write("head 1\n".format())
    fio.write("head 2\n".format())
    fio.write("{0:5d} {1[0]:12.6f} {1[1]:12.6f} {1[2]:12.6f}\n".format(
        nat, grid[0] / Bohr2Angstrom))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        nptx, dx / Bohr2Angstrom, 0, 0))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        npty, 0, dy / Bohr2Angstrom, 0))
    fio.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(
        nptz, 0, 0, dz / Bohr2Angstrom))
    for el in geom:
        if (el['label'].strip().lower() == "c"):
            nuclear_charge = 6
        elif (el['label'].strip().lower() == "h"):
            nuclear_charge = 1
        charge = 0.0
        fio.write(
            "{:5d} {:12.6f} {:12.6f} {:12.6f}  {:12.6f}\n".format(
                nuclear_charge,
                charge,
                el['x'] / Bohr2Angstrom,
                el['y'] / Bohr2Angstrom,
                el['z'] / Bohr2Angstrom))
    idx = 0
    for x in range(0, nptx + 1):
        for y in range(0, npty + 1):
            for z in range(0, nptz + 1):
                k = str(idx)
                val = 0
                if (k in grid_values.keys()):
                    val = grid_values[k]
                fio.write("{:12.6f}".format(val))
                if (z % 6 == 5):
                    fio.write("\n".format())
                idx = idx + 1
        fio.write("\n".format())
    return cubefile


def closest_node(node, nodes):
    # from
    # https://codereview.stackexchange.com/questions/28207/finding-the-closest-point-to-a-list-of-points
    return nodes[cdist([node], nodes).argmin()]


def generate_values_on_grid(geom, nics_grid, npts):
    """
    Put the nics_values on a cubic grid based on the proximity of a measured point and a grid point
    """
    xmin = min(a['x'] for a in nics_grid)
    xmax = max(a['x'] for a in nics_grid)
    ymin = min(a['y'] for a in nics_grid)
    ymax = max(a['y'] for a in nics_grid)
    zmin = min(a['z'] for a in nics_grid)
    zmax = max(a['z'] for a in nics_grid)
    #
    # Add 5% margin to handle volume plotting
    #
    xmin = xmin - xmin * .05
    ymin = ymin - ymin * .05
    zmin = zmin - zmin * .05
    xmax = xmax + xmax * .05
    ymax = ymax + ymax * .05
    zmax = zmax + zmax * .05
    logger.info(
        "xmin={} xmax={} ymin={} ymax={} zmin={} zmax={}".format(
            xmin,
            xmax,
            ymin,
            ymax,
            zmin,
            zmax))
    nptx = npty = nptz = npts
    tmp_lst = []
    # Generate grid : there is probably more efficient
    logger.info("Generating grid")
    for x in np.linspace(xmin, xmax, nptx):         # see dz
        for y in np.linspace(ymin, ymax, npty):     # see dz
            # dz = ( zmax - zmin ) / ( npts - 1 )
            for z in np.linspace(zmin, zmax, nptz):
                tmp_lst.append([x, y, z])
    grid = np.array(tmp_lst)
    grid_values = {}
    # remember dz = ( zmax - zmin ) / ( npts - 1 )
    dx = (xmax - xmin) / (nptx - 1)
    dy = (ymax - ymin) / (npty - 1)  # idem
    dz = (zmax - zmin) / (nptz - 1)  # idem
    for el in nics_grid:
        el_prox = closest_node([el['x'], el['y'], el['z']], grid)
        idz = (el_prox[2] - zmin) / dz
        idy = (el_prox[1] - ymin) / dy
        idx = (el_prox[0] - xmin) / dx
        index = idx * npty * nptz + idy * nptz + idz
        grid_values[str(int(index))] = el['nics']
        # store the value as a pair index (in str format), nics value: it is much better than storing a large number of
        # zero nics values on the grid
    return grid, grid_values, dx, dy, dz, nptx, npty, nptz


def readlogfile(logfile):
    """
    Read a guassian output file and store the geometry and the nics values (if any)
    """
    f = open(logfile, "r")
    store_geom = False
    store_nics = False
    index = 0
    geom = []
    nics_grid = []
    for l in f.readlines():
        if (("Charge" in l) and ("Multiplicity" in l)):
            store_geom = True
        if (store_geom and len(l) ==
                2):  # line with 1 space character and a carriage return symbol
            # end of geometry
            store_geom = False
        if (store_geom and not("Charge" in l)):
            atmp = l.split()
            geom.append({'label': str(atmp[0]),
                         'x': float(atmp[1]),
                         'y': float(atmp[2]),
                         'z': float(atmp[3])
                         })
        if ("Anisotropy" in l):
            atmp = l.split()
            geom[index]['nics'] = -float(atmp[7])
            index = index + 1
    # split data into two sparate lists
    # as one will process many log files and want only 1 geometry but the full
    # nics grid
    g = geom
    geom = []
    for el in g:
        if "Bq" in el['label']:  # it is a bq atom -> nics grid
            nics_grid.append(el)
        else:
            geom.append(el)
    return geom, nics_grid


def store_data(geom, nics_grid):
    geom_file = "geom.xyz"
    fio = open(geom_file, "w")
    fio.write("{}\n\n".format(len(geom)))
    for el in geom:
        fio.write(
            "{0[label]:s} {0[x]:16.10f}  {0[y]:16.10f}  {0[z]:16.10f}\n".format(el))
    fio.close()
    nics_file = "nics.dat"
    fio = open(nics_file, "w")
    for el in nics_grid:
        fio.write(
            "{0[x]:16.10f}  {0[y]:16.10f}  {0[z]:16.10f}  {0[nics]:16.10f}\n".format(el))
    fio.close()


def main():
    # 'application' code
#    logger.debug('debug message')
#    logger.info('info message')
#    logger.warning('warn message')
#    logger.error('error message')
#    logger.critical('critical message')
    #
    parser = argparse.ArgumentParser(
        description='Harvest the calcuated data of NICS calculations.')
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='More info')
    parser.add_argument(
        '--npts',
        '-n',
        type=int,
        default="10",
        help='Number of grid points in all 3 directions. default: %(default)s')
    parser.add_argument(
        'logfile',
        type=str,
        default="input_cycle_01_batch_01.log. default: %(default)s",
        help='Log filename of a series of calculations.')
    args = parser.parse_args()
    logfile = args.logfile
    npts = args.npts
    radical = re.sub(r'_cycle_\d*_batch_\d*.log$',
                     '', os.path.basename(logfile))
    dirname = os.path.dirname(logfile)

#
# Print for debugging
#
    logger.debug("logfile", logfile)
#
# Read the geometry stored in geom for all radical_###.log files
#  and the data for all these files
#
    logfiles = sorted(
        glob.glob(
            dirname +
            "/" +
            radical +
            "_cycle_[0-9]*_batch_[0-9]*.log"))
    geom = []
    nics_grid = []
    for f in logfiles:
        logger.info("Extracting from {0:s} ".format(f))
        geom_tmp, nics_grid_tmp = readlogfile(f)
        if len(geom) == 0:
            geom = geom_tmp
            logger.info("geometry and ".format())
        nics_grid.extend(nics_grid_tmp)
        logger.info("NICS values")
    store_data(geom, nics_grid)
    grid, grid_values, dx, dy, dz, nptx, npty, nptz = generate_values_on_grid(
        geom, nics_grid, npts)
    cubefile = generate_cubefile(
        geom,
        grid,
        grid_values,
        dx,
        dy,
        dz,
        nptx,
        npty,
        nptz)
    # Part 2 plot usin jmol
    planes = jsonUtils.load_state("state.json")["planes"]
    id_plane = 0
    for plane in planes:
        p = plane['plane']
        id_plane = id_plane + 1
        jmolfile = "commands_{}.jmol".format(id_plane)
        pngfile = "nics_{}.png".format(id_plane)
        jmol_interface.generate_jmolfile(
            jmolfile, cubefile, p, pngfile)


if __name__ == "__main__":
    main()
