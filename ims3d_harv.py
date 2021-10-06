import sys
import numpy as np
import re
import os
import glob
import argparse
import logging
import geometry.geometry

import grids.angular

# Create logger
logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)

# create console handler and set level to error
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)

# create file handler and set level to info
fh = logging.FileHandler("log_ims_harv")
fh.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)



def readlogfile(logfile):
    """
    Read a guassian output file and store the geometry and the ims values (if any)
    """
    f = open(logfile, "r")
    store_geom = False
    store_ims = False
    index = 0
    geom = []
    ims_grid = []
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
            geom[index]['ims'] = float(atmp[4])
            index = index + 1
    # split data into two sparate lists
    # as one will process many log files and want only 1 geometry but the full
    # ims grid
    g = geom
    geom = []
    for el in g:
        if "Bq" in el['label']:  # it is a bq atom -> ims grid
            ims_grid.append(el)
        else:
            geom.append(el)
    return geom, ims_grid


def store_data(geom, ims_grid, geode=True):
    geom_file = "geom.xyz"
    fio = open(geom_file, "w")
    fio.write("{}\n\n".format(len(geom)))
    for el in geom:
        fio.write(
            "{0[label]:s} {0[x]:16.10f}  {0[y]:16.10f}  {0[z]:16.10f}\n".format(el))
    fio.close()
    ims_file = "ims.dat"
    fio = open(ims_file, "w")
    for el in ims_grid:
        if geode:
            fio.write(
                "{0[x]:16.10f},  {0[y]:16.10f},  {0[z]:16.10f},  {0[x]:16.10f}, {0[y]:16.10f}, {0[z]:16.10f}, {0[ims]:16.10f}\n".format(el))
        else:
            fio.write(
                "{0[x]:16.10f},  {0[y]:16.10f},  {0[z]:16.10f},  {0[nx]:16.10f}, {0[ny]:16.10f}, {0[nz]:16.10f}, {0[ims]:16.10f}\n".format(el))
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
        description='Harvest the calcuated data of IMS calculations.')
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='More info')
    parser.add_argument(
        '--angular',
        '-a',
        action='store_true',
        help='Handle angular grid')
    parser.add_argument(
        '--npts',
        '-n',
        type=int,
        default="30",
        help='Number of grid points in all 3 directions. default: %(default)s')
    parser.add_argument(
        'logfile',
        type=str,
        default="input_cycle_01_batch_01.log",
        help='Log filename of a series of calculations. default: %(default)s')
    args = parser.parse_args()
    logfile = args.logfile
    npts = args.npts
    geode = not(args.angular)
    radical = re.sub(r'_cycle_\d*_batch_\d*.log$', '', os.path.basename(logfile))
    radical = re.sub(r'_batch_\d*.log$', '', os.path.basename(logfile))
    dirname = os.path.dirname(logfile)
    if len(dirname)==0:
        dirname="."

#
# Print for debugging
#
    logger.debug("logfile: {}".format(logfile))
#
# Read the geometry stored in geom for all radical_###.log files
#  and the data for all these files
#
    logfiles = sorted(
        glob.glob(
            dirname +
            "/" +
            radical +
            "_batch_[0-9]*.log"))
#            "_cycle_[0-9]*_batch_[0-9]*.log"))
    logger.debug("dirname : {}\nradical : {}\n".format(dirname, radical))
    logger.debug("logfiles: {} ...".format(logfiles))
    geom = []
    ims_grid = []
    for f in logfiles:
        logger.info("Extracting from {0:s} ".format(f))
        geom_tmp, ims_grid_tmp = readlogfile(f)
        if len(geom) == 0:
            geom = geom_tmp
            logger.info("geometry and ".format())
        ims_grid.extend(ims_grid_tmp)
        logger.info("IMS values")
    sym_ops = geometry.geometry.readSymmOps()
    ims_grid = geometry.geometry.applySymmOps_onGrid(sym_ops, ims_grid)
    if geode:
        store_data(geom, ims_grid)
    else:
        grid, normals = grids.angular.readgrid()
        grids.angular.addNormals(ims_grid, grid, normals)
        store_data(geom, ims_grid, geode=False)

if __name__ == "__main__":
    main()
