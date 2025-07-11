import sys
import numpy as np
import re
import os
import glob
import argparse
import logging
import geometry.geometry

import grids.angular
from interface.dalton import readdalfile
from interface.gaussian import readlogfile
from interface.orca import readorcafile

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
        '-f',
        '--format',
        choices=['com', 'dal', 'orca'],
        help='output format: %(default)s',
        default="com")
    parser.add_argument(
        'logfile',
        type=str,
        default="input_cycle_01_batch_01.log",
        help='Log filename of a series of calculations. default: %(default)s')
    args = parser.parse_args()
    logfile = args.logfile
    npts = args.npts
    output_format = args.format
    geode = not(args.angular)
    logfiles= []
    #check if logfile has a name of the form input_cycle_i_batch_j.log where i and j are decimal numbers
    if output_format=="com":
        radical = re.sub(r'_cycle_\d*_batch_\d*.log$', '', os.path.basename(logfile))
        radical = re.sub(r'_batch_\d*.log$', '', os.path.basename(logfile))
        if not (re.match(r'input_cycle_\d+_batch_\d+\.log$', logfile) or re.match(r'input_batch_\d+\.log$', logfile)):
            logger.info("log file name {} is not of the form *_cycle_i_batch_j.log or *_batch_i.log with i and j decimal numbers.".format(logfile))
            logger.info("Proceeding without looking for batches i.e. : only {} is processed.".format(logfile))
            logfiles.append(logfile)
    else:
        radical = re.sub(r'_cycle_\d*_batch_\d*.out$', '', os.path.basename(logfile))
        radical = re.sub(r'_batch_\d*.out$', '', os.path.basename(logfile))
        if not (re.match(r'input_cycle_\d+_batch_\d+\.out$', logfile) or re.match(r'input_batch_\d+\.out$', logfile)):
            logger.info("out file name {} is not of the form *_cycle_i_batch_j.out or *_batch_i.out with i and j decimal numbers.".format(logfile))
            logger.info("Proceeding without looking for batches i.e. : only {} is processed.".format(logfile))
            logfiles.append(logfile)
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
    if output_format=="com":
        if len(logfiles)==0:
            logfiles = sorted(
                glob.glob(
                    dirname +
                    "/" +
                    radical +
                    "_batch_[0-9]*.log"))
    else:
        if len(logfiles)==0:
            logfiles = sorted(
                    glob.glob(
                        dirname +
                        "/" +
                        radical +
                        "_batch_[0-9]*.out"))
    logger.debug("dirname : {}\nradical : {}\n".format(dirname, radical))
    logger.debug("logfiles: {} ...".format(logfiles))
    geom = []
    ims_grid = []
    for f in logfiles:
        logger.info("Extracting from {0:s} ".format(f))
        if output_format=="com":
            geom_tmp, ims_grid_tmp = readlogfile(f)
        elif output_format=="dal":
            geom_tmp, ims_grid_tmp = readdalfile(f)
        elif output_format=="orca":
            geom_tmp, ims_grid_tmp = readorcafile(f)
        else:
            raise SystemExit("Unknown format: {}".format(output_format))
        if len(geom) == 0:
            geom = geom_tmp
            logger.info("geometry and ".format())
        ims_grid.extend(ims_grid_tmp)
        logger.info("IMS values")
    sym_ops = geometry.geometry.readSymmOps()
    if not(sym_ops==None):
        ims_grid = geometry.geometry.applySymmOps_onGrid(sym_ops, ims_grid)
    if geode:
        store_data(geom, ims_grid)
    else:
        grid, normals = grids.angular.readgrid()
        grids.angular.addNormals(ims_grid, grid, normals)
        store_data(geom, ims_grid, geode=False)

if __name__ == "__main__":
    main()
