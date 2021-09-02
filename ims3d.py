import os
import sys
import argparse
import logging
import numpy as np
#import pymatgen

import geometry.geometry
import graph_theory.detect_cycle
import grids.angular
import grids.geode
import interface.gaussian

# Create logger
logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)

# create console handler and set level to error
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)

# create file handler and set level to info
fh = logging.FileHandler("log_ims_prep_angular", mode="w")
fh.setLevel(logging.INFO)

# create formatter
#formatter = logging.Formatter(
#    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
formatter = logging.Formatter(
    '%(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
logger.addHandler(fh)


def valtoRGB(values):
    """
    Returns RGB colors for each value of values
        arg: values[:]
    """
    min_val = np.min(values)
    max_val = np.max(values)
    rgb=[]
    for val in values:
        ratio = (val-min_val)/(max_val-min_val)
        if (ratio<0.5):
            R = 1
            B = 1 - 2 * ratio
            G = B
        else:
            B = 1
            R = 1 - ratio
            G = R
        rgb.append(np.asarray([R, G, B]))
    return rgb

def generate_command_line(args):
    command_line = "python3 ims3d.py "
    for arg in vars(args):
        command_line = command_line + " {} {}".format(arg, getattr(args, arg))
    return command_line

#def readgeom(f):
#    """ Store a geometry from a file into the geom list """
#    logger.debug("in readgeom")
#    fgeom = open(f, "r")
#    geom = []
#    for line in fgeom.readlines():
#        l = line.strip()
#        print(l)
#        geom.append(l)
#        logger.debug(l)
#    fgeom.close()
#    return geom

def main():

    #
    parser = argparse.ArgumentParser(
        description='Generate gaussian inputs for IMS calculations.')
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='More info')
    parser.add_argument(
        '-d',
        '--debug',
        action='store_true',
        help='Debug info')
    parser.add_argument(
        '-r',
        '--radius',
        type=float,
        help="Set the radius to 1 angstrom"
            )
    parser.add_argument(
        '-n',
        '--npts',
        type=int,
        help="Number of angular points by half circle. default: %(default)s",
        default=12)
    parser.add_argument(
        '--batch',
        '-b',
        type=int,
        help="Change the number of bq per batch. default: infinity",
        default=float('inf'))
    parser.add_argument(
        '--depth',
        type=int,
        help="Change the depth for geodesic grid generation: %(default)s",
        default=3)
    parser.add_argument(
        '-o',
        '--orient',
        action='store_true',
        help="Reorient the molecule along its principal symmetry axis",
        default=False)
    parser.add_argument(
        '-i',
        '--ignoreH',
        action='store_true',
        help="Ignore hydrogen atoms for the generation of the surface",
        default=False)
    parser.add_argument(
        '-p',
        '--preview',
        action='store_true',
        help="Preview the grid and the resulting surface",
        default=False)
    parser.add_argument(
        '-a',
        '--angular',
        action='store_true',
        help="Activate the deprecated angular grid",
        default=False)
    parser.add_argument(
        '-c',
        '--cycle-max-size',
        type=int,
        help='Auto detect cycles of max size: %(default)s',
        default=7)
    parser.add_argument(
        'geomfile',
        type=str,
        help="Geometry file in xyz format. default: %(default)s",
        default="geom.xyz")
    args = parser.parse_args()
    for arg in vars(args):
        print("{:} ... {:}".format(arg, getattr(args, arg)))
    if (args.debug):
        logger.setLevel(logging.DEBUG)
        fh.setLevel(logging.DEBUG)
    elif(args.verbose):
        logger.setLevel(logging.INFO)
    ignoreH = args.ignoreH
    preview = args.preview
    ntheta = args.npts
    orient = args.orient
    angular = args.angular
    depth = args.depth
    maxbq = args.batch
    cycle_max_size = args.cycle_max_size
    #
    # Read the geometry in the geom file
    #
    geomfile = args.geomfile
    geom = geometry.geometry.Geometry(geomfile, orient=orient)

    geomfile_atomsonly = geom.getgeomfilename_Atomsonly()
    cycles = []
    molecularGraph = graph_theory.detect_cycle.MolecularGraph(geomfile_atomsonly)
    for c in molecularGraph.getCycles():
        if len(c) <= cycle_max_size:
            cycles.append(list(c))

    os.remove(geomfile_atomsonly)
    if (len(cycles)>0):
        for cycle in cycles:
            atomlist = [int(str(i).replace('a', '')) - 0 for i in cycle]
            barycenter = geom.getBarycenter(atomlist)
            print(atomlist)
            print(barycenter)
            geom.addPseudoAtom(barycenter)

    #
    # Generate the full command_line
    #
    command_line = generate_command_line(args)
    print(command_line)
    logger.info(command_line)
    grid=[]
    if angular:
        if args.radius:
            radius_all = args.radius
            r_grid = grids.angular.angular_grid(ignoreH = ignoreH, ntheta = ntheta, radius_all = radius_all)
        else:
            r_grid = grids.angular.angular_grid(ignoreH = ignoreH, ntheta = ntheta, radius_all = None)
        angular_grid, angular_grid_normals = grids.angular.generate_angular_grid(geom, r_grid, logger)
        grids.angular.writegrid(angular_grid, angular_grid_normals)
        grid = angular_grid
    else:
        if args.radius:
            radius_all = args.radius
            geodesic_grid = grids.geode.geodesic_grid(ignoreH = ignoreH, depth = depth, radius_all = radius_all)
        else:
            geodesic_grid = grids.geode.geodesic_grid(ignoreH = ignoreH, depth = depth, radius_all = None)
        grid = grids.geode.generate_geodesic_grid(geom, geodesic_grid, logger)
        print(len(grid))
        grids.geode.writegrid(grid)
    interface.gaussian.generate_gaussianFile(geom, grid, logger, maxbq = maxbq)

    if preview==True:
        import open3d as o3d
        point_cloud = np.loadtxt("points_values.csv", delimiter=",", skiprows=1)
#        points_normals = np.loadtxt("normals.csv", delimiter=",", skiprows=1)
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(point_cloud[:,:3])
#        pcd.normals = o3d.utility.Vector3dVector(points_normals[:,:3])
#        point_rgb = valtoRGB(point_cloud[:,3])
#        pcd.colors = o3d.utility.Vector3dVector(np.asarray(point_rgb))
        o3d.visualization.draw_geometries([pcd])
        poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=9)[0]
        poisson_mesh.compute_vertex_normals()
        o3d.visualization.draw_geometries([poisson_mesh])
        o3d.io.write_triangle_mesh("./p_mesh_c.ply", poisson_mesh)

if __name__ == "__main__":
    main()
