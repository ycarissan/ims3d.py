#!/usr/bin/python3

import os
import sys
import numpy as np
import detect_cycle
import argparse
import logging
import jsonUtils
import geometry
import cubeUtils
import angularGrid
import gaussianUtils
import open3d as o3d
import colorsys

# Create logger
logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)

# create console handler and set level to error
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)

# create file handler and set level to info
fh = logging.FileHandler("log_nics_prep_angular", mode="w")
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
    command_line = "python3 nics_angular.py "
    for arg in vars(args):
        command_line = command_line + " {} {}".format(arg, getattr(args, arg))
    return command_line

def readgeom(f):
    """ Store a geometry from a file into the geom list """
    logger.debug("in readgeom")
    fgeom = open(f, "r")
    geom = []
    for line in fgeom.readlines():
        l = line.strip()
        print(l)
        geom.append(l)
        logger.debug(l)
    fgeom.close()
    return geom

def main():

    #
    parser = argparse.ArgumentParser(
        description='Generate gaussian inputs for NICS calculations.')
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
        '-i',
        '--ignoreH',
        action='store_true',
        help="Ignore hydrogen atoms for the generation of the surface",
        default=False)
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
    ntheta = args.npts
    #
    # Read the geometry in the geom file
    #
    geomfile = args.geomfile
    geom = geometry.Geometry(readgeom(geomfile))
    geomfile_atomsonly = geom.getgeomfilename_Atomsonly()
    cycles = detect_cycle.detect_cycles(geomfile_atomsonly)
    os.remove(geomfile_atomsonly)
    if (len(cycles)>0):
        for cycle in cycles:
            atomlist = [int(i.replace('a', '')) - 1 for i in cycle]
            barycenter = geom.getBarycenter(atomlist)
            print(atomlist)
            print(barycenter)
            geom.addPseudoAtom(barycenter)

    if args.radius:
        radius_all = args.radius
        r_grid = angularGrid.angular_grid(ignoreH = ignoreH, ntheta = ntheta, radius_all = radius_all)
    else:
        r_grid = angularGrid.angular_grid(ignoreH = ignoreH, ntheta = ntheta, radius_all = None)
    #
    # Generate the full command_line
    #
    command_line = generate_command_line(args)
    print(command_line)
    logger.info(command_line)
    angular_grid, angular_grid_normals = angularGrid.generate_angular_grid(geom, r_grid, logger)
    angularGrid.writegrid(angular_grid, angular_grid_normals)
    gaussianUtils.generate_gaussianFile(geom, angular_grid, logger)

    point_cloud = np.loadtxt("points_values.csv", delimiter=",", skiprows=1)
    points_normals = np.loadtxt("normals.csv", delimiter=",", skiprows=1)
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(point_cloud[:,:3])
    pcd.normals = o3d.utility.Vector3dVector(points_normals[:,:3])
#    pcd.colors = o3d.utility.Vector3dVector(point_cloud[:,3:6]/0.1)
    print(point_cloud[:,3:6])
    point_rgb = valtoRGB(point_cloud[:,3])
    pcd.colors = o3d.utility.Vector3dVector(np.asarray(point_rgb))
    o3d.visualization.draw_geometries([pcd])
#    o3d.visualization.draw_geometries([pcd],point_show_normal=True)
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=9)[0]
#    poisson_mesh.paint_uniform_color([1, 0.706, 0])
    poisson_mesh.compute_vertex_normals()
    density_mesh = o3d.geometry.TriangleMesh()
#    density_mesh.vertices = mesh.vertices
#    density_mesh.triangles = mesh.triangles
#    density_mesh.triangle_normals = mesh.triangle_normals
    o3d.visualization.draw_geometries([poisson_mesh])
    o3d.io.write_triangle_mesh("./p_mesh_c.ply", poisson_mesh)

if __name__ == "__main__":
    main()
