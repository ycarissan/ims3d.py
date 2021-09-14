import os
import sys
import argparse
import logging
import numpy as np

import geometry.geometry
import graph_theory.detect_cycle
import grids.angular
import grids.geode
import interface.gaussian

try :
    import pyvista as pv
    pyvista = True
except ModuleNotFoundError as error:
    pyvista = None
    print("pyvista module not found: preview not possible")

try :
    from ims3d_view import MyPlotter
    myplotter = True
except ModuleNotFoundError as error:
    print("ims3d_view not found in the current directory: preview not possible")
    myplotter = False

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
        description='Generate gaussian inputs for IMS calculations.',
        epilog="Make sure you use the python3 interpreter that comes with the conda environment. If this sentence makes no sense and you get error messages for missing librairies, please read the documentation.")
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
    if preview and (not pyvista or not myplotter):
        print("Preview impossible : overriding preview to False")
        preview = False
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
        values =  np.loadtxt("points_values.csv", delimiter=",", skiprows=1)
        points = values[:,:3]
        points = pv.pyvista_ndarray(points)
        point_cloud = pv.PolyData(points)
        cloud = pv.wrap(point_cloud)
        p = MyPlotter()
        p.subplot(0, 0)
        p.add_points(point_cloud, render_points_as_spheres=True)
        spheres = []
        cylinders = []
        if True:
            for at in geom.atoms:
                mesh_sphere = pv.Sphere(radius=0.5, center=[at['x'], at['y'], at['z']])
                if at['label'] == 'C':
                    color=[0.4, 0.4, 0.4]
                elif at['label'] == 'H':
                    color=[0.9, 0.9, 0.9]
                else:
                    color=[1.0, 0.0, 0.0]
                spheres.append(mesh_sphere)
            molecularGraph = graph_theory.detect_cycle.MolecularGraph(geomfile)
            for e in molecularGraph.getEdges():
                idx1 = e.GetBeginAtomIdx()
                idx2 = e.GetEndAtomIdx()
                at1 = geom.getAtom(idx1)
                at2 = geom.getAtom(idx2)
                pos1 = np.asarray([at1['x'], at1['y'], at1['z']])
                pos2 = np.asarray([at2['x'], at2['y'], at2['z']])
                vect_bond = pos2 - pos1
                middle_bond = 0.5 * (pos1 + pos2)

                mesh_cylinder = pv.Cylinder(center=middle_bond, direction=vect_bond, radius=.2, height=np.linalg.norm(vect_bond))
                cylinders.append(mesh_cylinder)
#
        for sphere in spheres:
            p.add_mesh(sphere, color="tan", show_edges=False)
        for cyl in cylinders:
            p.add_mesh(cyl, color="tan", show_edges=False)
        p.show()

if __name__ == "__main__":
    main()
