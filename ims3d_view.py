import sys
import numpy as np
import logging
import colorsys
import matplotlib.pyplot as plt
import argparse
from matplotlib.colors import ListedColormap

import geometry.geometry
import graph_theory.detect_cycle
import math_utils.trigonometry

try :
    import pyvista as pv
except ModuleNotFoundError as error:
    pyvista = None
    print("pyvista module not found")

try:
    import open3d as o3d
except ModuleNotFoundError as error:
    open3d = None


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

class MyPlotter(pv.Plotter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_key_event("t", self.setTransparency)

    def setTransparency():
        return

def getIntermediateColor(val, min_val, max_val, color_min, color_max):
    delta = (max_val-min_val)
    ratio = (val-min_val)/delta
    color_val = [0,0,0]
    for i in range(len(color_val)):
        color_val[i] = color_min[i] + ratio * (color_max[i]-color_min[i])
    return color_val

def valtoRGB(values, color_min=[0, 0, 1], color_avg=[.5, .5, .5], color_max=[1, 1, 0],
        limits = [float("inf"), -10, -15, -20, -25, -30, float("-inf")],
        colors = [[1,1,1], [.75,.75,1], [.50,.50,1], [.25,.25,1], [.25,.00,1], [.50,.00,1]], colormode="auto"):
    """
    Returns RGB colors for each value of values
        arg: values[:]
    """
    min_val = np.min(values)
    max_val = np.max(values)
    avg_val = np.average(values)
    rgb=[]
    if colormode=="auto":
        for val in values:
            ratio = (val-min_val)/(max_val-min_val)
            if (ratio<0.5):
                color = getIntermediateColor(val,min_val, avg_val, color_min, color_avg)
            else:
                color = getIntermediateColor(val,avg_val, max_val, color_avg, color_max)
            rgb.append(np.asarray(color))
    elif colormode=="iso":
        for val in values:
            color=[0,0,0]
            for i in range(len(limits)-1):
                if val<limits[i] and val>=limits[i+1]:
                    color = colors[i]
            rgb.append(np.asarray(color))
    return rgb

def main():
    # Begin trating arguments
    parser = argparse.ArgumentParser(
        description='Display the calcuated data of IMS calculations.')
    parser.add_argument(
        '--showstat',
        action='store_true',
        help='Show statistics')
    parser.add_argument(
        '--open3d',
        action='store_false',
        help='Turn on open3d rendering')
    parser.add_argument(
        '--mate',
        action='store_true',
        help='Turn on mate rendering')
    parser.add_argument(
        '--colormode',
        '-c',
        type=str,
        default="auto",
        help='Color mode: auto|iso. default: %(default)s')
    parser.add_argument(
        '--representation',
        '-r',
        type=str,
        default="cloud",
        help='Representation mode: cloud|surface. default: %(default)s')
    parser.add_argument(
        '--shot',
        action='store_true',
        help='Save a picture')
    parser.add_argument(
        '--twopanels',
        '-2',
        action='store_true',
        help='Show the molecule')
    parser.add_argument(
        '--molecule',
        '-m',
        action='store_true',
        help='Show the molecule')
    args = parser.parse_args()
    showstat = args.showstat
    colormode = args.colormode
    mate = args.mate
    representation = args.representation
    molecule = args.molecule
    open3d = args.open3d
    twopanels = args.twopanels
    shot = args.shot
    # End trating arguments

    # Read ims.dat
    # The ims.dat format is:
    #  x, y, z, nx, ny, nz, val
    #skirow does not read the first line
    values =  np.loadtxt("ims.dat", delimiter=",", skiprows=1)

    if showstat:
        a = np.hstack(values[:,6])
        _ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
        plt.title("Repartition of IMS values")
        plt.show()

    geom = geometry.geometry.Geometry("geom.xyz")
    if open3d:
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(values[:,:3])
        pcd.normals = o3d.utility.Vector3dVector(values[:,3:6])
        spheres = []
        mesh_molecule=None
        #generate molecule in open3D format
        if molecule:
            for at in geom.atoms:
                mat = [[1, 0, 0, at['x']],
                       [0, 1, 0, at['y']],
                       [0, 0, 1, at['z']],
                       [0, 0, 0, 1     ] ]
                mesh_sphere = o3d.geometry.TriangleMesh.create_sphere(radius=0.4)
                mesh_sphere.transform(mat)
                mesh_sphere.compute_vertex_normals()
                if at['label'] == 'C':
                    color=[0.4, 0.4, 0.4]
                elif at['label'] == 'H':
                    color=[0.9, 0.9, 0.9]
                else:
                    color=[1.0, 0.0, 0.0]
                mesh_sphere.paint_uniform_color(color)
                spheres.append(mesh_sphere)
            cylinders = []
            molecularGraph = graph_theory.detect_cycle.MolecularGraph("geom.xyz")
            for e in molecularGraph.getEdges():
                idx1 = e.GetBeginAtomIdx()
                idx2 = e.GetEndAtomIdx()
                at1 = geom.getAtom(idx1)
                at2 = geom.getAtom(idx2)
                pos1 = np.asarray([at1['x'], at1['y'], at1['z']])
                pos2 = np.asarray([at2['x'], at2['y'], at2['z']])
                vect_bond = pos2 - pos1
                middle_bond = 0.5 * (pos1 + pos2)
                mat_translation = np.asarray([
                       [1, 0, 0, middle_bond[0]],
                       [0, 1, 0, middle_bond[1]],
                       [0, 0, 1, middle_bond[2]],
                       [0, 0, 0, 1]])
    #            print(mat_translation)
                vect_axis = np.cross(vect_bond, [0, 0, 1]) #cylinders are aligne along z when created
                theta = np.arcsin(np.linalg.norm(vect_axis)/np.linalg.norm(vect_bond))
                vect_axis = vect_axis / np.linalg.norm(vect_axis)
                ux = vect_axis[0]
                uy = vect_axis[1]
                uz = vect_axis[2]
                c = np.cos(theta)
                s = np.sin(theta)
                mat_rotation=np.asarray([
                        [ux * ux * (1-c) + c     , ux * uy * (1-c) - uz * s, ux * uz * (1-c) + uy * s, 0],
                        [ux * uy * (1-c) + uz * s, uy * uy * (1-c) + c     , uy * uz * (1-c) - ux * s, 0],
                        [ux * uz * (1-c) - uy * s, uy * uz  *(1-c) + ux * s, uz * uz * (1-c) + c     , 0],
                        [0                       , 0                       , 0                       , 1]
                        ]
                        )

                mesh_cylinder = o3d.geometry.TriangleMesh.create_cylinder(radius=.2, height=np.linalg.norm(vect_bond))
                mesh_cylinder.transform(mat_rotation)
                mesh_cylinder.transform(mat_translation)
                mesh_cylinder.compute_vertex_normals()
                mesh_cylinder.paint_uniform_color([0.0, 0.6, 1.0])
                cylinders.append(mesh_cylinder)

            mesh_molecule = spheres
            mesh_molecule.extend(cylinders)
        # end generate molecule in open3D format

        point_rgb = valtoRGB(values[:,6], colormode=colormode)
        pcd.colors = o3d.utility.Vector3dVector(np.asarray(point_rgb))
        if representation=='cloud':
            if molecule:
                o3d.visualization.draw_geometries([pcd] + mesh_molecule)
            else:
                o3d.visualization.draw_geometries([pcd])
        poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=9)[0]
        if not(mate):
            poisson_mesh.compute_vertex_normals()
        density_mesh = o3d.geometry.TriangleMesh()

        if representation=='surface':
            if molecule:
                vis = o3d.visualization.draw_geometries([poisson_mesh] + mesh_molecule)
            else:
                vis = o3d.visualization.draw_geometries([poisson_mesh])
        o3d.io.write_triangle_mesh("./surface_mesh.ply", poisson_mesh)
    else:
        #begin pysta
        points = values[:,:3]
        data = values[:,6]
        points = pv.pyvista_ndarray(points)
        datac = pv.pyvista_ndarray(data)

        point_cloud = pv.PolyData(points)
        point_cloud["IMS"] = datac
        cloud = pv.wrap(point_cloud)
        cloud.save('test.vtk', binary=False)
#        cloud.save('test_asc.vtk', binary=False)
#        cloud.save('test_bin.vtk', binary=True)

        alpha = 1
        superred = np.array([1, 1, 0, alpha])
        red = np.array([1, 0, 0, alpha])
        darkred = np.array([1, .25, .25, alpha])
        lightred = np.array([1, .75, .75, alpha])
        lightblue = np.array([1, .75, .75, alpha])
        darkblue = np.array([1, .25, .25, alpha])
        blue = np.array([1, 0, 0, alpha])
        purple = np.array([1, 0, 1, alpha])

        mapping = np.linspace(datac.min(), datac.max(), 256)
        newcolors = np.empty((256, 4))
        newcolors[mapping < 40]  = superred
        newcolors[mapping < 30]  = red
        newcolors[mapping < 20]   = darkred
        newcolors[mapping < 10]   = lightred
        newcolors[mapping < -10]  = lightblue
        newcolors[mapping < -20]  = darkblue
        newcolors[mapping < -30]  = blue
        newcolors[mapping <= -40] = purple
        my_colormap = ListedColormap(newcolors)

        #Palette (1) ACIE Karadakov:
        p1_level_1 = np.array([1, 0.9, 0.8, alpha]) # Pale Brown
        p1_level_2 = np.array([0.6, 0.7, 0.95, alpha]) # Clear blue
        p1_level_3 = np.array([0, 0.4, 0.85, alpha]) # Clear blue but no so clear
        p1_level_4 = np.array([0, 0, 1, alpha]) #Plain blue
        p1_level_5 = np.array([0, 0, 0.4, alpha]) # Noir blue
        #Palette (2) Yellow to blue:
        p2_level_1 = np.array([1, 0.5, 0, alpha]) # Orange
        p2_level_2 = np.array([1, 1, 0, alpha]) # Yellow
        p2_level_3 = np.array([0, 1, 0, alpha]) # Green
        p2_level_4 = np.array([0, 0, 1, alpha]) # Blue
        p2_level_5 = np.array([1, 0, 1, alpha]) # "Noir blue
        #Palette (3) The Oranges:
        p3_level_1 = np.array([1, 0.9, 0.75, alpha]) # White Orange
        p3_level_2 = np.array([1, 0.75, 0.5, alpha]) # Clear Orange
        p3_level_3 = np.array([1, 0.6, 0.25, alpha]) # Pale Orange
        p3_level_4 = np.array([1, 0.5, 0, alpha]) # Orange
        p3_level_5 = np.array([0.5, 0.25, 0, alpha]) # Dark Orange
        #Palette (4) The Reds:
        p4_level_1 = np.array([1, 0.8, 0.8, alpha]) # White Red
        p4_level_2 = np.array([1, 0.6, 0.6, alpha]) # Clear Red
        p4_level_3 = np.array([1, 0.2, 0.2, alpha]) # Pale Red
        p4_level_4 = np.array([1, 0, 0, alpha]) # Red
        p4_level_5 = np.array([0.5, 0, 0, alpha]) # Dark Red
        #Palette (5) The Greens:
        p5_level_1 = np.array([0.8, 1, 0.8, alpha]) # White Green
        p5_level_2 = np.array([0.6, 1, 0.6, alpha]) # Clear Green
        p5_level_3 = np.array([0.4, 1, 0.4, alpha]) # Pale Green
        p5_level_4 = np.array([0, 1, 0, alpha]) # Green
        p5_level_5 = np.array([0, 0.25, 0, alpha]) # Dark Green
        #Palette (6) The Purples:
        p6_level_1 = np.array([1, 0.8, 1, alpha]) # White Purple
        p6_level_2 = np.array([1, 0.6, 1, alpha]) # Clear Purple
        p6_level_3 = np.array([1, 0.4, 1, alpha]) # Pale Purple
        p6_level_4 = np.array([1, 0, 1, alpha]) # Purple
        p6_level_5 = np.array([0.25, 0, 0.25, alpha]) # Dark Purple
        #Palette (7) The Blues:
        p7_level_1 = np.array([0.8, 0.8, 1, alpha]) # White Blue
        p7_level_2 = np.array([0.6, 0.6, 1, alpha]) # Clear Blue
        p7_level_3 = np.array([0.4, 0.4, 1, alpha]) # Pale Blue
        p7_level_4 = np.array([0, 0, 1, alpha]) # Blue
        p7_level_5 = np.array([0, 0, 0.25, alpha]) # Dark Blue
        #Palette (8) The Greyscale:
        p8_level_1 = np.array([1, 1, 1, alpha]) # White
        p8_level_2 = np.array([0.8, 0.8, 0.8, alpha]) # Clear Grey
        p8_level_3 = np.array([0.6, 0.6, 0.6, alpha]) # Pale Grey
        p8_level_4 = np.array([0.4, 0.4, 0.4, alpha]) # Grey
        p8_level_5 = np.array([0, 0, 0, alpha]) # Black
        #Colors used
        mapping = np.linspace(datac.min(), datac.max(), 256)
        newcolors = np.empty((256, 4))
        newcolors[mapping >16.5 ]  = p4_level_5
        newcolors[mapping <16.5 ]  = p4_level_3
        newcolors[mapping <11 ]  = p4_level_2
        newcolors[mapping <5.5 ]  = p1_level_1
        newcolors[mapping <-5.5 ] = p1_level_2
        newcolors[mapping <-11 ] = p1_level_3
        newcolors[mapping <-16.5 ] = p1_level_5
        my_colormap = ListedColormap(newcolors)

        if twopanels:
            p = MyPlotter(shape=(1,2))
        else:
            p = MyPlotter()
        p.subplot(0, 0)
        p.add_points(point_cloud, render_points_as_spheres=True, cmap=my_colormap)
        if twopanels:
            p.subplot(0, 1)
        else:
            p.subplot(0, 0)
#
        spheres = []
        cylinders = []
        if molecule:
            for at in geom.atoms:
                mesh_sphere = pv.Sphere(radius=0.5, center=[at['x'], at['y'], at['z']])
                if at['label'] == 'C':
                    color=[0.4, 0.4, 0.4]
                elif at['label'] == 'H':
                    color=[0.9, 0.9, 0.9]
                else:
                    color=[1.0, 0.0, 0.0]
                spheres.append(mesh_sphere)
            molecularGraph = graph_theory.detect_cycle.MolecularGraph("geom.xyz")
            for e in molecularGraph.getEdges():
                idx1 = int(lbl1.replace('a',''))-1
                idx2 = int(lbl2.replace('a',''))-1
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
        p.link_views()
        if shot:
            p.show(screenshot='airplane.png')
        else:
            p.show()

if __name__ == "__main__":
    main()
