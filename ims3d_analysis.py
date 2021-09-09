import sys
import numpy as np
import logging
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

# Create logger
logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)

# create console handler and set level to error
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)

# create file handler and set level to info
fh = logging.FileHandler("log_ims_analysis", mode="w")
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

    def setTransparency(self):
        return

def show(values, geom, graphicalElements):
    #begin pysta
    points = values[:,:3]
    data = values[:,6]
    points = pv.pyvista_ndarray(points)
    datac = pv.pyvista_ndarray(data)

    point_cloud = pv.PolyData(points)
    point_cloud["IMS"] = datac

    alpha = 0.5
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

    p = MyPlotter()
    p.subplot(0, 0)
    p.add_points(point_cloud, render_points_as_spheres=True, cmap=my_colormap)
#
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

            mesh_cylinder = pv.Cylinder(center=middle_bond, direction=vect_bond, radius=.2, height=np.linalg.norm(vect_bond))
            cylinders.append(mesh_cylinder)
#
        for sphere in spheres:
            p.add_mesh(sphere, color="tan", show_edges=False)
        for cyl in cylinders:
            p.add_mesh(cyl, color="tan", show_edges=False)
        for k in graphicalElements.keys():
            p.add_mesh(graphicalElements[k], color="yellow")
        p.link_views()
        p.show()

def main():
    # Begin treating arguments
    parser = argparse.ArgumentParser(
        description='Map analysis')
    parser.add_argument(
        '--showstat',
        action='store_true',
        help='Show statistics')
    args = parser.parse_args()
    showstat = args.showstat
    # End trating arguments

    # Read ims.dat
    # The ims.dat format is:
    #  x, y, z, nx, ny, nz, val
    #skirow does not read the first line
    values =  np.loadtxt("ims.dat", delimiter=",", skiprows=1)
    geom = geometry.geometry.Geometry("geom.xyz")

    if showstat:
        a = np.hstack(values[:,6])
        _ = plt.hist(a, bins='auto')  # arguments are passed to np.histogram
        plt.title("Repartition of IMS values")
        plt.show()

    graphicalElements = {}
    normal = geom.getNormalToThreeNearestNeighbougsPlane(0)
    normal_arrow = pv.Arrow(start=geom.getXYZ(0), direction=normal)

    graphicalElements["normal"]=normal_arrow
    show(values, geom, graphicalElements)

if __name__ == "__main__":
    main()
